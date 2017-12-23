#!/usr/bin/env python2.7
# Copyright 2016 Arjun Arkal Rao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : protect/ProTECT.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import print_function
from collections import defaultdict
from multiprocessing import cpu_count

from protect.addons.assess_car_t_validity import run_car_t_validity_assessment
from protect.addons.assess_immunotherapy_resistance import run_itx_resistance_assessment
from protect.addons.assess_mhc_pathway import run_mhc_gene_assessment
from protect.alignment.dna import align_dna
from protect.alignment.rna import align_rna
from protect.binding_prediction.common import merge_mhc_peptide_calls, spawn_antigen_predictors
from protect.common import (delete_bams,
                            delete_fastqs,
                            email_report,
                            get_file_from_gdc,
                            get_file_from_s3,
                            get_file_from_url,
                            ParameterError,
                            parse_chromosome_string,
                            untargz,
                            is_gzipfile,
                            gunzip)
from protect.expression_profiling.rsem import wrap_rsem
from protect.haplotyping.phlat import merge_phlat_calls, phlat_disk, run_phlat
from protect.mutation_annotation.snpeff import run_snpeff, snpeff_disk
from protect.mutation_calling.common import run_mutation_aggregator
from protect.mutation_calling.fusion import wrap_fusion
from protect.mutation_calling.indel import run_indel_caller
from protect.mutation_calling.muse import run_muse
from protect.mutation_calling.mutect import run_mutect
from protect.mutation_calling.radia import run_radia
from protect.mutation_calling.somaticsniper import run_somaticsniper
from protect.mutation_calling.strelka import run_strelka
from protect.mutation_translation import run_transgene, transgene_disk
from protect.qc.rna import cutadapt_disk, run_cutadapt
from protect.rankboost import wrap_rankboost
from toil.job import Job, PromisedRequirement

import argparse
import os
import pkg_resources
import re
import shutil
import sys
import yaml


def _ensure_set_contains(test_object, required_object, test_set_name=None):
    """
    Ensure that the required entries (set or keys of a dict) are present in the test set or keys
    of the test dict.

    :param set|dict test_object: The test set or dict
    :param set|dict required_object: The entries that need to be present in the test set (keys of
            input dict if input is dict)
    :param str test_set_name: Optional name for the set
    :raises ParameterError: If required entry doesnt exist
    """
    assert isinstance(test_object, (set, dict)), '%s,%s' % (test_object, test_set_name)
    assert isinstance(required_object, (set, dict)), '%s,%s' % (required_object, test_set_name)
    # set(dict) = set of keys of the dict
    test_set = set(test_object)
    required_set = set(required_object)
    set_name = ' ' if test_set_name is None else ' entry "%s" of ' % test_set_name
    missing_opts = required_set.difference(test_set)
    if len(missing_opts) > 0:
        raise ParameterError('The following entries are missing in%sthe config file:\n%s'
                             % (set_name, '\n'.join(missing_opts)))


def _ensure_patient_group_is_ok(patient_object, patient_name=None):
    """
    Ensure that the provided entries for the patient groups is formatted properly.

    :param set|dict patient_object: The values passed to the samples patient group
    :param str patient_name: Optional name for the set
    :raises ParameterError: If required entry doesnt exist
    """
    from protect.addons.common import TCGAToGTEx
    assert isinstance(patient_object, (set, dict)), '%s,%s' % (patient_object, patient_name)
    # set(dict) = set of keys of the dict
    test_set = set(patient_object)
    if 'tumor_type' not in patient_object:
        raise ParameterError(('The patient entry for sample %s ' % patient_name) +
                             'does not contain a Tumor type.')
    elif patient_object['tumor_type'] not in TCGAToGTEx:
        raise ParameterError(('The patient entry for sample %s ' % patient_name) +
                             'does contains an invalid Tumor type. Please use one of the '
                             'valid TCGA tumor types.')
    if {'tumor_dna_fastq_1', 'normal_dna_fastq_1', 'tumor_rna_fastq_1'}.issubset(test_set):
        # Best case scenario, we get all fastqs
        pass
    else:
        # We have less than 3 fastqs so we have to have a haplotype.
        if 'hla_haplotype_files' not in test_set:
            raise ParameterError(('The patient entry for sample %s ' % patient_name) +
                                 'does not contain a hla_haplotype_files entry.\nCannot haplotype '
                                 'patient if all the input sequence files are not fastqs.')
        # Either we have a fastq and/or bam for the tumor and normal, or we need to be given a vcf
        if (({re.search('tumor_dna_((bam)|(fastq_1)).*', x) for x in test_set} == {None} or
                {re.search('normal_dna_((bam)|(fastq_1)).*', x) for x in test_set} == {None}) and
                'mutation_vcf' not in test_set):
            raise ParameterError(('The patient entry for sample %s ' % patient_name) +
                                 'does not contain a mutation_vcf entry. If\nboth tumor and normal '
                                 'DNA sequences (fastqs or bam) are not provided, a pre-computed '
                                 'vcf must\nbe provided.')
        # We have to be given a tumor rna fastq or bam
        if {re.search('tumor_rna_((bam)|(fastq_1)).*', x) for x in test_set} == {None}:
            raise ParameterError(('The patient entry for sample %s ' % patient_name) +
                                 'does not contain a tumor rna sequence data entry. We require '
                                 'either tumor_rna_fastq_1 or tumor_rna_bam.')
        # If we are given an RNA bam then it needs to have a corresponding transcriptome bam
        if 'tumor_rna_bam' in test_set and 'tumor_rna_transcriptome_bam' not in test_set:
            raise ParameterError(('The patient entry for sample %s ' % patient_name +
                                  'was provided a tumor rna bam with sequences mapped to the '
                                  'genome but was not provided a matching rna bam for the '
                                  'transcriptome. We require both.'))


def _add_default_entries(input_dict, defaults_dict):
    """
    Add the entries in defaults dict into input_dict if they don't exist in input_dict

    This is based on the accepted answer at
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    :param dict input_dict: The dict to be updated
    :param dict defaults_dict: Dict containing the defaults for entries in input_dict
    :return: updated dict
    :rtype: dict
    """
    for key, value in defaults_dict.iteritems():
        if key == 'patients':
            print('Cannot default `patients`.')
            continue
        if isinstance(value, dict):
            if key not in input_dict or input_dict[key] is None:
                # User didn't specify anython for the tool, but the entry was still in there so we
                # just copy over the whole defaults dict
                input_dict[key] = value
            else:
                r = _add_default_entries(input_dict.get(key, {}), value)
                input_dict[key] = r
        else:
            # Only write if not in input_dict
            if key not in input_dict or input_dict[key] is None:
                # Either the user didn't have the entry, or had it without a value
                input_dict[key] = value
    return input_dict


def _process_group(input_group, required_group, groupname, append_subgroups=None):
    """
    Process one group from the input yaml.  Ensure it has the required entries.  If there is a
    subgroup that should be processed and then appended to the rest of the subgroups in that group,
    handle it accordingly.

    :param dict input_group: The dict of values of the input group
    :param dict required_group: The dict of required values for the input group
    :param str groupname: The name of the group being processed
    :param list append_subgroups: list of subgroups to append to each, other subgroup in this group
    :return: processed dict of entries for the group
    :rtype: dict
    """
    if append_subgroups is None:
        append_subgroups = []
    tool_options = {}
    for key in input_group:
        _ensure_set_contains(input_group[key], required_group.get(key, {}), groupname + '::' + key)
        if key in append_subgroups:
            continue
        else:
            tool_options[key] = input_group[key]
    for key in input_group:
        if key in append_subgroups:
            continue
        else:
            for yek in append_subgroups:
                tool_options[key].update(input_group[yek])
    return tool_options


def get_fastq_2(job, patient_id, sample_type, fastq_1):
    """
    For a path to a fastq_1 file, return a fastq_2 file with the same prefix and naming scheme.

    :param str patient_id: The patient_id
    :param str sample_type: The sample type of the file
    :param str fastq_1: The path to the fastq_1 file
    :return: The path to the fastq_2 file
    :rtype: str
    """
    prefix, extn = fastq_1, 'temp'
    final_extn = ''
    while extn:
        prefix, extn = os.path.splitext(prefix)
        final_extn = extn + final_extn
        if prefix.endswith('1'):
            prefix = prefix[:-1]
            job.fileStore.logToMaster('"%s" prefix for "%s" determined to be %s'
                                      % (sample_type, patient_id, prefix))
            break
    else:
        raise ParameterError('Could not determine prefix from provided fastq (%s). Is it '
                             'of the form <fastq_prefix>1.[fq/fastq][.gz]?' % fastq_1)
    if final_extn not in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']:
        raise ParameterError('If and _2 fastq path is not specified, only .fastq, .fq or '
                             'their gzippped extensions are accepted. Could not process '
                             '%s:%s.' % (patient_id, sample_type + '_fastq_1'))
    return ''.join([prefix, '2', final_extn])


def parse_patients(job, patient_dict, skip_fusions=False):
    """
    Parse a dict of patient entries to retain only the useful entries (The user may provide more
    than we need and we don't want to download redundant things)

    :param dict patient_dict: The dict of patient entries parsed from the input config
    :param bool skip_fusions: A flag to identify if we're skipping fusions
    :return: A parsed dict of items
    :rtype: dict
    """
    output_dict = {'ssec_encrypted': patient_dict.get('ssec_encrypted') in (True, 'True', 'true'),
                   'patient_id': patient_dict['patient_id'],
                   'tumor_type': patient_dict['tumor_type'],
                   'filter_for_OxoG': patient_dict.get('filter_for_OxoG') in (True, 'True', 'true')}
    patient_keys = set(patient_dict)
    out_keys = []

    if {'mutation_vcf', 'hla_haplotype_files'}.issubset(patient_dict):
        # We have a haplotype and a vcf so don't need the DNA fastqs/bams
        out_keys.extend(['mutation_vcf', 'hla_haplotype_files'])
        # We just add the rna files to the dict here and resolve expression and fusions later
        if 'tumor_rna_bam' in patient_keys:
            out_keys.append('tumor_rna_bam')
            out_keys.append('tumor_rna_transcriptome_bam')
            if 'tumor_rna_bai' in patient_keys:
                out_keys.append('tumor_rna_bai')
        else:
            out_keys.extend([x for x in patient_keys if x.startswith('tumor_rna_fastq')])
    elif 'hla_haplotype_files' in patient_keys:
        out_keys.extend(['hla_haplotype_files'])
        # We don't have a vcf.  Use the bams if possible.
        for stype in 'tumor_dna', 'normal_dna', 'tumor_rna':
            if stype + '_bam' in patient_keys:
                out_keys.append(stype + '_bam')
                if stype + '_bai' in patient_keys:
                    out_keys.append(stype + '_bai')
                if stype == 'tumor_rna':
                    out_keys.append(stype + '_transcriptome_bam')
            else:
                out_keys.extend([x for x in patient_keys if x.startswith(stype + '_fastq')])
    else:
        if 'mutation_vcf' in patient_keys:
            # nesting this under the else since both cases need all fastqs.
            out_keys.extend(['mutation_vcf'])
        for stype in 'tumor_dna', 'normal_dna', 'tumor_rna':
            out_keys.extend([x for x in patient_keys if x.startswith(stype + '_fastq')])
    # If there's expression files, add them now
    if 'expression_files' in patient_keys:
        out_keys.append('expression_files')
    # If there's a fusion file, add it now
    if 'fusion_bedpe' in patient_keys:
        out_keys.append('fusion_bedpe')

    for key in out_keys:
        output_dict[key] = patient_dict[key]
    fastq1s = [x for x in output_dict if x.endswith('fastq_1')]
    for f in fastq1s:
        f = f[:-8]
        if f + '_fastq_2' not in output_dict:
            output_dict[f + '_fastq_2'] = get_fastq_2(job, patient_dict['patient_id'], f,
                                                      output_dict[f + '_fastq_1'])
    output_dict['gdc_inputs'] = [k for k, v in output_dict.items() if str(v).startswith('gdc')]
    return output_dict


def _parse_config_file(job, config_file, max_cores=None):
    """
    Parse the input yaml config file and download all tool inputs into the file store.

    :param str config_file: Path to the input config file
    :param int max_cores: The maximum cores to use for any single high-compute job.
    :return: tuple of dict of sample dicts, dict of universal options, dict of dicts of tool options
    :rtype: tuple(dict, dict, dict)
    """
    job.fileStore.logToMaster('Parsing config file')
    config_file = os.path.abspath(config_file)
    if not os.path.exists(config_file):
        raise ParameterError('The config file was not found at specified location. Please verify ' +
                             'and retry.')
    # Initialize variables to hold the sample sets, the universal options, and the per-tool options
    sample_set = {}
    univ_options = {}
    tool_options = {}

    # Read in the input yaml
    with open(config_file, 'r') as conf:
        input_config = yaml.load(conf.read())

    # Read the set of all required keys
    required_entries_file = pkg_resources.resource_filename(__name__, "required_entries.yaml")
    with open(required_entries_file) as ref:
        required_keys = yaml.load(ref.read())

    # Ensure every required entry is present in the file
    all_keys = set(input_config.keys())
    _ensure_set_contains(all_keys, required_keys)

    # Now that we are sure that mutation_callers is an entry in the file we can initialise this
    mutation_caller_list = []

    # Get the default values from the defaults file
    protect_defaults_file = pkg_resources.resource_filename(__name__, "defaults.yaml")
    with open(protect_defaults_file) as pdf:
        protect_defaults = yaml.load(pdf.read())

    # Update the input yaml with defaults
    input_config = _add_default_entries(input_config, protect_defaults)

    # Flags to check for presence of encryption keys if required
    gdc_inputs = ssec_encrypted = False
    for key in input_config.keys():
        if key == 'patients':
            # Ensure each patient contains the required entries
            for sample_name in input_config[key]:
                patient_keys = input_config[key][sample_name]
                _ensure_patient_group_is_ok(patient_keys, sample_name)
                # Add options to the sample_set dictionary
                input_config['patients'][sample_name]['patient_id'] = str(sample_name)
                sample_set[sample_name] = parse_patients(
                    job, input_config['patients'][sample_name],
                    not input_config['mutation_calling']['star_fusion']['run'])
                # If any of the samples requires gdc or ssec encrypted input, make the flag True
                if sample_set[sample_name]['ssec_encrypted']:
                    ssec_encrypted = True
                if sample_set[sample_name]['gdc_inputs']:
                    if 'tumor_rna_bam' in sample_set[sample_name]['gdc_inputs']:
                        raise ParameterError('Cannot run ProTECT using GDC RNA bams. Please fix '
                                             'sample %s' % sample_name)
                    gdc_inputs = True
        else:
            # Ensure the required entries exist for this key
            _ensure_set_contains(input_config[key], required_keys[key], key)
            if key == 'Universal_Options':
                univ_options.update(input_config['Universal_Options'])
                if univ_options['reference_build'].lower() in ['hg19', 'grch37']:
                    univ_options['ref'] = 'hg19'
                elif univ_options['reference_build'].lower() in ['hg38', 'grch38']:
                    univ_options['ref'] = 'hg38'
                else:
                    raise ParameterError('reference_build can only be hg19, hg38, GRCh37 or GRCh38')
                assert univ_options['storage_location'].startswith(('Local', 'local', 'aws'))
                if univ_options['storage_location'] in ('Local', 'local'):
                    assert os.path.isabs(univ_options['output_folder']), ('Needs to be absolute if '
                                                                          'storage_location is '
                                                                          'Local.')
                    assert univ_options['output_folder'] != 'NA', ('Cannot have NA as output '
                                                                   'folder if storage location is '
                                                                   'Local.')
                    univ_options['storage_location'] = 'local'
                if univ_options['storage_location'].startswith('aws'):
                    if 'sse_key' not in univ_options:
                        raise ParameterError('Cannot write results to aws without an sse key.')
                    if 'sse_key_is_master' in univ_options:
                        if univ_options['sse_key_is_master'] not in (True, False):
                            raise ParameterError('sse_key_is_master must be True or False')
                    else:
                        univ_options['sse_key_is_master'] = False
                univ_options['max_cores'] = cpu_count() if max_cores is None else max_cores
            else:
                if key == 'reports':
                    # The reporting group doesn't have any sub-dicts
                    tool_options[key] = input_config[key]
                else:
                    if key == 'alignment':
                        append_subgroup = ['post']
                    elif key == 'mutation_calling':
                        mutation_caller_list = input_config[key].keys()
                        append_subgroup = []
                    else:
                        append_subgroup = []
                    tool_options.update(_process_group(input_config[key], required_keys[key],
                                                       key, append_subgroup))
    # netmhciipan needs to be handled separately
    tool_options['mhcii']['netmhciipan'] = tool_options['netmhciipan']
    tool_options.pop('netmhciipan')
    # Check for encryption related issues before we download files.
    if ssec_encrypted:
        assert univ_options['sse_key'] is not None, 'Cannot read ssec encrypted data without a key.'
    if gdc_inputs:
        assert univ_options['gdc_download_token'] is not None, ('Cannot read from GDC without a '
                                                                'token.')
    # Get all the tool inputs
    job.fileStore.logToMaster('Obtaining tool inputs')
    process_tool_inputs = job.addChildJobFn(get_all_tool_inputs, tool_options,
                                            mutation_caller_list=mutation_caller_list)
    job.fileStore.logToMaster('Obtained tool inputs')
    return sample_set, univ_options, process_tool_inputs.rv()


def parse_config_file(job, config_file, max_cores=None):
    """
    Parse the config file and spawn a ProTECT job for every input sample.

    :param str config_file: Path to the input config file
    :param int max_cores: The maximum cores to use for any single high-compute job.
    """
    sample_set, univ_options, processed_tool_inputs = _parse_config_file(job, config_file,
                                                                         max_cores)
    # Start a job for each sample in the sample set
    for patient_id in sample_set.keys():
        job.addFollowOnJobFn(launch_protect, sample_set[patient_id], univ_options,
                             processed_tool_inputs)
    return None


def ascertain_cpu_share(max_cores=None):
    """
    Ascertain the number of cpus allowed for each high-compte job instance (bwa, star, rsem, phlat).

    :param max_cores: The user-specified max
    :return: The number of cpus allowed for each high-compte job instance
    :rtype: int
    """
    # Ascertain the number of available CPUs. Jobs will be given fractions of this value.
    num_cores = cpu_count()
    # The minimum number of cpus should be at least 6 if possible
    min_cores = min(num_cores, 6)
    cpu_share = max(num_cores / 2, min_cores)
    if max_cores is not None:
        cpu_share = min(cpu_share, max_cores)
    return cpu_share


def launch_protect(job, patient_data, univ_options, tool_options):
    """
    The launchpad for ProTECT. The DAG for ProTECT can be viewed in Flowchart.txt.

    :param dict patient_data: Dict of information regarding the input sequences for the patient
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict tool_options: Options for the various tools
    """
    # Add Patient id to univ_options as is is passed to every major node in the DAG and can be used
    # as a prefix for the logfile.
    univ_options['patient'] = patient_data['patient_id']
    univ_options['tumor_type'] = patient_data['tumor_type']
    # Ascertain number of cpus to use per job
    for tool in tool_options:
        tool_options[tool]['n'] = ascertain_cpu_share(univ_options['max_cores'])
    # Define the various nodes in the DAG
    # Need a logfile and a way to send it around
    sample_prep = job.wrapJobFn(prepare_samples, patient_data, univ_options, disk='40G')
    job.addChild(sample_prep)
    # Define the fastq deletion step
    fastq_deletion_1 = job.wrapJobFn(delete_fastqs, sample_prep.rv(), disk='100M', memory='100M')
    sample_prep.addChild(fastq_deletion_1)
    # Get all the input files
    haplotype_patient = get_mutations = None
    fastq_files = defaultdict(lambda: None)
    bam_files = defaultdict(lambda: None)
    delete_bam_files = defaultdict(lambda: None)
    phlat_files = defaultdict(lambda: None)
    for sample_type in 'tumor_dna', 'normal_dna', 'tumor_rna':
        if sample_type + '_fastq_1' in patient_data:
            fastq_files[sample_type] = job.wrapJobFn(get_patient_fastqs, sample_prep.rv(),
                                                     sample_type, disk='10M')
            sample_prep.addChild(fastq_files[sample_type])
            fastq_files[sample_type].addChild(fastq_deletion_1)
        elif sample_type + '_bam' in patient_data:
            bam_files[sample_type] = job.wrapJobFn(get_patient_bams, sample_prep.rv(), sample_type,
                                                   univ_options, tool_options['bwa'],
                                                   tool_options['mutect'],
                                                   disk='10M').encapsulate()
            sample_prep.addChild(bam_files[sample_type])

    # define the haplotyping subgraph of the DAG
    if 'hla_haplotype_files' in patient_data:
        haplotype_patient = job.wrapJobFn(get_patient_mhc_haplotype, sample_prep.rv())
        sample_prep.addChild(haplotype_patient)
    else:
        assert None not in fastq_files.values()
        # We are guaranteed to have fastqs here
        for sample_type in 'tumor_dna', 'normal_dna', 'tumor_rna':
            phlat_files[sample_type] = job.wrapJobFn(
                run_phlat, fastq_files[sample_type].rv(), sample_type, univ_options,
                tool_options['phlat'], cores=tool_options['phlat']['n'],
                disk=PromisedRequirement(phlat_disk, fastq_files[sample_type].rv()))
            fastq_files[sample_type].addChild(phlat_files[sample_type])
            phlat_files[sample_type].addChild(fastq_deletion_1)
        haplotype_patient = job.wrapJobFn(merge_phlat_calls,
                                          phlat_files['tumor_dna'].rv(),
                                          phlat_files['normal_dna'].rv(),
                                          phlat_files['tumor_rna'].rv(),
                                          univ_options, disk='100M', memory='100M', cores=1)
        phlat_files['tumor_dna'].addChild(haplotype_patient)
        phlat_files['normal_dna'].addChild(haplotype_patient)
        phlat_files['tumor_rna'].addChild(haplotype_patient)

    # Define the RNA-Seq Alignment subgraph if needed
    if bam_files['tumor_rna'] is None:
        assert fastq_files['tumor_rna'] is not None
        cutadapt = job.wrapJobFn(run_cutadapt, fastq_files['tumor_rna'].rv(), univ_options,
                                 tool_options['cutadapt'], cores=1,
                                 disk=PromisedRequirement(cutadapt_disk,
                                                          fastq_files['tumor_rna'].rv()))
        bam_files['tumor_rna'] = job.wrapJobFn(align_rna, cutadapt.rv(), univ_options,
                                               tool_options['star'], cores=1,
                                               disk='100M').encapsulate()
        fastq_deletion_2 = job.wrapJobFn(delete_fastqs, {'cutadapted_rnas': cutadapt.rv()},
                                         disk='100M', memory='100M')
        fastq_files['tumor_rna'].addChild(cutadapt)
        cutadapt.addChild(fastq_deletion_1)
        cutadapt.addChild(fastq_deletion_2)
        cutadapt.addChild(bam_files['tumor_rna'])
        bam_files['tumor_rna'].addChild(fastq_deletion_2)
        # Define the fusion calling node
        if 'fusion_bedpe' not in patient_data:
            tool_options['star_fusion']['index'] = tool_options['star']['index']
            tool_options['fusion_inspector']['index'] = tool_options['star']['index']
            fusions = job.wrapJobFn(wrap_fusion,
                                    cutadapt.rv(),
                                    bam_files['tumor_rna'].rv(),
                                    univ_options,
                                    tool_options['star_fusion'],
                                    tool_options['fusion_inspector'],
                                    disk='100M', memory='100M', cores=1).encapsulate()
        else:
            fusions = job.wrapJobFn(get_patient_bedpe, sample_prep.rv())
            sample_prep.addChild(fusions)
        bam_files['tumor_rna'].addChild(fusions)
        fusions.addChild(fastq_deletion_1)
        fusions.addChild(fastq_deletion_2)
    else:
        # Define the fusion calling node
        if 'fusion_bedpe' in patient_data:
            fusions = job.wrapJobFn(get_patient_bedpe, sample_prep.rv())
            sample_prep.addChild(fusions)
        else:
            if tool_options['star_fusion']['run'] is True:
                job.fileStore.logToMaster('Input RNA bams were provided for sample %s. Fusion '
                                          'detection can only be run with input '
                                          'fastqs.' % univ_options['patient'])
            fusions = None

    # Define the Expression estimation node
    if 'expression_files' in patient_data:
        rsem = job.wrapJobFn(get_patient_expression, sample_prep.rv())
        sample_prep.addChild(rsem)
    else:
        rsem = job.wrapJobFn(wrap_rsem, bam_files['tumor_rna'].rv(), univ_options,
                             tool_options['rsem'], cores=1, disk='100M').encapsulate()
        bam_files['tumor_rna'].addChild(rsem)
    # Define the bam deletion node
    delete_bam_files['tumor_rna'] = job.wrapJobFn(delete_bams,
                                                  bam_files['tumor_rna'].rv(),
                                                  univ_options['patient'], disk='100M',
                                                  memory='100M')
    bam_files['tumor_rna'].addChild(delete_bam_files['tumor_rna'])
    rsem.addChild(delete_bam_files['tumor_rna'])
    if fusions:
        fusions.addChild(delete_bam_files['tumor_rna'])
    # Define the reporting leaves
    if phlat_files['tumor_rna'] is not None:
        mhc_pathway_assessment = job.wrapJobFn(run_mhc_gene_assessment, rsem.rv(),
                                               phlat_files['tumor_rna'].rv(), univ_options,
                                               tool_options['reports'], disk='100M',
                                               memory='100M', cores=1)
        rsem.addChild(mhc_pathway_assessment)
        phlat_files['tumor_rna'].addChild(mhc_pathway_assessment)
    else:
        mhc_pathway_assessment = job.wrapJobFn(run_mhc_gene_assessment, rsem.rv(), None,
                                               univ_options, tool_options['reports'],
                                               disk='100M', memory='100M', cores=1)
        rsem.addChild(mhc_pathway_assessment)
    itx_resistance_assessment = job.wrapJobFn(run_itx_resistance_assessment, rsem.rv(),
                                              univ_options, tool_options['reports'],
                                              disk='100M', memory='100M', cores=1)
    rsem.addChild(itx_resistance_assessment)
    car_t_validity_assessment = job.wrapJobFn(run_car_t_validity_assessment, rsem.rv(),
                                              univ_options, tool_options['reports'],
                                              disk='100M', memory='100M', cores=1)
    rsem.addChild(car_t_validity_assessment)
    # Define the DNA-Seq alignment and mutation calling subgraphs if necessary
    if 'mutation_vcf' in patient_data:
        get_mutations = job.wrapJobFn(get_patient_vcf, sample_prep.rv())
        sample_prep.addChild(get_mutations)
    else:
        assert (None, None) not in zip(fastq_files.values(), bam_files.values())
        for sample_type in 'tumor_dna', 'normal_dna':
            if bam_files[sample_type] is None:
                assert fastq_files[sample_type] is not None
                bam_files[sample_type] = job.wrapJobFn(align_dna, fastq_files[sample_type].rv(),
                                                       sample_type, univ_options,
                                                       tool_options['bwa'], cores=1,
                                                       disk='100M').encapsulate()
                fastq_files[sample_type].addChild(bam_files[sample_type])
                bam_files[sample_type].addChild(fastq_deletion_1)
            else:
                # We already have the bam ready to go
                pass
            delete_bam_files[sample_type] = job.wrapJobFn(delete_bams,
                                                          bam_files[sample_type].rv(),
                                                          univ_options['patient'], disk='100M',
                                                          memory='100M')
            bam_files[sample_type].addChild(delete_bam_files[sample_type])
        # Time to call mutations
        mutations = {
            'radia': job.wrapJobFn(run_radia, bam_files['tumor_rna'].rv(),
                                   bam_files['tumor_dna'].rv(), bam_files['normal_dna'].rv(),
                                   univ_options, tool_options['radia'],
                                   disk='100M').encapsulate(),
            'mutect': job.wrapJobFn(run_mutect, bam_files['tumor_dna'].rv(),
                                    bam_files['normal_dna'].rv(), univ_options,
                                    tool_options['mutect'], disk='100M').encapsulate(),
            'muse': job.wrapJobFn(run_muse, bam_files['tumor_dna'].rv(),
                                  bam_files['normal_dna'].rv(), univ_options,
                                  tool_options['muse']).encapsulate(),
            'somaticsniper': job.wrapJobFn(run_somaticsniper, bam_files['tumor_dna'].rv(),
                                           bam_files['normal_dna'].rv(), univ_options,
                                           tool_options['somaticsniper']).encapsulate(),
            'strelka': job.wrapJobFn(run_strelka, bam_files['tumor_dna'].rv(),
                                     bam_files['normal_dna'].rv(), univ_options,
                                     tool_options['strelka']).encapsulate(),
            'indels': job.wrapJobFn(run_indel_caller, bam_files['tumor_dna'].rv(),
                                    bam_files['normal_dna'].rv(), univ_options, 'indel_options',
                                    disk='100M', memory='100M', cores=1)}
        for sample_type in 'tumor_dna', 'normal_dna':
            for caller in mutations:
                bam_files[sample_type].addChild(mutations[caller])
        bam_files['tumor_rna'].addChild(mutations['radia'])
        get_mutations = job.wrapJobFn(run_mutation_aggregator,
                                      {caller: cjob.rv() for caller, cjob in mutations.items()},
                                      univ_options, disk='100M', memory='100M',
                                      cores=1).encapsulate()
        for caller in mutations:
            mutations[caller].addChild(get_mutations)
        # We don't need the normal dna bam any more
        get_mutations.addChild(delete_bam_files['normal_dna'])
        # We may need the tumor one depending on OxoG
        if not patient_data['filter_for_OxoG']:
            get_mutations.addChild(delete_bam_files['tumor_dna'])

    # The rest of the subgraph should be unchanged
    snpeff = job.wrapJobFn(run_snpeff, get_mutations.rv(), univ_options, tool_options['snpeff'],
                           disk=PromisedRequirement(snpeff_disk,
                                                    tool_options['snpeff']['index']))
    get_mutations.addChild(snpeff)
    tumor_dna_bam = bam_files['tumor_dna'].rv() if patient_data['filter_for_OxoG'] else None
    fusion_calls = fusions.rv() if fusions else None
    transgene = job.wrapJobFn(run_transgene, snpeff.rv(), bam_files['tumor_rna'].rv(), univ_options,
                              tool_options['transgene'],
                              disk=PromisedRequirement(transgene_disk, bam_files['tumor_rna'].rv()),
                              memory='100M', cores=1, tumor_dna_bam=tumor_dna_bam,
                              fusion_calls=fusion_calls)
    snpeff.addChild(transgene)
    bam_files['tumor_rna'].addChild(transgene)
    transgene.addChild(delete_bam_files['tumor_rna'])
    if patient_data['filter_for_OxoG']:
        bam_files['tumor_dna'].addChild(transgene)
        transgene.addChild(delete_bam_files['tumor_dna'])
    if fusions:
        fusions.addChild(transgene)

    spawn_mhc = job.wrapJobFn(spawn_antigen_predictors, transgene.rv(), haplotype_patient.rv(),
                              univ_options, (tool_options['mhci'], tool_options['mhcii']),
                              disk='100M', memory='100M', cores=1).encapsulate()
    haplotype_patient.addChild(spawn_mhc)
    transgene.addChild(spawn_mhc)

    merge_mhc = job.wrapJobFn(merge_mhc_peptide_calls, spawn_mhc.rv(), transgene.rv(), univ_options,
                              disk='100M', memory='100M', cores=1)
    spawn_mhc.addFollowOn(merge_mhc)
    transgene.addChild(merge_mhc)

    rankboost = job.wrapJobFn(wrap_rankboost, rsem.rv(), merge_mhc.rv(), transgene.rv(),
                              univ_options, tool_options['rankboost'], disk='100M', memory='100M',
                              cores=1)
    rsem.addChild(rankboost)
    merge_mhc.addChild(rankboost)
    transgene.addChild(rankboost)
    report_success = job.wrapJobFn(email_report, univ_options)
    rankboost.addChild(report_success)
    return None


def get_all_tool_inputs(job, tools, outer_key='', mutation_caller_list=None):
    """
    Iterate through all the tool options and download required files from their remote locations.

    :param dict tools: A dict of dicts of all tools, and their options
    :param str outer_key: If this is being called recursively, what was the outer dict called?
    :param list mutation_caller_list: A list of mutation caller keys to append the indexes to.
    :return: The fully resolved tool dictionary
    :rtype: dict
    """
    for tool in tools:
        for option in tools[tool]:
            if isinstance(tools[tool][option], dict):
                tools[tool][option] = get_all_tool_inputs(
                    job, {option: tools[tool][option]},
                    outer_key=':'.join([outer_key, tool]).lstrip(':'))[option]
            else:
                # If a file is of the type file, vcf, tar or fasta, it needs to be downloaded from
                # S3 if reqd, then written to job store.
                if option.split('_')[-1] in ['file', 'vcf', 'index', 'fasta', 'fai', 'idx', 'dict',
                                             'tbi', 'beds', 'gtf', 'config']:
                    tools[tool][option] = job.addChildJobFn(
                        get_pipeline_inputs, ':'.join([outer_key, tool, option]).lstrip(':'),
                        tools[tool][option]).rv()
                elif option == 'version':
                    tools[tool][option] = str(tools[tool][option])
    if mutation_caller_list is not None:
        # Guaranteed to occur only in the outermost loop
        indexes = tools.pop('indexes')
        indexes['chromosomes'] = parse_chromosome_string(job, indexes['chromosomes'])
        for mutation_caller in mutation_caller_list:
            if mutation_caller == 'indexes':
                continue
            tools[mutation_caller].update(indexes)
    return tools


def get_pipeline_inputs(job, input_flag, input_file, encryption_key=None, per_file_encryption=False,
                        gdc_download_token=None):
    """
    Get the input file from s3 or disk and write to file store.

    :param str input_flag: The name of the flag
    :param str input_file: The value passed in the config file
    :param str encryption_key: Path to the encryption key if encrypted with sse-c
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :param str gdc_download_token: The download token to obtain files from the GDC
    :return: fsID for the file
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    job.fileStore.logToMaster('Obtaining file (%s) to the file job store' % input_flag)
    if input_file.startswith(('http', 'https', 'ftp')):
        input_file = get_file_from_url(job, input_file, encryption_key=encryption_key,
                                       per_file_encryption=per_file_encryption,
                                       write_to_jobstore=True)
    elif input_file.startswith(('S3', 's3')):
        input_file = get_file_from_s3(job, input_file, encryption_key=encryption_key,
                                      per_file_encryption=per_file_encryption,
                                      write_to_jobstore=True)
    elif input_file.startswith(('GDC', 'gdc')):
        input_file = get_file_from_gdc(job, input_file, gdc_download_token=gdc_download_token,
                                       write_to_jobstore=True)
    else:
        assert os.path.exists(input_file), 'Bogus Input : ' + input_file
        input_file = job.fileStore.writeGlobalFile(input_file)
    return input_file


def prepare_samples(job, patient_dict, univ_options):
    """
    Obtain the input files for the patient and write them to the file store.

    :param dict patient_dict: The input fastq dict
           patient_dict:
               |- 'tumor_dna_fastq_[12]' OR 'tumor_dna_bam': str
               |- 'tumor_rna_fastq_[12]' OR 'tumor_rna_bam': str
               |- 'normal_dna_fastq_[12]' OR 'normal_dna_bam': str
               |- 'mutation_vcf': str
               |- 'hla_haplotype_files': str
               +- 'patient_id': str
    :param dict univ_options: Dict of universal options used by almost all tools
    :return: Updated fastq dict
             output_dict:
                 |- 'tumor_dna_fastq_[12]' OR 'tumor_dna_bam': fsID
                 |- 'tumor_rna_fastq_[12]' OR 'tumor_rna_bam': fsID
                 |- 'normal_dna_fastq_[12]' OR 'normal_dna_bam': fsID
                 |- 'mutation_vcf': fsID
                 |- 'hla_haplotype_files': fsId
                 +- 'patient_id': str
    :rtype: dict
    """
    job.fileStore.logToMaster('Downloading Inputs for %s' % univ_options['patient'])
    # For each sample type, check if the prefix is an S3 link or a regular file
    # Download S3 files.
    output_dict = {}
    for input_file in patient_dict:
        if not input_file.endswith(('bam', 'bai', '_1', '_2', 'files', 'vcf', 'bedpe')):
            output_dict[input_file] = patient_dict[input_file]
            continue
        output_dict[input_file] = get_pipeline_inputs(
            job, ':'.join([univ_options['patient'], input_file]), patient_dict[input_file],
            encryption_key=(univ_options['sse_key'] if patient_dict['ssec_encrypted'] else None),
            per_file_encryption=univ_options['sse_key_is_master'],
            gdc_download_token=univ_options['gdc_download_token'])
    return output_dict


def get_patient_fastqs(job, patient_dict, sample_type):
    """
    Convenience function to return only a list of fq_1, fq_2 for a sample type.

    :param dict patient_dict: dict of patient info
    :param str sample_type: key in sample_type to return
    :return: fastqs[sample_type]
    :rtype: list[toil.fileStore.FileID]
    """
    return [patient_dict[sample_type + '_fastq_1'], patient_dict[sample_type + '_fastq_2']]


def get_patient_bams(job, patient_dict, sample_type, univ_options, bwa_options, mutect_options):
    """
    Convenience function to return the bam and its index in the correct format for a sample type.

    :param dict patient_dict: dict of patient info
    :param str sample_type: 'tumor_rna', 'tumor_dna', 'normal_dna'
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict bwa_options: Options specific to bwa
    :param dict bwa_options: Options specific to mutect
    :return: formatted dict of bam and bai
    :rtype: dict
    """
    output_dict = {}
    if 'dna' in sample_type:
        sample_info = 'fix_pg_sorted'
        prefix = sample_type + '_' + sample_info
    else:
        sample_info = 'genome_sorted'
        prefix = 'rna_' + sample_info
    if sample_type + '_bam' in patient_dict['gdc_inputs']:
        output_dict[prefix + '.bam'] = patient_dict[sample_type + '_bam'][0]
        output_dict[prefix + '.bam.bai'] = patient_dict[sample_type + '_bam'][1]
    elif sample_type + '_bai' in patient_dict:
        output_dict[prefix + '.bam'] = patient_dict[sample_type + '_bam']
        output_dict[prefix + '.bam.bai'] = patient_dict[sample_type + '_bai']
    else:
        from protect.alignment.dna import index_bamfile, index_disk
        output_job = job.wrapJobFn(index_bamfile, patient_dict[sample_type + '_bam'],
                                   'rna' if sample_type == 'tumor_rna' else sample_type,
                                   univ_options, bwa_options['samtools'],
                                   sample_info=sample_info, export=False,
                                   disk=PromisedRequirement(index_disk,
                                                            patient_dict[sample_type + '_bam']))
        job.addChild(output_job)
        output_dict = output_job.rv()
    if sample_type == 'tumor_rna':
        return{'rna_genome': output_dict,
               'rna_transcriptome.bam': patient_dict['tumor_rna_transcriptome_bam']}
    else:
        return output_dict


def get_patient_vcf(job, patient_dict):
    """
    Convenience function to get the vcf from the patient dict

    :param dict patient_dict: dict of patient info
    :return: The vcf
    :rtype: toil.fileStore.FileID
    """
    temp = job.fileStore.readGlobalFile(patient_dict['mutation_vcf'],
                                        os.path.join(os.getcwd(), 'temp.gz'))
    if is_gzipfile(temp):
        outfile = job.fileStore.writeGlobalFile(gunzip(temp))
        job.fileStore.deleteGlobalFile(patient_dict['mutation_vcf'])
    else:
        outfile = patient_dict['mutation_vcf']
    return outfile


def get_patient_bedpe(job, patient_dict):
    """
    Convenience function to get the bedpe from the patient dict

    :param dict patient_dict: dict of patient info
    :return: The bedpe
    :rtype: toil.fileStore.FileID
    """
    temp = job.fileStore.readGlobalFile(patient_dict['fusion_bedpe'],
                                        os.path.join(os.getcwd(), 'temp.gz'))
    if is_gzipfile(temp):
        outfile = job.fileStore.writeGlobalFile(gunzip(temp))
        job.fileStore.deleteGlobalFile(patient_dict['fusion_bedpe'])
    else:
        outfile = patient_dict['fusion_bedpe']
    return outfile


def get_patient_mhc_haplotype(job, patient_dict):
    """
    Convenience function to get the mhc haplotype from the patient dict

    :param dict patient_dict: dict of patient info
    :return: The MHCI and MHCII haplotypes
    :rtype: toil.fileStore.FileID
    """
    haplotype_archive = job.fileStore.readGlobalFile(patient_dict['hla_haplotype_files'])
    haplotype_archive = untargz(haplotype_archive, os.getcwd())
    output_dict = {}
    for filename in 'mhci_alleles.list', 'mhcii_alleles.list':
        output_dict[filename] = job.fileStore.writeGlobalFile(os.path.join(haplotype_archive,
                                                                           filename))
    return output_dict


def get_patient_expression(job, patient_dict):
    """
    Convenience function to get the expression from the patient dict

    :param dict patient_dict: dict of patient info
    :return: The gene and isoform expression
    :rtype: toil.fileStore.FileID
    """
    expression_archive = job.fileStore.readGlobalFile(patient_dict['expression_files'])
    expression_archive = untargz(expression_archive, os.getcwd())
    output_dict = {}
    for filename in 'rsem.genes.results', 'rsem.isoforms.results':
        output_dict[filename] = job.fileStore.writeGlobalFile(os.path.join(expression_archive,
                                                                           filename))
    return output_dict


def generate_config_file():
    """
    Generate a config file for a ProTECT run on hg19.

    :return: None
    """
    shutil.copy(os.path.join(os.path.dirname(__file__), 'input_parameters.yaml'),
                os.path.join(os.getcwd(), 'ProTECT_config.yaml'))


def main():
    """
    This is the main function for ProTECT.
    """
    parser = argparse.ArgumentParser(prog='ProTECT',
                                     description='Prediction of T-Cell Epitopes for Cancer Therapy',
                                     epilog='Contact Arjun Rao (aarao@ucsc.edu) if you encounter '
                                     'any problems while running ProTECT')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--config_file', dest='config_file', help='Config file to be used in the '
                        'run.', type=str, default=None)
    inputs.add_argument('--generate_config', dest='generate_config', help='Generate a config file '
                        'in the current directory that is pre-filled with references and flags for '
                        'an hg19 run.', action='store_true', default=False)
    parser.add_argument('--max-cores-per-job', dest='max_cores', help='Maximum cores to use per '
                        'job. Aligners and Haplotypers ask for cores dependent on the machine that '
                        'the launchpad gets assigned to -- In a heterogeneous cluster, this can '
                        'lead to problems. This value should be set to the number of cpus on the '
                        'smallest node in a cluster.',
                        type=int, required=False, default=None)
    # We parse the args once to see if the user has asked for a config file to be generated.  In
    # this case, we don't need a jobstore.  To handle the case where Toil arguments are passed to
    # ProTECT, we parse known args, and if the used specified config_file instead of generate_config
    # we re-parse the arguments with the added Toil parser.
    params, others = parser.parse_known_args()
    if params.generate_config:
        generate_config_file()
    else:
        Job.Runner.addToilOptions(parser)
        params = parser.parse_args()
        params.config_file = os.path.abspath(params.config_file)
        if params.maxCores:
            if not params.max_cores:
                params.max_cores = int(params.maxCores)
            else:
                if params.max_cores > int(params.maxCores):
                    print("The value provided to max-cores-per-job (%s) was greater than that "
                          "provided to maxCores (%s). Setting max-cores-per-job = maxCores." %
                          (params.max_cores, params.maxCores), file=sys.stderr)
                    params.max_cores = int(params.maxCores)
        start = Job.wrapJobFn(parse_config_file, params.config_file, params.max_cores)
        Job.Runner.startToil(start, params)
    return None


if __name__ == '__main__':
    main()
