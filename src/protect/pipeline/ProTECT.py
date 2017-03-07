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
from multiprocessing import cpu_count

from protect.addons import run_mhc_gene_assessment
from protect.alignment.dna import align_dna
from protect.alignment.rna import align_rna
from protect.binding_prediction.common import merge_mhc_peptide_calls, spawn_antigen_predictors
from protect.common import delete_fastqs, get_file_from_s3, get_file_from_url, ParameterError
from protect.expression_profiling.rsem import wrap_rsem
from protect.haplotyping.phlat import merge_phlat_calls, phlat_disk, run_phlat
from protect.mutation_annotation.snpeff import run_snpeff, snpeff_disk
from protect.mutation_calling.common import run_mutation_aggregator
from protect.mutation_calling.fusion import run_fusion_caller
from protect.mutation_calling.indel import run_indel_caller
from protect.mutation_calling.muse import run_muse
from protect.mutation_calling.mutect import run_mutect
from protect.mutation_calling.radia import run_radia
from protect.mutation_calling.somaticsniper import run_somaticsniper
from protect.mutation_calling.strelka import run_strelka
from protect.mutation_translation import run_transgene
from protect.qc.rna import cutadapt_disk, run_cutadapt
from protect.rankboost import wrap_rankboost
from toil.job import Job, PromisedRequirement

import argparse
import os
import pkg_resources
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

    # Read the set of all requried keys
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

    for key in input_config.keys():
        if key == 'patients':
            # Ensure each patient contains the required entries
            for sample_name in input_config[key]:
                patient_keys = input_config[key][sample_name]
                _ensure_set_contains(patient_keys, required_keys[key]['test'], sample_name)
                if 'ssec_encrypted' not in patient_keys:
                    input_config[key][sample_name]['ssec_encrypted'] = False
                else:
                    input_config[key][sample_name]['ssec_encrypted'] = \
                        bool(input_config[key][sample_name]['ssec_encrypted'])
            # Add options to the sample_set dictionary
            sample_set.update(input_config['patients'])
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
                if key == 'mhc_pathway_assessment':
                    # netmhciipan needs to be handled separately
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
    # Get all the tool inputs
    job.fileStore.logToMaster('Obtaining tool inputs')
    process_tool_inputs = job.addChildJobFn(get_all_tool_inputs, tool_options)
    for mutation_caller in mutation_caller_list:
        if mutation_caller == 'indexes':
            continue
        tool_options[mutation_caller].update(tool_options['indexes'])
    tool_options.pop('indexes')
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
        # Add the patient id to the sample set.  Typecast to str so cat operations work later on.
        sample_set[patient_id]['patient_id'] = str(patient_id)
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


def launch_protect(job, fastqs, univ_options, tool_options):
    """
    The launchpad for ProTECT. The DAG for ProTECT can be viewed in Flowchart.txt.

    :param dict fastqs: Dict of lists of fastq files and the patient ID
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict tool_options: Options for the various tools
    """
    # Add Patient id to univ_options as is is passed to every major node in the DAG and can be used
    # as a prefix for the logfile.
    univ_options['patient'] = fastqs['patient_id']
    # Ascertin number of cpus to use per job
    tool_options['star']['n'] = tool_options['bwa']['n'] = tool_options['phlat']['n'] = \
        tool_options['rsem']['n'] = ascertain_cpu_share(univ_options['max_cores'])
    # Define the various nodes in the DAG
    # Need a logfile and a way to send it around
    sample_prep = job.wrapJobFn(prepare_samples, fastqs, univ_options, disk='40G')
    tumor_dna_fqs = job.wrapJobFn(get_fqs, sample_prep.rv(), 'tumor_dna', disk='10M')
    normal_dna_fqs = job.wrapJobFn(get_fqs, sample_prep.rv(), 'normal_dna', disk='10M')
    tumor_rna_fqs = job.wrapJobFn(get_fqs, sample_prep.rv(), 'tumor_rna', disk='10M')
    cutadapt = job.wrapJobFn(run_cutadapt, tumor_rna_fqs.rv(), univ_options,
                             tool_options['cutadapt'], cores=1,
                             disk=PromisedRequirement(cutadapt_disk, tumor_rna_fqs.rv()
                                                      ))
    star = job.wrapJobFn(align_rna, cutadapt.rv(), univ_options, tool_options['star'],
                         cores=1, disk='100M').encapsulate()
    bwa_tumor = job.wrapJobFn(align_dna, tumor_dna_fqs.rv(), 'tumor_dna', univ_options,
                              tool_options['bwa'], cores=1, disk='100M'
                              ).encapsulate()
    bwa_normal = job.wrapJobFn(align_dna, normal_dna_fqs.rv(), 'normal_dna', univ_options,
                               tool_options['bwa'], cores=1, disk='100M'
                               ).encapsulate()
    phlat_tumor_dna = job.wrapJobFn(run_phlat, tumor_dna_fqs.rv(), 'tumor_dna', univ_options,
                                    tool_options['phlat'], cores=tool_options['phlat']['n'],
                                    disk=PromisedRequirement(phlat_disk, tumor_dna_fqs.rv()))
    phlat_normal_dna = job.wrapJobFn(run_phlat, normal_dna_fqs.rv(), 'normal_dna', univ_options,
                                     tool_options['phlat'], cores=tool_options['phlat']['n'],
                                     disk=PromisedRequirement(phlat_disk, normal_dna_fqs.rv()))
    phlat_tumor_rna = job.wrapJobFn(run_phlat, tumor_rna_fqs.rv(), 'tumor_rna', univ_options,
                                    tool_options['phlat'], cores=tool_options['phlat']['n'],
                                    disk=PromisedRequirement(phlat_disk, tumor_rna_fqs.rv()))
    fastq_deletion_1 = job.wrapJobFn(delete_fastqs, sample_prep.rv(), disk='100M', memory='100M')
    fastq_deletion_2 = job.wrapJobFn(delete_fastqs, {'cutadapted_rnas': cutadapt.rv()},
                                     disk='100M', memory='100M')
    rsem = job.wrapJobFn(wrap_rsem, star.rv(), univ_options, tool_options['rsem'], cores=1,
                         disk='100M').encapsulate()
    mhc_pathway_assessment = job.wrapJobFn(run_mhc_gene_assessment, rsem.rv(), phlat_tumor_rna.rv(),
                                           univ_options, tool_options['mhc_pathway_assessment'],
                                           disk='100M', memory='100M', cores=1)
    fusions = job.wrapJobFn(run_fusion_caller, star.rv(), univ_options, 'fusion_options',
                            disk='100M', memory='100M', cores=1)
    radia = job.wrapJobFn(run_radia, star.rv(), bwa_tumor.rv(),
                          bwa_normal.rv(), univ_options, tool_options['radia'],
                          disk='100M').encapsulate()
    mutect = job.wrapJobFn(run_mutect, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                           tool_options['mutect'], disk='100M').encapsulate()
    muse = job.wrapJobFn(run_muse, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                         tool_options['muse']).encapsulate()
    somaticsniper = job.wrapJobFn(run_somaticsniper, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                                  tool_options['somaticsniper']).encapsulate()
    strelka = job.wrapJobFn(run_strelka, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                            tool_options['strelka']).encapsulate()
    indels = job.wrapJobFn(run_indel_caller, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                           'indel_options', disk='100M', memory='100M', cores=1)
    merge_mutations = job.wrapJobFn(run_mutation_aggregator,
                                    {'fusions': fusions.rv(),
                                     'radia': radia.rv(),
                                     'mutect': mutect.rv(),
                                     'strelka': strelka.rv(),
                                     'indels': indels.rv(),
                                     'muse': muse.rv(),
                                     'somaticsniper': somaticsniper.rv()}, univ_options,
                                    disk='100M', memory='100M',
                                    cores=1).encapsulate()
    snpeff = job.wrapJobFn(run_snpeff, merge_mutations.rv(), univ_options, tool_options['snpeff'],
                           disk=PromisedRequirement(snpeff_disk,
                                                    tool_options['snpeff']['index']))
    transgene = job.wrapJobFn(run_transgene, snpeff.rv(), star.rv(), univ_options,
                              tool_options['transgene'], disk='100M', memory='100M', cores=1)
    merge_phlat = job.wrapJobFn(merge_phlat_calls, phlat_tumor_dna.rv(), phlat_normal_dna.rv(),
                                phlat_tumor_rna.rv(), univ_options, disk='100M', memory='100M',
                                cores=1)
    spawn_mhc = job.wrapJobFn(spawn_antigen_predictors, transgene.rv(), merge_phlat.rv(),
                              univ_options, (tool_options['mhci'], tool_options['mhcii']),
                              disk='100M', memory='100M', cores=1).encapsulate()
    merge_mhc = job.wrapJobFn(merge_mhc_peptide_calls, spawn_mhc.rv(), transgene.rv(), univ_options,
                              disk='100M', memory='100M', cores=1)
    rankboost = job.wrapJobFn(wrap_rankboost, rsem.rv(), merge_mhc.rv(), transgene.rv(),
                              univ_options, tool_options['rankboost'], disk='100M', memory='100M',
                              cores=1)
    # Define the DAG in a static form
    job.addChild(sample_prep)
    # A. The first step is running the alignments and the MHC haplotypers
    sample_prep.addChild(tumor_dna_fqs)
    sample_prep.addChild(normal_dna_fqs)
    sample_prep.addChild(tumor_rna_fqs)

    tumor_rna_fqs.addChild(cutadapt)
    tumor_dna_fqs.addChild(bwa_tumor)
    normal_dna_fqs.addChild(bwa_normal)

    tumor_dna_fqs.addChild(phlat_tumor_dna)
    normal_dna_fqs.addChild(phlat_normal_dna)
    tumor_rna_fqs.addChild(phlat_tumor_rna)
    # B. cutadapt will be followed by star
    cutadapt.addChild(star)
    # Ci.  gene expression and fusion detection follow start alignment
    star.addChild(rsem)
    star.addChild(fusions)
    # Cii.  Radia depends on all 3 alignments
    star.addChild(radia)
    bwa_tumor.addChild(radia)
    bwa_normal.addChild(radia)
    # Ciii. mutect and indel calling depends on dna to have been aligned
    bwa_tumor.addChild(mutect)
    bwa_normal.addChild(mutect)
    bwa_tumor.addChild(muse)
    bwa_normal.addChild(muse)
    bwa_tumor.addChild(somaticsniper)
    bwa_normal.addChild(somaticsniper)
    bwa_tumor.addChild(strelka)
    bwa_normal.addChild(strelka)
    bwa_tumor.addChild(indels)
    bwa_normal.addChild(indels)
    # D. MHC haplotypes will be merged once all 3 samples have been PHLAT-ed
    phlat_tumor_dna.addChild(merge_phlat)
    phlat_normal_dna.addChild(merge_phlat)
    phlat_tumor_rna.addChild(merge_phlat)
    # E. Delete the fastqs from the job store since all alignments are complete
    sample_prep.addChild(fastq_deletion_1)
    cutadapt.addChild(fastq_deletion_1)
    bwa_normal.addChild(fastq_deletion_1)
    bwa_tumor.addChild(fastq_deletion_1)
    phlat_normal_dna.addChild(fastq_deletion_1)
    phlat_tumor_dna.addChild(fastq_deletion_1)
    phlat_tumor_rna.addChild(fastq_deletion_1)
    star.addChild(fastq_deletion_2)
    # F. Mutation calls need to be merged before they can be used
    # G. All mutations get aggregated when they have finished running
    fusions.addChild(merge_mutations)
    radia.addChild(merge_mutations)
    mutect.addChild(merge_mutations)
    muse.addChild(merge_mutations)
    somaticsniper.addChild(merge_mutations)
    strelka.addChild(merge_mutations)
    indels.addChild(merge_mutations)
    # H. Aggregated mutations will be translated to protein space
    merge_mutations.addChild(snpeff)
    # I. snpeffed mutations will be converted into peptides.
    # Transgene also accepts the RNA-seq bam and bai so that it can be rna-aware
    snpeff.addChild(transgene)
    star.addChild(transgene)
    # J. Merged haplotypes and peptides will be converted into jobs and submitted for mhc:peptide
    # binding prediction
    merge_phlat.addChild(spawn_mhc)
    transgene.addChild(spawn_mhc)
    # K. The results from all the predictions will be merged. This is a follow-on job because
    # spawn_mhc will spawn an undetermined number of children.
    spawn_mhc.addFollowOn(merge_mhc)
    # L. Finally, the merged mhc along with the gene expression will be used for rank boosting
    rsem.addChild(rankboost)
    merge_mhc.addChild(rankboost)
    # M. Assess the status of the MHC genes in the patient
    phlat_tumor_rna.addChild(mhc_pathway_assessment)
    rsem.addChild(mhc_pathway_assessment)
    return None


def get_all_tool_inputs(job, tools):
    """
    Iterate through all the tool options and download required files from their remote locations.

    :param dict tools: A dict of dicts of all tools, and their options
    :return: The fully resolved tool dictionary
    :rtype: dict
    """
    for tool in tools:
        for option in tools[tool]:
            if isinstance(tools[tool][option], dict):
                tools[tool][option] = get_all_tool_inputs(job,
                                                          {option: tools[tool][option]})[option]
            else:
                # If a file is of the type file, vcf, tar or fasta, it needs to be downloaded from
                # S3 if reqd, then written to job store.
                if option.split('_')[-1] in ['file', 'vcf', 'index', 'fasta', 'fai', 'idx', 'dict',
                                             'tbi', 'beds']:
                    tools[tool][option] = job.addChildJobFn(get_pipeline_inputs, option,
                                                            tools[tool][option]).rv()
                elif option == 'version':
                    tools[tool][option] = str(tools[tool][option])
    return tools


def get_pipeline_inputs(job, input_flag, input_file, encryption_key=None,
                        per_file_encryption=False):
    """
    Get the input file from s3 or disk and write to file store.

    :param str input_flag: The name of the flag
    :param str input_file: The value passed in the config file
    :param str encryption_key: Path to the encryption key if encrypted with sse-c
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :return: fsID for the file
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    job.fileStore.logToMaster('Obtaining file (%s) to the file job store' %
                              os.path.basename(input_file))
    if input_file.startswith(('http', 'https', 'ftp')):
        input_file = get_file_from_url(job, input_file, encryption_key=encryption_key,
                                       per_file_encryption=per_file_encryption,
                                       write_to_jobstore=False)
    elif input_file.startswith(('S3', 's3')):
        input_file = get_file_from_s3(job, input_file, write_to_jobstore=False,
                                      encryption_key=encryption_key,
                                      per_file_encryption=per_file_encryption)
    else:
        assert os.path.exists(input_file), 'Bogus Input : ' + input_file
    return job.fileStore.writeGlobalFile(input_file)


def prepare_samples(job, fastqs, univ_options):
    """
    Obtain the 6 input fastqs and write them to the file store.  The input files must satisfy the
    following conditions:
    1. Files must be on the form
                    <prefix_for_file>1.<fastq/fq>[.gz]
                    <prefix_for_file>2.<fastq/fq>[.gz]
    2. Forward and reverse reads MUST be in the same folder with the same prefix

    :param dict fastqs: The input fastq dict
           fastqs:
               |- 'tumor_dna': str
               |- 'tumor_rna': str
               |- 'normal_dna': str
               +- 'patient_id': str
    :param dict univ_options: Dict of universal options used by almost all tools
    :return: Updated fastq dict
             sample_fastqs:
                 |-  'tumor_dna': [fsID, fsID]
                 |-  'normal_dna': [fsID, fsID]
                 |-  'tumor_rna': [fsID, fsID]
                 +- 'gzipped': bool
                 +- 'patient_id': str
    :rtype: dict
    """
    job.fileStore.logToMaster('Downloading Inputs for %s' % univ_options['patient'])
    # For each sample type, check if the prefix is an S3 link or a regular file
    # Download S3 files.
    sample_fastqs = {}
    if fastqs['ssec_encrypted']:
        assert univ_options['sse_key'] is not None, 'Cannot read ssec encrypted data without a key.'
    for sample_type in ['tumor_dna', 'tumor_rna', 'normal_dna']:
        if sample_type + '_fastq_2' in fastqs.keys():
            # The user has specified paths to both the _1 and _2 files.
            file_path_1 = fastqs[sample_type + '_fastq_1']
            file_path_2 = fastqs[sample_type + '_fastq_2']
        else:
            prefix, extn = fastqs[''.join([sample_type, '_fastq_1'])], 'temp'
            final_extn = ''
            while extn:
                prefix, extn = os.path.splitext(prefix)
                final_extn = extn + final_extn
                if prefix.endswith('1'):
                    prefix = prefix[:-1]
                    job.fileStore.logToMaster('"%s" prefix for "%s" determined to be %s'
                                              % (sample_type, univ_options['patient'], prefix))
                    break
            else:
                raise ParameterError('Could not determine prefix from provided fastq (%s). Is it '
                                     'of the form <fastq_prefix>1.[fq/fastq][.gz]?'
                                     % fastqs[''.join([sample_type, '_fastq_1'])])
            if final_extn not in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']:
                raise ParameterError('If and _2 fastq path is not specified, only .fastq, .fq or '
                                     'their gzippped extensions are accepted. Could not process '
                                     '%s:%s.' % (univ_options['patient'], sample_type + '_fastq_1'))
            file_path_1 = ''.join([prefix, '1', final_extn])
            file_path_2 = ''.join([prefix, '2', final_extn])
        sample_fastqs[sample_type] = [
            get_pipeline_inputs(job, univ_options['patient'] + ':' + sample_type + '_fastq_1',
                                file_path_1,
                                encryption_key=(univ_options['sse_key']
                                                if fastqs['ssec_encrypted'] else None),
                                per_file_encryption=univ_options['sse_key_is_master']),
            get_pipeline_inputs(job, univ_options['patient'] + ':' + sample_type + '_fastq_2',
                                file_path_2,
                                encryption_key=(univ_options['sse_key']
                                                if fastqs['ssec_encrypted'] else None),
                                per_file_encryption=univ_options['sse_key_is_master'])]
    return sample_fastqs


def get_fqs(job, fastqs, sample_type):
    """
    Convenience function to return only a list of tumor_rna, tumor_dna or normal_dna fqs

    :param dict fastqs: dict of list of fq files
    :param str sample_type: key in sample_type to return
    :return: fastqs[sample_type]
    :rtype: list[toil.fileStore.FileID]
    """
    return fastqs[sample_type]


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
