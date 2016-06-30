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

from urlparse import urlparse

from protect.addons import run_mhc_gene_assessment
from protect.alignment.dna import align_dna
from protect.alignment.rna import align_rna
from protect.binding_prediction.common import spawn_antigen_predictors, merge_mhc_peptide_calls
from protect.common import delete_fastqs, ParameterError
from protect.expression_profiling.rsem import wrap_rsem
from protect.haplotyping.phlat import (merge_phlat_calls,
                                       run_phlat)
from protect.mutation_annotation.snpeff import run_snpeff
from protect.mutation_calling.common import run_mutation_aggregator
from protect.mutation_calling.fusion import run_fusion_caller
from protect.mutation_calling.indel import run_indel_caller
from protect.mutation_calling.muse import run_muse
from protect.mutation_calling.mutect import run_mutect
from protect.mutation_calling.radia import run_radia
from protect.mutation_translation import run_transgene
from protect.qc.rna import run_cutadapt
from protect.rankboost import wrap_rankboost
from toil.job import Job

import argparse
import os
import yaml
import shutil
import subprocess


def parse_config_file(job, config_file, max_cores=None):
    """
    This module will parse the config file withing params and set up the variables that will be
    passed to the various tools in the pipeline.

    ARGUMENTS
    config_file: string containing path to a config file.  An example config
                 file is available at
                        https://s3-us-west-2.amazonaws.com/pimmuno-references
                        /input_parameters.list

    RETURN VALUES
    None
    """
    job.fileStore.logToMaster('Parsing config file')
    config_file = os.path.abspath(config_file)
    if not os.path.exists(config_file):
        raise ParameterError('The config file was not found at specified location. Please verify ' +
                             'and retry.')
    # Initialize variables to hold the sample sets, the universal options, and the per-tool options
    sample_set = defaultdict()
    univ_options = defaultdict()
    tool_options = defaultdict()
    # Read through the notes and the
    with open(config_file, 'r') as conf:
        parsed_config_file = yaml.load(conf.read())

    all_keys = parsed_config_file.keys()
    assert 'patients' in all_keys
    assert 'Universal_Options' in all_keys

    for key in parsed_config_file.keys():
        if key == 'patients':
            sample_set = parsed_config_file['patients']
        elif key == 'Universal_Options':
            univ_options = parsed_config_file['Universal_Options']
            required_options = {'java_Xmx', 'output_folder', 'storage_location'}
            missing_opts = required_options.difference(set(univ_options.keys()))
            if len(missing_opts) > 0:
                raise ParameterError(' The following options have no arguments in the config file:'
                                     '\n' + '\n'.join(missing_opts))
            assert univ_options['storage_location'].startswith(('Local', 'aws'))
            if univ_options['storage_location'].startswith('aws') and 'sse_key' not in univ_options:
                raise ParameterError('Cannot write results to aws without an sse key.')
            if 'sse_key_is_master' in univ_options:
                if univ_options['sse_key_is_master'] not in (True, False):
                    print(univ_options['sse_key_is_master'])
                    print(type(univ_options['sse_key_is_master']))
                    raise ParameterError('sse_key_is_master must be True or False')
            else:
                univ_options['sse_key_is_master'] = False
            if 'sse_key' not in univ_options:
                univ_options['sse_key'] = None
            univ_options['max_cores'] = max_cores
        else:
            tool_options[key] = parsed_config_file[key]
    # Ensure that all tools have been provided options.
    required_tools = {'cutadapt', 'bwa', 'star', 'rsem', 'phlat', 'mut_callers', 'snpeff',
                      'transgene', 'mhci', 'mhcii', 'rank_boost', 'mhc_pathway_assessment'}
    # TODO: Fusions and Indels
    missing_tools = required_tools.difference(set(tool_options.keys()))
    if len(missing_tools) > 0:
        raise ParameterError(' The following tools have no arguments in the config file : \n' +
                             '\n'.join(missing_tools))

    # Get all the tool inputs
    job.fileStore.logToMaster('Obtaining tool inputs')
    process_tool_inputs = job.addChildJobFn(get_all_tool_inputs, tool_options)
    job.fileStore.logToMaster('Obtained tool inputs')
    # Start a job for each sample in the sample set
    for patient_id in sample_set.keys():
        # Add the patient id to the sample set
        sample_set[patient_id]['patient_id'] = patient_id
        job.addFollowOnJobFn(pipeline_launchpad, sample_set[patient_id], univ_options,
                             process_tool_inputs.rv())
    return None


def ascertain_cpu_share(max_cores=None):
    # Ascertain the number of available CPUs. Jobs will be given fractions of this value.
    num_cores = cpu_count()
    # The minimum number of cpus should be at least 6 if possible
    min_cores = min(num_cores, 6)
    cpu_share = max(num_cores / 2, min_cores)
    if max_cores is not None:
        cpu_share = min(cpu_share, max_cores)
    return cpu_share


def pipeline_launchpad(job, fastqs, univ_options, tool_options):
    """
    The precision immuno pipeline begins at this module. The DAG can be viewed in Flowchart.txt

    :param job job: job
    :param dict fastqs: Dict of lists of fastq files
    :param univ_options: Universal Options
    :param tool_options: Options for the various tools
    :return: None
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
    tumor_dna_fqs = job.wrapJobFn(get_fqs, sample_prep.rv(), 'tumor_dna')
    normal_dna_fqs = job.wrapJobFn(get_fqs, sample_prep.rv(), 'normal_dna')
    tumor_rna_fqs = job.wrapJobFn(get_fqs, sample_prep.rv(), 'tumor_rna')
    cutadapt = job.wrapJobFn(run_cutadapt, tumor_rna_fqs.rv(), univ_options,
                             tool_options['cutadapt'], cores=1, disk='80G')
    star = job.wrapJobFn(align_rna, cutadapt.rv(), univ_options, tool_options['star'],
                         cores=tool_options['star']['n'], memory='40G', disk='120G').encapsulate()
    bwa_tumor = job.wrapJobFn(align_dna, tumor_dna_fqs.rv(), 'tumor_dna', univ_options,
                              tool_options['bwa'], cores=tool_options['bwa']['n'],
                              disk='120G').encapsulate()
    bwa_normal = job.wrapJobFn(align_dna, normal_dna_fqs.rv(), 'normal_dna', univ_options,
                               tool_options['bwa'], cores=tool_options['bwa']['n'],
                               disk='120G').encapsulate()
    phlat_tumor_dna = job.wrapJobFn(run_phlat, tumor_dna_fqs.rv(), 'tumor_dna', univ_options,
                                    tool_options['phlat'], cores=tool_options['phlat']['n'],
                                    disk='60G')
    phlat_normal_dna = job.wrapJobFn(run_phlat, normal_dna_fqs.rv(), 'normal_dna', univ_options,
                                     tool_options['phlat'], cores=tool_options['phlat']['n'],
                                     disk='60G')
    phlat_tumor_rna = job.wrapJobFn(run_phlat, tumor_rna_fqs.rv(), 'tumor_rna', univ_options,
                                    tool_options['phlat'], cores=tool_options['phlat']['n'],
                                    disk='60G')
    fastq_deletion_1 = job.wrapJobFn(delete_fastqs, sample_prep.rv())
    fastq_deletion_2 = job.wrapJobFn(delete_fastqs, {'cutadapted_rnas': cutadapt.rv()})
    rsem = job.wrapJobFn(wrap_rsem, star.rv(), univ_options, tool_options['rsem'],
                         cores=tool_options['rsem']['n'], disk='80G').encapsulate()
    mhc_pathway_assessment = job.wrapJobFn(run_mhc_gene_assessment, rsem.rv(), phlat_tumor_rna.rv(),
                                           univ_options, tool_options['mhc_pathway_assessment'])
    fusions = job.wrapJobFn(run_fusion_caller, star.rv(), univ_options, 'fusion_options')
    radia = job.wrapJobFn(run_radia, star.rv(), bwa_tumor.rv(),
                          bwa_normal.rv(), univ_options, tool_options['mut_callers']).encapsulate()
    mutect = job.wrapJobFn(run_mutect, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                           tool_options['mut_callers']).encapsulate()
    muse = job.wrapJobFn(run_muse, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                         tool_options['mut_callers']).encapsulate()
    indels = job.wrapJobFn(run_indel_caller, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                           'indel_options')
    merge_mutations = job.wrapJobFn(run_mutation_aggregator,
                                    {'fusions': fusions.rv(),
                                     'radia': radia.rv(),
                                     'mutect': mutect.rv(),
                                     'indels': indels.rv(),
                                     'muse': muse.rv()}, univ_options).encapsulate()
    snpeff = job.wrapJobFn(run_snpeff, merge_mutations.rv(), univ_options, tool_options['snpeff'],
                           disk='30G')
    transgene = job.wrapJobFn(run_transgene, snpeff.rv(), univ_options, tool_options['transgene'],
                              disk='5G')
    merge_phlat = job.wrapJobFn(merge_phlat_calls, phlat_tumor_dna.rv(), phlat_normal_dna.rv(),
                                phlat_tumor_rna.rv(), univ_options, disk='5G')
    spawn_mhc = job.wrapJobFn(spawn_antigen_predictors, transgene.rv(), merge_phlat.rv(),
                              univ_options, (tool_options['mhci'],
                                             tool_options['mhcii'])).encapsulate()
    merge_mhc = job.wrapJobFn(merge_mhc_peptide_calls, spawn_mhc.rv(), transgene.rv(), univ_options,
                              disk='5G')
    rank_boost = job.wrapJobFn(wrap_rankboost, rsem.rv(), merge_mhc.rv(), transgene.rv(),
                               univ_options, tool_options['rank_boost'], disk='5G')
    # Define the DAG in a static form
    job.addChild(sample_prep)  # Edge  0->1
    # A. The first step is running the alignments and the MHC haplotypers
    sample_prep.addChild(tumor_dna_fqs)  # Edge  1->2
    sample_prep.addChild(normal_dna_fqs)  # Edge  1->2
    sample_prep.addChild(tumor_rna_fqs)  # Edge  1->2

    tumor_rna_fqs.addChild(cutadapt)  # Edge  1->2
    tumor_dna_fqs.addChild(bwa_tumor)  # Edge  1->3
    normal_dna_fqs.addChild(bwa_normal)  # Edge  1->4

    tumor_dna_fqs.addChild(phlat_tumor_dna)  # Edge  1->5
    normal_dna_fqs.addChild(phlat_normal_dna)  # Edge  1->6
    tumor_rna_fqs.addChild(phlat_tumor_rna)  # Edge  1->7
    # B. cutadapt will be followed by star
    cutadapt.addChild(star)  # Edge 2->9
    # Ci.  gene expression and fusion detection follow start alignment
    star.addChild(rsem)  # Edge  9->10
    star.addChild(fusions)  # Edge  9->11
    # Cii.  Radia depends on all 3 alignments
    star.addChild(radia)  # Edge  9->12
    bwa_tumor.addChild(radia)  # Edge  3->12
    bwa_normal.addChild(radia)  # Edge  4->12
    # Ciii. mutect and indel calling depends on dna to have been aligned
    bwa_tumor.addChild(mutect)  # Edge  3->13
    bwa_normal.addChild(mutect)  # Edge  4->13
    bwa_tumor.addChild(muse)  # Edge  3->13
    bwa_normal.addChild(muse)  # Edge  4->13
    bwa_tumor.addChild(indels)  # Edge  3->14
    bwa_normal.addChild(indels)  # Edge  4->14
    # D. MHC haplotypes will be merged once all 3 samples have been PHLAT-ed
    phlat_tumor_dna.addChild(merge_phlat)  # Edge  5->15
    phlat_normal_dna.addChild(merge_phlat)  # Edge  6->15
    phlat_tumor_rna.addChild(merge_phlat)  # Edge  7->15
    # E. Delete the fastqs from the job store since all alignments are complete
    sample_prep.addChild(fastq_deletion_1)  # Edge 1->8
    cutadapt.addChild(fastq_deletion_1)  # Edge 2->8
    bwa_normal.addChild(fastq_deletion_1)  # Edge 3->8
    bwa_tumor.addChild(fastq_deletion_1)  # Edge 4->8
    phlat_normal_dna.addChild(fastq_deletion_1)  # Edge 5->8
    phlat_tumor_dna.addChild(fastq_deletion_1)  # Edge 6>8
    phlat_tumor_rna.addChild(fastq_deletion_1)  # Edge 7->8
    star.addChild(fastq_deletion_2)
    # F. Mutation calls need to be merged before they can be used
    # G. All mutations get aggregated when they have finished running
    fusions.addChild(merge_mutations)  # Edge 11->18
    radia.addChild(merge_mutations)  # Edge 16->18
    mutect.addChild(merge_mutations)  # Edge 17->18
    muse.addChild(merge_mutations)  # Edge 17->18
    indels.addChild(merge_mutations)  # Edge 14->18
    # H. Aggregated mutations will be translated to protein space
    merge_mutations.addChild(snpeff)  # Edge 18->19
    # I. snpeffed mutations will be converted into peptides
    snpeff.addChild(transgene)  # Edge 19->20
    # J. Merged haplotypes and peptides will be converted into jobs and submitted for mhc:peptide
    # binding prediction
    merge_phlat.addChild(spawn_mhc)  # Edge 15->21
    transgene.addChild(spawn_mhc)  # Edge 20->21
    # K. The results from all the predictions will be merged. This is a follow-on job because
    # spawn_mhc will spawn an undetermined number of children.
    spawn_mhc.addFollowOn(merge_mhc)  # Edges 21->XX->22 and 21->YY->22
    # L. Finally, the merged mhc along with the gene expression will be used for rank boosting
    rsem.addChild(rank_boost)  # Edge  10->23
    merge_mhc.addChild(rank_boost)  # Edge 22->23
    # M. Assess the status of the MHC genes in the patient
    phlat_tumor_rna.addChild(mhc_pathway_assessment)  # Edge 7->24
    rsem.addChild(mhc_pathway_assessment)  # Edge 10->24
    return None


def get_all_tool_inputs(job, tools):
    """
    This function will iterate through all the tool options and download the rquired file from their
    remote locations.

    :param dict tools: A dict of dicts of all tools, and their options
    :returns dict: The fully resolved tool dictionary
    """
    for tool in tools:
        for option in tools[tool]:
            # If a file is of the type file, vcf, tar or fasta, it needs to be downloaded from S3 if
            # reqd, then written to job store.
            if option.split('_')[-1] in ['file', 'vcf', 'index', 'fasta', 'fai', 'idx', 'dict',
                                         'tbi']:
                tools[tool][option] = job.addChildJobFn(get_pipeline_inputs, option,
                                                        tools[tool][option]).rv()
    return tools


def get_pipeline_inputs(job, input_flag, input_file):
    """
    Get the input file from s3 or disk, untargz if necessary and then write to file job store.
    :param job: job
    :param str input_flag: The name of the flag
    :param str input_file: The value passed in the config file
    :return: The jobstore ID for the file
    """
    work_dir = os.getcwd()
    job.fileStore.logToMaster('Obtaining file (%s) to the file job store' %
                              os.path.basename(input_file))
    if input_file.startswith('http'):
        assert input_file.startswith('https://s3'), input_file + ' is not an S3 file'
        input_file = get_file_from_s3(job, input_file, write_to_jobstore=False)
    elif input_file.startswith('S3'):
        input_file = get_file_from_s3(job, input_file, write_to_jobstore=False)
    else:
        assert os.path.exists(input_file), 'Bogus Input : ' + input_file
    return job.fileStore.writeGlobalFile(input_file)


def prepare_samples(job, fastqs, univ_options):
    """
    This module will accept a dict object holding the 3 input prefixes and the patient id and will
    attempt to store the fastqs to the jobstore.  The input files must satisfy the following
    1. File extensions can only be fastq or fq (.gz is also allowed)
    2. Forward and reverse reads MUST be in the same folder with the same prefix
    3. The files must be on the form
                    <prefix_for_file>1.<fastq/fq>[.gz]
                    <prefix_for_file>2.<fastq/fq>[.gz]
    The input  dict is:
    tumor_dna: prefix_to_tumor_dna
    tumor_rna: prefix_to_tumor_rna
    normal_dna: prefix_to_normal_dna
    patient_id: patient ID

    The input dict is updated to
    tumor_dna: [jobStoreID_for_fastq1, jobStoreID_for_fastq2]
    tumor_rna: [jobStoreID_for_fastq1, jobStoreID_for_fastq2]
    normal_dna: [jobStoreID_for_fastq1, jobStoreID_for_fastq2]
    patient_id: patient ID
    gzipped: True/False

    This module corresponds to node 1 in the tree
    """
    job.fileStore.logToMaster('Downloading Inputs for %s' % univ_options['patient'])
    allowed_samples = {'tumor_dna_fastq_1', 'tumor_rna_fastq_1', 'normal_dna_fastq_1'}
    if set(fastqs.keys()).difference(allowed_samples) != {'patient_id'}:
        raise ParameterError('Sample with the following parameters has an error:\n' +
                             '\n'.join(fastqs.values()))
    # For each sample type, check if the prefix is an S3 link or a regular file
    # Download S3 files.
    for sample_type in ['tumor_dna', 'tumor_rna', 'normal_dna']:
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
            raise ParameterError('Could not determine prefix from provided fastq (%s). Is it of '
                                 'The form <fastq_prefix>1.[fq/fastq][.gz]?'
                                 % fastqs[''.join([sample_type, '_fastq_1'])])

        # If it is a weblink, it HAS to be from S3
        if prefix.startswith('https://s3') or prefix.lower().startswith('s3://'):
            fastqs[sample_type] = [
                get_file_from_s3(job, ''.join([prefix, '1', final_extn]), univ_options['sse_key'],
                                 per_file_encryption=univ_options['sse_key_is_master']),
                get_file_from_s3(job, ''.join([prefix, '2', final_extn]), univ_options['sse_key'],
                                 per_file_encryption=univ_options['sse_key_is_master'])]
        else:
            # Relies heavily on the assumption that the pair will be in the same
            # folder
            assert os.path.exists(''.join([prefix, '1', final_extn])), \
                'Bogus input: %s' % ''.join([prefix, '1', final_extn])
            # Python lists maintain order hence the values are always guaranteed to be
            # [fastq1, fastq2]
            fastqs[sample_type] = [
                job.fileStore.writeGlobalFile(''.join([prefix, '1', final_extn])),
                job.fileStore.writeGlobalFile(''.join([prefix, '2', final_extn]))]
    [fastqs.pop(x) for x in allowed_samples]
    return fastqs


def get_fqs(job, fastqs, sample_type):
    """
    Convenience function to return only a list of tumor_rna, tumor_dna or normal_dna fqs

    :param job: job
    :param dict fastqs: dict of list of fq files
    :param str sample_type: key in sample_type to return
    :return: fastqs[sample_type]
    """
    return fastqs[sample_type]


def get_file_from_s3(job, s3_url, encryption_key=None, per_file_encryption=True,
                     write_to_jobstore=True):
    """
    Downloads a supplied URL that points to an unencrypted, unprotected file on Amazon S3. The file
    is downloaded and a subsequently written to the jobstore and the return value is a the path to
    the file in the jobstore.

    :param str s3_url: URL for the file (can be s3 or https)
    :param str encryption_key: Path to the master key
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :param bool write_to_jobstore: Should the file be written to the job store?
    """
    work_dir = job.fileStore.getLocalTempDir()

    parsed_url = urlparse(s3_url)
    if parsed_url.scheme == 'https':
        download_url = 'S3:/' + parsed_url.path  # path contains the second /
    elif parsed_url.scheme == 's3':
        download_url = s3_url
    else:
        raise RuntimeError('Unexpected url scheme: %s' % s3_url)

    filename = '/'.join([work_dir, os.path.basename(s3_url)])
    # This is common to encrypted and unencrypted downloads
    download_call = ['s3am', 'download', '--download-exists', 'resume']
    # If an encryption key was provided, use it.
    if encryption_key:
        download_call.extend(['--sse-key-file', encryption_key])
        if per_file_encryption:
            download_call.append('--sse-key-is-master')
    # This is also common to both types of downloads
    download_call.extend([download_url, filename])
    attempt = 0
    exception = ''
    while True:
        try:
            with open(work_dir + '/stderr', 'w') as stderr_file:
                subprocess.check_call(download_call, stderr=stderr_file)
        except subprocess.CalledProcessError:
            # The last line of the stderr will have the error
            with open(stderr_file.name) as stderr_file:
                for line in stderr_file:
                    line = line.strip()
                    if line:
                        exception = line
            if exception.startswith('boto'):
                exception = exception.split(': ')
                if exception[-1].startswith('403'):
                    raise RuntimeError('s3am failed with a "403 Forbidden" error  while obtaining '
                                       '(%s). Did you use the correct credentials?' % s3_url)
                elif exception[-1].startswith('400'):
                    raise RuntimeError('s3am failed with a "400 Bad Request" error while obtaining '
                                       '(%s). Are you trying to download an encrypted file without '
                                       'a key, or an unencrypted file with one?' % s3_url)
                else:
                    raise RuntimeError('s3am failed with (%s) while downloading (%s)' %
                                       (': '.join(exception), s3_url))
            elif exception.startswith('AttributeError'):
                exception = exception.split(': ')
                if exception[-1].startswith("'NoneType'"):
                    raise RuntimeError('Does (%s) exist on s3?' % s3_url)
                else:
                    raise RuntimeError('s3am failed with (%s) while downloading (%s)' %
                                       (': '.join(exception), s3_url))
            else:
                if attempt < 3:
                    attempt += 1
                    continue
                else:
                    raise RuntimeError('Could not diagnose the error while downloading (%s)' %
                                       s3_url)
        except OSError:
            raise RuntimeError('Failed to find "s3am". Install via "apt-get install --pre s3am"')
        else:
            break
        finally:
            os.remove(stderr_file.name)
    assert os.path.exists(filename)
    if write_to_jobstore:
        filename = job.fileStore.writeGlobalFile(filename)
    return filename


def generate_config_file():
    """
    Generate a config file for a ProTECT run on hg19.
    :return: None
    """
    shutil.copy(os.path.join(os.path.dirname(__file__), 'input_parameters.yaml'),
                os.path.join(os.getcwd(), 'ProTECT_config.yaml'))


def main():
    """
    This is the main function for the UCSC Precision Immuno pipeline.
    """
    parser = argparse.ArgumentParser(prog='ProTECT',
                                     description='Prediction of T-Cell Epitopes for Cancer Therapy',
                                     epilog='Contact Arjun Rao (aarao@ucsc.edu) if you encounter '
                                     'any problems while running ProTECT')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--config_file', dest='config_file', help='Config file to be used in the '
                        'run.', type=str, default=None)
    inputs.add_argument('--max-cores-per-job', dest='max_cores', help='Maximum cores to use per '
                        'job. Aligners and Haplotypers ask for cores dependent on the machine that '
                        'the launchpad gets assigned to -- In a heterogeneous cluster, this can '
                        'lead to problems. This value should be set to the number of cpus on the '
                        'smallest node in a cluster.', type=int, required=False, default=None)
    inputs.add_argument('--generate_config', dest='generate_config', help='Generate a config file '
                        'in the current directory that is pre-filled with references and flags for '
                        'an hg19 run.', action='store_true', default=False)
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
        start = Job.wrapJobFn(parse_config_file, params.config_file, params.max_cores)
        Job.Runner.startToil(start, params)
    return None


if __name__ == '__main__':
    main()
