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
File : docker_scripts/precision_immuno.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import print_function
from collections import defaultdict
from encrypt_files_in_dir_to_s3 import write_to_s3
from multiprocessing import cpu_count
from pysam import Samfile
from toil.job import Job

import argparse
import base64
import errno
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time


def parse_config_file(job, config_file):
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
        for line in conf:
            line = line.strip()
            if line.startswith('##') or len(line) == 0:
                continue
            if line.startswith('BEGIN'):
                break
        # The generator function tool_specific_param_generator will yield one group name at a time
        # along with it's parameters.
        for groupname, group_params in tool_specific_param_generator(job, conf):
            if groupname == 'patient':
                if 'patient_id' not in group_params.keys():
                    raise ParameterError('A patient group is missing the patient_id flag.')
                sample_set[group_params['patient_id']] = group_params
            elif groupname == 'Universal_Options':
                univ_options = group_params
                required_options = {'java_Xmx', 'output_folder', 'storage_location'}
                missing_opts = required_options.difference(set(univ_options.keys()))
                if len(missing_opts) > 0:
                    raise ParameterError(' The following options have no arguments in the config '
                                         'file :\n' + '\n'.join(missing_opts))
            # If it isn't any of the above, it's a tool group
            else:
                tool_options[groupname] = group_params
    # Ensure that all tools have been provided options.
    required_tools = {'cutadapt', 'bwa', 'star', 'phlat', 'transgene', 'mut_callers', 'rsem',
                      'mhci', 'mhcii', 'snpeff', 'rank_boost'}
    #                'fusion', 'indels'}
    missing_tools = required_tools.difference(set(tool_options.keys()))
    if len(missing_tools) > 0:
        raise ParameterError(' The following tools have no arguments in the config file : \n' +
                             '\n'.join(missing_tools))
    # Start a job for each sample in the sample set
    for patient_id in sample_set.keys():
        job.addFollowOnJobFn(pipeline_launchpad, sample_set[patient_id], univ_options, tool_options)
    return None


def pipeline_launchpad(job, fastqs, univ_options, tool_options):
    """
    The precision immuno pipeline begins at this module.  The DAG for the pipeline is as follows

                                   0
                                   |
                                   1
                        ___________|______________
                       /        |     |    |   |  \
                      2--+     +3     4+  +5  +6  7+
                      |  |     ||_____||  ||__||__||
                      9  |_ _ _|_ _| _ |_ |_ _|| _ |_ _,8
                  ____|____   _____|__         |
                 /    |    \ /   /    \        |
                |     |     |   |     |        |
                10    11   *12  13   14        15
                |     |     |   |     |        |
                |     |     16  17    |        |
                |     |     |___|     |        |
                |     |       |       |        |
                |     +------18-------+        |
                |             |                |
                |            19                |
                |             |                |
                |            20                |
                |             |________________|
                |                      |
                |                     *21
                |                 _____|_____
                |                /           \
                |               XX           YY
                |               |_____________|
                |                      |
                |                      22
                |                      |
                +----------------------23


     0 = Start Node
     1 = Prepare sample (Download if necessary)
     2 = Process RNA for Adapters (CUTADAPT)
     3 = Align Tumor DNA (BWA)
     4 = Align Normal DNA (BWA)
     5 = Decipher MHC Haplotype Tumor DNA (PHLAT)
     6 = Decipher MHC Haplotype Normal DNA (PHLAT)
     7 = Decipher MHC Haplotype Tumor RNA (PHLAT)
     8 = Delete Fastqs to cleanup space for future jobs
     9 = Align RNA (STAR)
    10 = Calculate Gene Expression (RSEM)
    11 = Identify Fusion Genes (CURRENTLY NOT IMPLEMENTED)
    12 = Mutation Calling 1 (RADIA)
    13 = Mutation Calling 2 (Mutect)
    14 = INDEL Calling (CURRENTLY NOT IMPLEMENTED)
    15 = Merge PHLAT outputs (PYTHON SCRIPT)
    16 = Merge radia mutation calls
    17 = Merge mutect mutation calls
    18 = Merge Mutation calls (PYTHON SCRIPT)
    19 = translate to Protein space (SnpEff, Future Translator)
    20 = Convert AA change to peptides (TRANSGENE)
    21 = Dynamically spawn mhci prediction on n mhci alleles (XX)
         and mhcii prediction on m mhcii alleles (YY)
    XX = Predict MHCI peptides for n predicted alleles
    YY = Predict MHCII peptides for m predicted alleles
    22 = merge MHC:peptide binding predictions
    23 = Rank Boost

     * = Nodes will have dynamically allocated children

    This module corresponds to node 0 on the tree
    """
    # Add Patient id to univ_options as is is passed to every major node in the DAG and can be used
    # as a prefix for the logfile.
    univ_options['patient'] = fastqs['patient_id']
    # Ascertain the number of available CPUs. Jobs will be given fractions of this value.
    ncpu = cpu_count()
    tool_options['star']['n'] = tool_options['bwa']['n'] = tool_options['phlat']['n'] = \
        tool_options['rsem']['n'] = ncpu / 3
    # Define the various nodes in the DAG
    # Need a logfile and a way to send it around
    sample_prep = job.wrapJobFn(prepare_samples, fastqs, univ_options, disk='140G')
    cutadapt = job.wrapJobFn(run_cutadapt, sample_prep.rv(), univ_options, tool_options['cutadapt'],
                             cores=1, disk='80G')
    star = job.wrapJobFn(run_star, cutadapt.rv(), univ_options, tool_options['star'],
                         cores=tool_options['star']['n'], memory='40G', disk='120G').encapsulate()
    bwa_tumor = job.wrapJobFn(run_bwa, sample_prep.rv(), 'tumor_dna', univ_options,
                              tool_options['bwa'], cores=tool_options['bwa']['n'],
                              disk='120G').encapsulate()
    bwa_normal = job.wrapJobFn(run_bwa, sample_prep.rv(), 'normal_dna', univ_options,
                               tool_options['bwa'], cores=tool_options['bwa']['n'],
                               disk='120G').encapsulate()
    phlat_tumor_dna = job.wrapJobFn(run_phlat, sample_prep.rv(), 'tumor_dna', univ_options,
                                    tool_options['phlat'], cores=tool_options['phlat']['n'],
                                    disk='60G')
    phlat_normal_dna = job.wrapJobFn(run_phlat, sample_prep.rv(), 'normal_dna', univ_options,
                                     tool_options['phlat'], cores=tool_options['phlat']['n'],
                                     disk='60G')
    phlat_tumor_rna = job.wrapJobFn(run_phlat, sample_prep.rv(), 'tumor_rna', univ_options,
                                    tool_options['phlat'], cores=tool_options['phlat']['n'],
                                    disk='60G')
    fastq_deletion = job.wrapJobFn(delete_fastqs, sample_prep.rv())
    rsem = job.wrapJobFn(run_rsem, star.rv(), univ_options, tool_options['rsem'],
                         cores=tool_options['rsem']['n'], disk='80G')
    fusions = job.wrapJobFn(run_fusion_caller, star.rv(), univ_options, 'fusion_options')
    Sradia = job.wrapJobFn(spawn_radia, star.rv(), bwa_tumor.rv(),
                           bwa_normal.rv(), univ_options, tool_options['mut_callers']).encapsulate()
    Mradia = job.wrapJobFn(merge_radia, Sradia.rv())
    Smutect = job.wrapJobFn(spawn_mutect, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                            tool_options['mut_callers']).encapsulate()
    Mmutect = job.wrapJobFn(merge_mutect, Smutect.rv())
    indels = job.wrapJobFn(run_indel_caller, bwa_tumor.rv(), bwa_normal.rv(), univ_options,
                           'indel_options')
    merge_mutations = job.wrapJobFn(run_mutation_aggregator, fusions.rv(), Mradia.rv(),
                                    Mmutect.rv(), indels.rv(), univ_options)
    snpeff = job.wrapJobFn(run_snpeff, merge_mutations.rv(), univ_options, tool_options['snpeff'],
                           disk='30G')
    transgene = job.wrapJobFn(run_transgene, snpeff.rv(), univ_options, tool_options['transgene'],
                              disk='5G')
    merge_phlat = job.wrapJobFn(merge_phlat_calls, phlat_tumor_dna.rv(), phlat_normal_dna.rv(),
                                phlat_tumor_rna.rv(), disk='5G')
    spawn_mhc = job.wrapJobFn(spawn_antigen_predictors, transgene.rv(), merge_phlat.rv(),
                              univ_options, (tool_options['mhci'],
                                             tool_options['mhcii'])).encapsulate()
    merge_mhc = job.wrapJobFn(merge_mhc_peptide_calls, spawn_mhc.rv(), transgene.rv(), disk='5G')
    rank_boost = job.wrapJobFn(boost_ranks, rsem.rv(), merge_mhc.rv(), transgene.rv(), univ_options,
                               tool_options['rank_boost'], disk='5G')
    # Define the DAG in a static form
    job.addChild(sample_prep)  # Edge  0->1
    # A. The first step is running the alignments and the MHC haplotypers
    sample_prep.addChild(cutadapt)  # Edge  1->2
    sample_prep.addChild(bwa_tumor)  # Edge  1->3
    sample_prep.addChild(bwa_normal)  # Edge  1->4
    sample_prep.addChild(phlat_tumor_dna)  # Edge  1->5
    sample_prep.addChild(phlat_normal_dna)  # Edge  1->6
    sample_prep.addChild(phlat_tumor_rna)  # Edge  1->7
    # B. cutadapt will be followed by star
    cutadapt.addChild(star)  # Edge 2->8
    # Ci.  gene expression and fusion detection follow start alignment
    star.addChild(rsem)  # Edge  8->9
    star.addChild(fusions)  # Edge  8->10
    # Cii.  Radia depends on all 3 alignments
    star.addChild(Sradia)  # Edge  8->11
    bwa_tumor.addChild(Sradia)  # Edge  3->11
    bwa_normal.addChild(Sradia)  # Edge  4->11
    # Ciii. mutect and indel calling depends on dna to have been aligned
    bwa_tumor.addChild(Smutect)  # Edge  3->12
    bwa_normal.addChild(Smutect)  # Edge  4->12
    bwa_tumor.addChild(indels)  # Edge  3->13
    bwa_normal.addChild(indels)  # Edge  4->13
    # D. MHC haplotypes will be merged once all 3 samples have been PHLAT-ed
    phlat_tumor_dna.addChild(merge_phlat)  # Edge  5->14
    phlat_normal_dna.addChild(merge_phlat)  # Edge  6->14
    phlat_tumor_rna.addChild(merge_phlat)  # Edge  7->14
    # E. Delete the fastqs from the job store since all alignments are complete
    sample_prep.addChild(fastq_deletion)
    cutadapt.addChild(fastq_deletion)
    bwa_normal.addChild(fastq_deletion)
    bwa_tumor.addChild(fastq_deletion)
    phlat_normal_dna.addChild(fastq_deletion)
    phlat_tumor_dna.addChild(fastq_deletion)
    phlat_tumor_rna.addChild(fastq_deletion)
    # F. Mutation calls need to be merged before they can be used
    Sradia.addChild(Mradia)  # Edge 11->15
    Smutect.addChild(Mmutect)  # Edge 12->16
    # G. All mutations get aggregated when they have finished running
    fusions.addChild(merge_mutations)  # Edge 10->17
    Mradia.addChild(merge_mutations)  # Edge 15->17
    Mmutect.addChild(merge_mutations)  # Edge 16->17
    indels.addChild(merge_mutations)  # Edge 13->17
    # H. Aggregated mutations will be translated to protein space
    merge_mutations.addChild(snpeff)  # Edge 17->18
    # I. snpeffed mutations will be converted into peptides
    snpeff.addChild(transgene)  # Edge 18->19
    # J. Merged haplotypes and peptides will be converted into jobs and submitted for mhc:peptide
    # binding prediction
    merge_phlat.addChild(spawn_mhc)  # Edge 14->20
    transgene.addChild(spawn_mhc)  # Edge 19->20
    # K. The results from all the predictions will be merged. This is a follow-on job because
    # spawn_mhc will spawn an undetermined number of children.
    spawn_mhc.addFollowOn(merge_mhc)  # Edges 20->XX->21 and 20->YY->21
    # L. Finally, the merged mhc along with the gene expression will be used for rank boosting
    rsem.addChild(rank_boost)  # Edge  9->22
    merge_mhc.addChild(rank_boost)  # Edge 21->22
    return None


def delete_fastqs(job, fastqs):
    """
    This module will delete the fastqs from teh job Store once their purpose has been achieved (i.e.
    after all mapping steps)

    ARGUMENTS
    1. fastqs: Dict of list of input fastqs
         fastqs
            +- 'tumor_rna': [<JSid for 1.fastq> , <JSid for 2.fastq>]
            +- 'tumor_dna': [<JSid for 1.fastq> , <JSid for 2.fastq>]
            +- 'normal_dna': [<JSid for 1.fastq> , <JSid for 2.fastq>]
    """
    for fq_type in ['tumor_rna', 'tumor_dna', 'normal_dna']:
        for i in xrange(0,2):
            job.fileStore.deleteGlobalFile(fastqs[fq_type][i])
    return None


def run_cutadapt(job, fastqs, univ_options, cutadapt_options):
    """
    This module runs cutadapt on the input RNA fastq files and then calls the RNA aligners.

    ARGUMENTS
    1. fastqs: Dict of list of input RNA-Seq fastqs
         fastqs
            +- 'tumor_rna': [<JSid for 1.fastq> , <JSid for 2.fastq>]
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
              +- 'dockerhub': <dockerhub to use>
    3. cutadapt_options: Dict of parameters specific to cutadapt
         cutadapt_options
              |- 'a': <sequence of 3' adapter to trim from fwd read>
              +- 'A': <sequence of 3' adapter to trim from rev read>
    RETURN VALUES
    1. output_files: Dict of cutadapted fastqs
         output_files
             |- 'rna_cutadapt_1.fastq': <JSid>
             +- 'rna_cutadapt_2.fastq': <JSid>

    This module corresponds to node 2 on the tree
    """
    job.fileStore.logToMaster('Running cutadapt on %s' %univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    fq_extn = '.gz' if fastqs['gzipped'] else ''
    input_files = {
        'rna_1.fastq' + fq_extn: fastqs['tumor_rna'][0],
        'rna_2.fastq' + fq_extn: fastqs['tumor_rna'][1]}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['-a', cutadapt_options['a'],  # Fwd read 3' adapter
                  '-A', cutadapt_options['A'],  # Rev read 3' adapter
                  '-m', '35',  # Minimum size of read
                  '-o', docker_path('rna_cutadapt_1.fastq'),  # Output for R1
                  '-p', docker_path('rna_cutadapt_2.fastq'),  # Output for R2
                  input_files['rna_1.fastq'],
                  input_files['rna_2.fastq']]
    docker_call(tool='cutadapt', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_files = defaultdict()
    for fastq_file in ['rna_cutadapt_1.fastq', 'rna_cutadapt_2.fastq']:
        output_files[fastq_file] = job.fileStore.writeGlobalFile('/'.join([work_dir, fastq_file]))
    return output_files


def run_star(job, fastqs, univ_options, star_options):
    """
    This module uses STAR to align the RNA fastqs to the reference

    ARGUMENTS
    1. fastqs: REFER RETURN VALUE of run_cutadapt()
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
              +- 'dockerhub': <dockerhub to use>
    3. star_options: Dict of parameters specific to STAR
         star_options
             |- 'index_tar': <JSid for the STAR index tarball>
             +- 'n': <number of threads to allocate>
    RETURN VALUES
    1. output_files: Dict of aligned bams
         output_files
             |- 'rnaAligned.toTranscriptome.out.bam': <JSid>
             +- 'rnaAligned.sortedByCoord.out.bam': Dict of genome bam + bai
                                |- 'rna_fix_pg_sorted.bam': <JSid>
                                +- 'rna_fix_pg_sorted.bam.bai': <JSid>

    This module corresponds to node 9 on the tree
    """
    assert star_options['type'] in ('star', 'starlong')
    job.fileStore.logToMaster('Running STAR on %s' %univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'rna_cutadapt_1.fastq': fastqs['rna_cutadapt_1.fastq'],
        'rna_cutadapt_2.fastq': fastqs['rna_cutadapt_2.fastq'],
        'star_index.tar.gz': star_options['index_tar']}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    parameters = ['--runThreadN', str(star_options['n']),
                  '--genomeDir', input_files['star_index'],
                  '--outFileNamePrefix', 'rna',
                  '--readFilesIn',
                  input_files['rna_cutadapt_1.fastq'],
                  input_files['rna_cutadapt_2.fastq'],
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  '--quantMode', 'TranscriptomeSAM',
                  '--outSAMunmapped', 'Within']
    if star_options['type'] == 'star':
        docker_call(tool='star', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'])
    else:
        docker_call(tool='starlong', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'])
    output_files = defaultdict()
    for bam_file in ['rnaAligned.toTranscriptome.out.bam',
                     'rnaAligned.sortedByCoord.out.bam']:
        output_files[bam_file] = job.fileStore.writeGlobalFile('/'.join([
            work_dir, bam_file]))
    job.fileStore.deleteGlobalFile(fastqs['rna_cutadapt_1.fastq'])
    job.fileStore.deleteGlobalFile(fastqs['rna_cutadapt_2.fastq'])
    index_star = job.wrapJobFn(index_bamfile,
                               output_files['rnaAligned.sortedByCoord.out.bam'],
                               'rna', univ_options, disk='120G')
    job.addChild(index_star)
    output_files['rnaAligned.sortedByCoord.out.bam'] = index_star.rv()
    return output_files


def run_bwa(job, fastqs, sample_type, univ_options, bwa_options):
    """
    This module aligns the SAMPLE_TYPE dna fastqs to the reference

    ARGUMENTS -- <ST> depicts the sample type. Substitute with 'tumor'/'normal'
    1. fastqs: Dict of list of input WGS/WXS fastqs
         fastqs
              +- '<ST>_dna': [<JSid for 1.fastq> , <JSid for 2.fastq>]
    2. sample_type: string of 'tumor_dna' or 'normal_dna'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. bwa_options: Dict of parameters specific to bwa
         bwa_options
              |- 'index_tar': <JSid for the bwa index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_files: Dict of aligned bam + reference (nested return)
         output_files
             |- '<ST>_fix_pg_sorted.bam': <JSid>
             +- '<ST>_fix_pg_sorted.bam.bai': <JSid>

    This module corresponds to nodes 3 and 4 on the tree
    """
    job.fileStore.logToMaster('Running bwa on %s:%s' % (univ_options['patient'], sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    fq_extn = '.gz' if fastqs['gzipped'] else ''
    input_files = {
        'dna_1.fastq' + fq_extn: fastqs[sample_type][0],
        'dna_2.fastq' + fq_extn: fastqs[sample_type][1],
        'bwa_index.tar.gz': bwa_options['index_tar']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['mem',
                  '-t', str(bwa_options['n']),
                  '-v', '1',  # Don't print INFO messages to the stderr
                  '/'.join([input_files['bwa_index'], 'hg19.fa']),
                  input_files['dna_1.fastq'],
                  input_files['dna_2.fastq']]
    with open(''.join([work_dir, '/', sample_type, '_aligned.sam']), 'w') as samfile:
        docker_call(tool='bwa', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=samfile)
    # samfile.name retains the path info
    output_file = job.fileStore.writeGlobalFile(samfile.name)
    samfile_processing = job.wrapJobFn(bam_conversion, output_file, sample_type, univ_options,
                                       disk='60G')
    job.addChild(samfile_processing)
    # Return values get passed up the chain to here.  The return value will be a dict with
    # SAMPLE_TYPE_fix_pg_sorted.bam: jobStoreID
    # SAMPLE_TYPE_fix_pg_sorted.bam.bai: jobStoreID
    return samfile_processing.rv()


def bam_conversion(job, samfile, sample_type, univ_options):
    """
    This module converts SAMFILE from sam to bam

    ARGUMENTS
    1. samfile: <JSid for a sam file>
    2. sample_type: string of 'tumor_dna' or 'normal_dna'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    RETURN VALUES
    1. output_files: REFER output_files in run_bwa()
    """
    job.fileStore.logToMaster('Running sam2bam on %s:%s' % (univ_options['patient'], sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'aligned.sam': samfile}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    bamfile = '/'.join([work_dir, 'aligned.bam'])
    parameters = ['view',
                  '-bS',
                  '-o', docker_path(bamfile),
                  input_files['aligned.sam']
                  ]
    docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_file = job.fileStore.writeGlobalFile(bamfile)
    job.fileStore.deleteGlobalFile(samfile)
    reheader_bam = job.wrapJobFn(fix_bam_header, output_file, sample_type, univ_options, disk='60G')
    job.addChild(reheader_bam)
    return reheader_bam.rv()


def fix_bam_header(job, bamfile, sample_type, univ_options):
    """
    This module modified the header in BAMFILE

    ARGUMENTS
    1. bamfile: <JSid for a bam file>
    2. sample_type: string of 'tumor_dna' or 'normal_dna'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    RETURN VALUES
    1. output_files: REFER output_files in run_bwa()
    """
    job.fileStore.logToMaster('Running reheader on %s:%s' % (univ_options['patient'], sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'aligned.bam': bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['view',
                  '-H',
                  input_files['aligned.bam']]
    with open('/'.join([work_dir, 'aligned_bam.header']), 'w') as headerfile:
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=headerfile)
    with open(headerfile.name, 'r') as headerfile, \
            open('/'.join([work_dir, 'output_bam.header']), 'w') as outheaderfile:
        for line in headerfile:
            if line.startswith('@PG'):
                line = '\t'.join([x for x in line.strip().split('\t') if not x.startswith('CL')])
            print(line.strip(), file=outheaderfile)
    parameters = ['reheader',
                  docker_path(outheaderfile.name),
                  input_files['aligned.bam']]
    with open('/'.join([work_dir, 'aligned_fixPG.bam']), 'w') as fixpg_bamfile:
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=fixpg_bamfile)
    output_file = job.fileStore.writeGlobalFile(fixpg_bamfile.name)
    job.fileStore.deleteGlobalFile(bamfile)
    add_rg = job.wrapJobFn(add_readgroups, output_file, sample_type, univ_options, disk='60G')
    job.addChild(add_rg)
    return add_rg.rv()


def add_readgroups(job, bamfile, sample_type, univ_options):
    """
    This module adds the appropriate read groups to the bam file
    ARGUMENTS
    1. bamfile: <JSid for a bam file>
    2. sample_type: string of 'tumor_dna' or 'normal_dna'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                |- 'dockerhub': <dockerhub to use>
                +- 'java_Xmx': value for max heap passed to java
    RETURN VALUES
    1. output_files: REFER output_files in run_bwa()
    """
    job.fileStore.logToMaster('Running add_read_groups on %s:%s' % (univ_options['patient'],
                                                                    sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'aligned_fixpg.bam': bamfile}
    get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['AddOrReplaceReadGroups',
                  'CREATE_INDEX=false',
                  'I=/data/aligned_fixpg.bam',
                  'O=/data/aligned_fixpg_sorted_reheader.bam',
                  'SO=coordinate',
                  'ID=1',
                  ''.join(['LB=', univ_options['patient']]),
                  'PL=ILLUMINA',
                  'PU=12345',
                  ''.join(['SM=', sample_type.rstrip('_dna')])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=univ_options['java_Xmx'])
    output_file = job.fileStore.writeGlobalFile('/'.join([work_dir,
                                                          'aligned_fixpg_sorted_reheader.bam']))
    job.fileStore.deleteGlobalFile(bamfile)
    bam_index = job.wrapJobFn(index_bamfile, output_file, sample_type, univ_options, disk='60G')
    job.addChild(bam_index)
    return bam_index.rv()


def index_bamfile(job, bamfile, sample_type, univ_options):
    """
    This module indexes BAMFILE
    ARGUMENTS
    1. bamfile: <JSid for a bam file>
    2. sample_type: string of 'tumor_dna' or 'normal_dna'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    RETURN VALUES
    1. output_files: REFER output_files in run_bwa(). This module is the one is
                     the one that generates the files.
    """
    job.fileStore.logToMaster('Running samtools-index on %s:%s' % (univ_options['patient'],
                                                                   sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    in_bamfile = '_'.join([sample_type, 'fix_pg_sorted.bam'])
    input_files = {
        in_bamfile: bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    parameters = ['index',
                  input_files[in_bamfile]]
    docker_call(tool='samtools', tool_parameters=parameters,
                work_dir=work_dir, dockerhub=univ_options['dockerhub'])
    output_files = {in_bamfile: bamfile,
                    in_bamfile + '.bai': job.fileStore.writeGlobalFile('/'.join([work_dir,
                                                                                 in_bamfile +
                                                                                 '.bai']))}
    return output_files


def run_rsem(job, star_bams, univ_options, rsem_options):
    """
    This module will run rsem on the RNA Bam file.

    ARGUMENTS
    1. star_bams: Dict of input STAR bams
         star_bams
              +- 'rnaAligned.toTranscriptome.out.bam': <JSid>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. rsem_options: Dict of parameters specific to rsem
         rsem_options
              |- 'index_tar': <JSid for the rsem index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_file: <Jsid of rsem.isoforms.results>

    This module corresponds to node 9 on the tree
    """
    job.fileStore.logToMaster('Running rsem index on %s' % univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'star_transcriptome.bam': star_bams['rnaAligned.toTranscriptome.out.bam'],
        'rsem_index.tar.gz': rsem_options['index_tar']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['--paired-end',
                  '-p', str(rsem_options['n']),
                  '--bam',
                  input_files['star_transcriptome.bam'],
                  '--no-bam-output',
                  '/'.join([input_files['rsem_index'], 'hg19']),
                  'rsem']
    docker_call(tool='rsem', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_file = \
        job.fileStore.writeGlobalFile('/'.join([work_dir, 'rsem.isoforms.results']))
    return output_file


def spawn_radia(job, rna_bam, tumor_bam, normal_bam, univ_options, radia_options):
    """
    This module will spawn a radia job for each chromosome, on the RNA and DNA.

    ARGUMENTS
    1. rna_bam: Dict of input STAR bams
         rna_bam
              |- 'rnaAligned.sortedByCoord.out.bam': REFER run_star()
                                |- 'rna_fix_pg_sorted.bam': <JSid>
                                +- 'rna_fix_pg_sorted.bam.bai': <JSid>
    2. tumor_bam: Dict of input tumor WGS/WSQ bam + bai
         tumor_bam
              |- 'tumor_fix_pg_sorted.bam': <JSid>
              +- 'tumor_fix_pg_sorted.bam.bai': <JSid>
    3. normal_bam: Dict of input normal WGS/WSQ bam + bai
         normal_bam
              |- 'normal_fix_pg_sorted.bam': <JSid>
              +- 'normal_fix_pg_sorted.bam.bai': <JSid>
    4. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    5. radia_options: Dict of parameters specific to radia
         radia_options
              |- 'genome_fasta': <JSid for genome fasta file>
              +- 'genome_fai': <JSid for genome fai file>

    RETURN VALUES
    1. perchrom_radia: Dict of results of radia per chromosome
         perchrom_radia
              |- 'chr1'
              |   |- 'radia_filtered_chr1.vcf': <JSid>
              |   +- 'radia_filtered_chr1_radia.log': <JSid>
              |- 'chr2'
              |   |- 'radia_filtered_chr2.vcf': <JSid>
              |   +- 'radia_filtered_chr2_radia.log': <JSid>
             etc...

    This module corresponds to node 11 on the tree
    """
    job.fileStore.logToMaster('Running spawn_radia on %s' % univ_options['patient'])
    rna_bam_key = 'rnaAligned.sortedByCoord.out.bam' # to reduce next line size
    bams = {'tumor_rna': rna_bam[rna_bam_key]['rna_fix_pg_sorted.bam'],
            'tumor_rnai': rna_bam[rna_bam_key]['rna_fix_pg_sorted.bam.bai'],
            'tumor_dna': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
            'tumor_dnai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
            'normal_dna': normal_bam['normal_dna_fix_pg_sorted.bam'],
            'normal_dnai': normal_bam['normal_dna_fix_pg_sorted.bam.bai']}
    # Make a dict object to hold the return values for each of the chromosome jobs.  Then run radia
    # on each chromosome.
    chromosomes = [''.join(['chr', str(x)]) for x in range(1, 23) + ['X', 'Y']]
    perchrom_radia = defaultdict()
    for chrom in chromosomes:
        perchrom_radia[chrom] = job.addChildJobFn(run_radia, bams, univ_options, radia_options,
                                                  chrom, disk='60G').rv()
    return perchrom_radia


def merge_radia(job, perchrom_rvs):
    """
    This module will merge the per-chromosome radia files created by spawn_radia into a genome vcf.
    It will make 2 vcfs, one for PASSing non-germline calls, and one for all calls.

    ARGUMENTS
    1. perchrom_rvs: REFER RETURN VALUE of spawn_radia()

    RETURN VALUES
    1. output_files: Dict of outputs
            output_files
                |- radia_calls.vcf: <JSid>
                +- radia_parsed_filter_passing_calls.vcf: <JSid>

    This module corresponds to node 11 on the tree
    """
    job.fileStore.logToMaster('Running merge_radia')
    work_dir = job.fileStore.getLocalTempDir()
    # We need to squash the input dict of dicts to a single dict such that it can be passed to
    # get_files_from_filestore
    input_files = {filename: jsid for perchrom_files in perchrom_rvs.values()
                   for filename, jsid in perchrom_files.items()}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=False)
    chromosomes = [''.join(['chr', str(x)]) for x in range(1, 23) + ['X', 'Y']]
    with open('/'.join([work_dir, 'radia_calls.vcf']), 'w') as radfile, \
            open('/'.join([work_dir, 'radia_filter_passing_calls.vcf']), 'w') as radpassfile:
        for chrom in chromosomes:
            with open(input_files[''.join(['radia_filtered_', chrom, '.vcf'])], 'r') as filtradfile:
                for line in filtradfile:
                    line = line.strip()
                    if line.startswith('#'):
                        if chrom == 'chr1':
                            print(line, file=radfile)
                            print(line, file=radpassfile)
                        continue
                    else:
                        print(line, file=radfile)
                        line = line.split('\t')
                        if line[6] == 'PASS' and 'MT=GERM' not in line[7]:
                            print('\t'.join(line), file=radpassfile)
    # parse the PASS radia vcf for multiple alt alleles
    with open(radpassfile.name, 'r') as radpassfile, \
            open('/'.join([work_dir, 'radia_parsed_filter_passing_calls.vcf']),
                 'w') as parsedradfile:
        parse_radia_multi_alt(radpassfile, parsedradfile)
    output_files = defaultdict()
    for radia_file in [radfile.name, parsedradfile.name]:
        output_files[os.path.basename(radia_file)] = job.fileStore.writeGlobalFile(radia_file)
    return output_files


def run_radia(job, bams, univ_options, radia_options, chrom):
    """
    This module will run radia on the RNA and DNA bams

    ARGUMENTS
    1. bams: Dict of bams and their indexes
        bams
         |- 'tumor_rna': <JSid>
         |- 'tumor_rnai': <JSid>
         |- 'tumor_dna': <JSid>
         |- 'tumor_dnai': <JSid>
         |- 'normal_dna': <JSid>
         +- 'normal_dnai': <JSid>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. radia_options: Dict of parameters specific to radia
         radia_options
              |- 'dbsnp_vcf': <JSid for dnsnp vcf file>
              +- 'genome': <JSid for genome fasta file>
    4. chrom: String containing chromosome name with chr appended

    RETURN VALUES
    1. Dict of filtered radia output vcf and logfile (Nested return)
        |- 'radia_filtered_CHROM.vcf': <JSid>
        +- 'radia_filtered_CHROM_radia.log': <JSid>
    """
    job.fileStore.logToMaster('Running radia on %s:%s' %(univ_options['patient'], chrom))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'rna.bam': bams['tumor_rna'],
        'rna.bam.bai': bams['tumor_rnai'],
        'tumor.bam': bams['tumor_dna'],
        'tumor.bam.bai': bams['tumor_dnai'],
        'normal.bam': bams['normal_dna'],
        'normal.bam.bai': bams['normal_dnai'],
        'genome.fasta': radia_options['genome_fasta'],
        'genome.fasta.fai': radia_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    radia_output = ''.join([work_dir, '/radia_', chrom, '.vcf'])
    radia_log = ''.join([work_dir, '/radia_', chrom, '_radia.log'])
    parameters = [univ_options['patient'],  # shortID
                  chrom,
                  '-n', input_files['normal.bam'],
                  '-t', input_files['tumor.bam'],
                  '-r', input_files['rna.bam'],
                  ''.join(['--rnaTumorFasta=', input_files['genome.fasta']]),
                  '-f', input_files['genome.fasta'],
                  '-o', docker_path(radia_output),
                  '-i', 'hg19_M_rCRS',
                  '-m', input_files['genome.fasta'],
                  '-d', 'aarjunrao@soe.ucsc.edu',
                  '-q', 'Illumina',
                  '--disease', 'CANCER',
                  '-l', 'INFO',
                  '-g', docker_path(radia_log)]
    docker_call(tool='radia', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_files = defaultdict()
    for radia_file in [radia_output, radia_log]:
        output_files[os.path.basename(radia_file)] = \
            job.fileStore.writeGlobalFile(radia_file)
    filterradia = job.wrapJobFn(run_filter_radia, bams,
                                output_files[os.path.basename(radia_output)],
                                univ_options, radia_options, chrom, disk='60G')
    job.addChild(filterradia)
    return filterradia.rv()


def run_filter_radia(job, bams, radia_file, univ_options, radia_options, chrom):
    """
    This module will run filterradia on the RNA and DNA bams.

    ARGUMENTS
    1. bams: REFER ARGUMENTS of run_radia()
    2. univ_options: REFER ARGUMENTS of run_radia()
    3. radia_file: <JSid of vcf generated by run_radia()>
    3. radia_options: REFER ARGUMENTS of run_radia()
    4. chrom: REFER ARGUMENTS of run_radia()

    RETURN VALUES
    1. Dict of filtered radia output vcf and logfile
        |- 'radia_filtered_CHROM.vcf': <JSid>
        +- 'radia_filtered_CHROM_radia.log': <JSid>
    """
    job.fileStore.logToMaster('Running filter-radia on %s:%s' % (univ_options['patient'], chrom))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'rna.bam': bams['tumor_rna'],
        'rna.bam.bai': bams['tumor_rnai'],
        'tumor.bam': bams['tumor_dna'],
        'tumor.bam.bai': bams['tumor_dnai'],
        'normal.bam': bams['normal_dna'],
        'normal.bam.bai': bams['normal_dnai'],
        'radia.vcf': radia_file,
        'genome.fasta': radia_options['genome_fasta'],
        'genome.fasta.fai': radia_options['genome_fai']
        }
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    filterradia_output = ''.join(['radia_filtered_', chrom, '.vcf'])
    filterradia_log = ''.join([work_dir, '/radia_filtered_', chrom, '_radia.log'
                              ])
    parameters = [univ_options['patient'],  # shortID
                  chrom.lstrip('chr'),
                  input_files['radia.vcf'],
                  '/data',
                  '/home/radia/scripts',
                  '-b', '/home/radia/data/hg19/blacklists/1000Genomes/phase1/',
                  '-d', '/home/radia/data/hg19/snp135',
                  '-r', '/home/radia/data/hg19/retroGenes/',
                  '-p', '/home/radia/data/hg19/pseudoGenes/',
                  '-c', '/home/radia/data/hg19/cosmic/',
                  '-t', '/home/radia/data/hg19/gaf/2_1',
                  '--noSnpEff',
                  '--rnaGeneBlckFile', '/home/radia/data/rnaGeneBlacklist.tab',
                  '--rnaGeneFamilyBlckFile',
                  '/home/radia/data/rnaGeneFamilyBlacklist.tab',
                  '-f', input_files['genome.fasta'],
                  '--log=INFO',
                  '-g', docker_path(filterradia_log)]
    docker_call(tool='filterradia', tool_parameters=parameters,
                work_dir=work_dir, dockerhub=univ_options['dockerhub'])
    output_files = defaultdict()
    output_files[filterradia_output] = \
        job.fileStore.writeGlobalFile(''.join([work_dir, '/',
                                               univ_options['patient'], '_',
                                               chrom, '.vcf']))
    output_files[os.path.basename(filterradia_log)] = \
        job.fileStore.writeGlobalFile(filterradia_log)
    return output_files


def spawn_mutect(job, tumor_bam, normal_bam, univ_options, mutect_options):
    """
    This module will spawn a mutect job for each chromosome on the DNA bams.

    ARGUMENTS
    1. tumor_bam: Dict of input tumor WGS/WSQ bam + bai
         tumor_bam
              |- 'tumor_fix_pg_sorted.bam': <JSid>
              +- 'tumor_fix_pg_sorted.bam.bai': <JSid>
    2. normal_bam: Dict of input normal WGS/WSQ bam + bai
         normal_bam
              |- 'normal_fix_pg_sorted.bam': <JSid>
              +- 'normal_fix_pg_sorted.bam.bai': <JSid>
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. mutect_options: Dict of parameters specific to mutect
         mutect_options
              |- 'dbsnp_vcf': <JSid for dnsnp vcf file>
              |- 'dbsnp_idx': <JSid for dnsnp vcf index file>
              |- 'cosmic_vcf': <JSid for cosmic vcf file>
              |- 'cosmic_idx': <JSid for cosmic vcf index file>
              +- 'genome_fasta': <JSid for genome fasta file>

    RETURN VALUES
    1. perchrom_mutect: Dict of results of mutect per chromosome
         perchrom_mutect
              |- 'chr1'
              |   +- 'mutect_chr1.vcf': <JSid>
              |   +- 'mutect_chr1.out': <JSid>
              |- 'chr2'
              |   |- 'mutect_chr2.vcf': <JSid>
              |   +- 'mutect_chr2.out': <JSid>
             etc...

    This module corresponds to node 11 on the tree
    """
    job.fileStore.logToMaster('Running spawn_mutect on %s' % univ_options['patient'])
    # Make a dict object to hold the return values for each of the chromosome
    # jobs.  Then run mutect on each chromosome.
    chromosomes = [''.join(['chr', str(x)]) for x in range(1, 23) + ['X', 'Y']]
    perchrom_mutect = defaultdict()
    for chrom in chromosomes:
        perchrom_mutect[chrom] = job.addChildJobFn(run_mutect, tumor_bam, normal_bam, univ_options,
                                                   mutect_options, chrom, disk='60G',
                                                   memory='3.5G').rv()
    return perchrom_mutect


def merge_mutect(job, perchrom_rvs):
    """
    This module will merge the per-chromosome mutect files created by spawn_mutect into a genome
    vcf.  It will make 2 vcfs, one for PASSing non-germline calls, and one for all calls.

    ARGUMENTS
    1. perchrom_rvs: REFER RETURN VALUE of spawn_mutect()

    RETURN VALUES
    1. output_files: <JSid for mutect_passing_calls.vcf>

    This module corresponds to node 11 on the tree
    """
    job.fileStore.logToMaster('Running merge_mutect')
    work_dir = job.fileStore.getLocalTempDir()
    # We need to squash the input dict of dicts to a single dict such that it can be passed to
    # get_files_from_filestore
    input_files = {filename: jsid for perchrom_files in perchrom_rvs.values()
                   for filename, jsid in perchrom_files.items()}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    chromosomes = [''.join(['chr', str(x)]) for x in range(1, 23) + ['X', 'Y']]
    with open('/'.join([work_dir, 'mutect_calls.vcf']), 'w') as mutvcf, \
            open('/'.join([work_dir, 'mutect_calls.out']), 'w') as mutout, \
            open('/'.join([work_dir, 'mutect_passing_calls.vcf']), 'w') as mutpassvcf:
        out_header_not_printed = True
        for chrom in chromosomes:
            with open(input_files[''.join(['mutect_', chrom, '.vcf'])], 'r') as mutfile:
                for line in mutfile:
                    line = line.strip()
                    if line.startswith('#'):
                        if chrom == 'chr1':
                            print(line, file=mutvcf)
                            print(line, file=mutpassvcf)
                        continue
                    else:
                        print(line, file=mutvcf)
                        line = line.split('\t')
                        if line[6] != 'REJECT':
                            print('\t'.join(line), file=mutpassvcf)
            with open(input_files[''.join(['mutect_', chrom,
                                           '.out'])], 'r') as mutfile:
                for line in mutfile:
                    line = line.strip()
                    if line.startswith('#'):
                        if chrom == 'chr1':
                            print(line, file=mutout)
                        continue
                    elif out_header_not_printed:
                        print(line, file=mutout)
                        out_header_not_printed = False
                    else:
                        print(line, file=mutout)
    output_file = job.fileStore.writeGlobalFile(mutpassvcf.name)
    return output_file


def run_mutect(job, tumor_bam, normal_bam, univ_options, mutect_options, chrom):
    """
    This module will run mutect on the DNA bams

    ARGUMENTS
    1. tumor_bam: REFER ARGUMENTS of spawn_mutect()
    2. normal_bam: REFER ARGUMENTS of spawn_mutect()
    3. univ_options: REFER ARGUMENTS of spawn_mutect()
    4. mutect_options: REFER ARGUMENTS of spawn_mutect()
    5. chrom: String containing chromosome name with chr appended

    RETURN VALUES
    1. output_files: Dict of results of mutect for chromosome
            output_files
              |- 'mutect_CHROM.vcf': <JSid>
              +- 'mutect_CHROM.out': <JSid>

    This module corresponds to node 12 on the tree
    """
    job.fileStore.logToMaster('Running mutect on %s:%s' % (univ_options['patient'], chrom))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'normal.bam': normal_bam['normal_dna_fix_pg_sorted.bam'],
        'normal.bam.bai': normal_bam['normal_dna_fix_pg_sorted.bam.bai'],
        'genome.fa': mutect_options['genome_fasta'],
        'genome.fa.fai': mutect_options['genome_fai'],
        'genome.dict': mutect_options['genome_dict'],
        'cosmic.vcf': mutect_options['cosmic_vcf'],
        'cosmic.vcf.idx': mutect_options['cosmic_idx'],
        'dbsnp.vcf': mutect_options['dbsnp_vcf'],
        'dbsnp.vcf.idx': mutect_options['dbsnp_idx']}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    mutout = ''.join([work_dir, '/mutect_', chrom, '.out'])
    mutvcf = ''.join([work_dir, '/mutect_', chrom, '.vcf'])
    parameters = ['-R', input_files['genome.fa'],
                  '--cosmic', input_files['cosmic.vcf'],
                  '--dbsnp', input_files['dbsnp.vcf'],
                  '--input_file:normal', input_files['normal.bam'],
                  '--input_file:tumor', input_files['tumor.bam'],
                  #'--tumor_lod', str(10),
                  #'--initial_tumor_lod', str(4.0),
                  '-L', chrom,
                  '--out', docker_path(mutout),
                  '--vcf', docker_path(mutvcf)
                 ]
    Xmx = mutect_options['java_Xmx'] if mutect_options['java_Xmx'] else univ_options['java_Xmx']
    docker_call(tool='mutect:1.1.7', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=Xmx)
    output_files = defaultdict()
    for mutect_file in [mutout, mutvcf]:
        output_files[os.path.basename(mutect_file)] = job.fileStore.writeGlobalFile(mutect_file)
    return output_files


def run_indel_caller(job, tumor_bam, normal_bam, univ_options, indel_options):
    """
    This module will run an indel caller on the DNA bams.  This module will be
    implemented in the future.

    This module corresponds to node 13 on the tree
    """
    job.fileStore.logToMaster('Running INDEL on %s' % univ_options['patient'])
    indel_file = job.fileStore.getLocalTempFile()
    output_file = job.fileStore.writeGlobalFile(indel_file)
    return output_file


def run_fusion_caller(job, star_bam, univ_options, fusion_options):
    """
    This module will run a fusion caller on DNA bams.  This module will be
    implemented in the future.

    This module corresponds to node 10 on the tree
    """
    job.fileStore.logToMaster('Running FUSION on %s' % univ_options['patient'])
    fusion_file = job.fileStore.getLocalTempFile()
    output_file = job.fileStore.writeGlobalFile(fusion_file)
    return output_file


def run_mutation_aggregator(job, fusion_output, radia_output, mutect_output, indel_output,
                            univ_options):
    """
    This module will aggregate all the mutations called in the previous steps and will then call
    snpeff on the results.

    ARGUMENTS
    1. fusion_output: <JSid for vcf generated by the fusion caller>
    2. radia_output: <JSid for vcf generated by radia>
    3. mutect_output: <JSid for vcf generated by mutect>
    4. indel_output: <JSid for vcf generated by the indel caller>

    RETURN VALUES
    1. output_file: <JSid for merged vcf>

    This module corresponds to node 15 on the tree
    """
    job.fileStore.logToMaster('Aggregating mutations for %s' % univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'mutect.vcf': mutect_output,
        'radia.vcf': radia_output['radia_parsed_filter_passing_calls.vcf'],
        'indel.vcf': indel_output,
        'fusion.vcf': fusion_output}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    # Modify these once INDELs and Fusions are implemented
    input_files.pop('indel.vcf')
    input_files.pop('fusion.vcf')
    # read files into memory
    vcf_file = defaultdict()
    mutcallers = input_files.keys()
    with open(''.join([work_dir, '/', univ_options['patient'], '_merged_mutations.vcf']),
              'w') as merged_mut_file:
        for mut_caller in mutcallers:
            caller = mut_caller.rstrip('.vcf')
            vcf_file[caller] = defaultdict()
            with open(input_files[mut_caller], 'r') as mutfile:
                for line in mutfile:
                    if line.startswith('#'):
                        if caller == 'radia':
                            print(line.strip(), file=merged_mut_file)
                        continue
                    line = line.strip().split()
                    vcf_file[caller][(line[0], line[1], line[3], line[4])] = line
    # This method can be changed in the future to incorporate more callers and
    # fancier integration methods
    merge_vcfs(vcf_file, merged_mut_file.name)
    export_results(merged_mut_file.name, univ_options)
    output_file = job.fileStore.writeGlobalFile(merged_mut_file.name)
    return output_file


def run_snpeff(job, merged_mutation_file, univ_options, snpeff_options):
    """
    This module will run snpeff on the aggregated mutation calls.  Currently the only mutations
    called are SNPs hence SnpEff suffices. This node will be replaced in the future with another
    translator.

    ARGUMENTS
    1. merged_mutation_file: <JSid for merged vcf>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. snpeff_options: Dict of parameters specific to snpeff
         snpeff_options
                +- 'index_tar': <JSid for the snpEff index tarball>

    RETURN VALUES
    1. output_file: <JSid for the snpeffed vcf>

    This node corresponds to node 16 on the tree
    """
    job.fileStore.logToMaster('Running snpeff on %s' % univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'merged_mutations.vcf': merged_mutation_file,
        'snpeff_index.tar.gz': snpeff_options['index_tar']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['eff',
                  '-dataDir', input_files['snpeff_index'],
                  '-c', '/'.join([input_files['snpeff_index'], 'snpEff_hg19_gencode.config']),
                  '-no-intergenic',
                  '-no-downstream',
                  '-no-upstream',
                  #'-canon',
                  '-noStats',
                  'hg19_gencode',
                  input_files['merged_mutations.vcf']]
    Xmx = snpeff_options['java_Xmx'] if snpeff_options['java_Xmx'] else univ_options['java_Xmx']
    with open('/'.join([work_dir, 'snpeffed_mutations.vcf']), 'w') as snpeff_file:
        docker_call(tool='snpeff', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], java_opts=Xmx, outfile=snpeff_file)
    output_file = job.fileStore.writeGlobalFile(snpeff_file.name)
    return output_file


def run_transgene(job, snpeffed_file, univ_options, transgene_options):
    """
    This module will run transgene on the input vcf file from the aggregator and produce the
    peptides for MHC prediction

    ARGUMENTS
    1. snpeffed_file: <JSid for snpeffed vcf>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. transgene_options: Dict of parameters specific to transgene
         transgene_options
                +- 'gencode_peptide_fasta': <JSid for the gencode protein fasta>

    RETURN VALUES
    1. output_files: Dict of transgened n-mer peptide fastas
         output_files
                |- 'transgened_tumor_9_mer_snpeffed.faa': <JSid>
                |- 'transgened_tumor_10_mer_snpeffed.faa': <JSid>
                +- 'transgened_tumor_15_mer_snpeffed.faa': <JSid>

    This module corresponds to node 17 on the tree
    """
    job.fileStore.logToMaster('Running transgene on %s' % univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'snpeffed_muts.vcf': snpeffed_file,
        'pepts.fa': transgene_options['gencode_peptide_fasta']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['--peptides', input_files['pepts.fa'],
                  '--snpeff', input_files['snpeffed_muts.vcf'],
                  '--prefix', 'transgened',
                  '--pep_lens', '9,10,15']
    docker_call(tool='transgene', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_files = defaultdict()
    for peplen in ['9', '10', '15']:
        peptfile = '_'.join(['transgened_tumor', peplen, 'mer_snpeffed.faa'])
        mapfile = '_'.join(['transgened_tumor', peplen, 'mer_snpeffed.faa.map'])
        output_files[peptfile] = job.fileStore.writeGlobalFile(os.path.join(work_dir, peptfile))
        output_files[mapfile] = job.fileStore.writeGlobalFile(os.path.join(work_dir, mapfile))
    return output_files


def run_phlat(job, fastqs, sample_type, univ_options, phlat_options):
    """
    This module will run PHLAT on SAMPLE_TYPE fastqs.

    ARGUMENTS -- <ST> depicts the sample type. Substitute with 'tumor_dna',
                 'normal_dna', or 'tumor_rna'
    1. fastqs: Dict of list of input WGS/WXS fastqs
         fastqs
              +- '<ST>': [<JSid for 1.fastq> , <JSid for 2.fastq>]
    2. sample_type: string of 'tumor' or 'normal'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. phlat_options: Dict of parameters specific to phlat
         phlat_options
              |- 'index_tar': <JSid for the PHLAT index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_file: <JSid for the allele predictions for ST>

    This module corresponds to nodes 5, 6 and 7 on the tree
    """
    job.fileStore.logToMaster('Running phlat on %s:%s' % (univ_options['patient'], sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    fq_extn = '.gz' if fastqs['gzipped'] else ''
    input_files = {
        'input_1.fastq' + fq_extn: fastqs[sample_type][0],
        'input_2.fastq' + fq_extn: fastqs[sample_type][1],
        'phlat_index.tar.gz': phlat_options['index_tar']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['-1', input_files['input_1.fastq'],
                  '-2', input_files['input_2.fastq'],
                  '-index', input_files['phlat_index'],
                  '-b2url', '/usr/local/bin/bowtie2',
                  '-tag', sample_type,
                  '-e', '/home/phlat-1.0',  # Phlat directory home
                  '-o', '/data',  # Output directory
                  '-p', str(phlat_options['n'])]  # Number of threads
    docker_call(tool='phlat', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_file = job.fileStore.writeGlobalFile(''.join([work_dir, '/', sample_type, '_HLA.sum']))
    return output_file


def merge_phlat_calls(job, tumor_phlat, normal_phlat, rna_phlat):
    """
    This module will merge the results form running PHLAT on the 3 input fastq
    pairs.

    ARGUMENTS
    1. tumor_phlat: <JSid for tumor DNA called alleles>
    2. normal_phlat: <JSid for normal DNA called alleles>
    3. rna_phlat: <JSid for tumor RNA called alleles>

    RETURN VALUES
    1. output_files: Dict of JSids for consensus MHCI and MHCII alleles
             output_files
                    |- 'mhci_alleles.list': <JSid>
                    +- 'mhcii_alleles.list': <JSid>

    This module corresponds to node 14 on the tree
    """
    job.fileStore.logToMaster('Merging Phlat calls')
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'tumor_dna': tumor_phlat,
        'normal_dna': normal_phlat,
        'tumor_rna': rna_phlat}
    input_files = get_files_from_filestore(job, input_files, work_dir)
    with open(input_files['tumor_dna'], 'r') as td_file, \
            open(input_files['normal_dna'], 'r') as nd_file, \
            open(input_files['tumor_rna'], 'r') as tr_file:
        # TODO: Could this be a defautdict?
        mhc_alleles = {'HLA_A': [], 'HLA_B': [], 'HLA_C': [], 'HLA_DPA': [], 'HLA_DQA': [],
                       'HLA_DPB': [], 'HLA_DQB': [], 'HLA_DRB': []}
        for phlatfile in td_file, nd_file, tr_file:
            for line in phlatfile:
                if line.startswith('Locus'):
                    continue
                line = line.strip().split()
                if line[0].startswith('HLA_D'):
                    line[0] = line[0][:-1]  # strip the last character
                # Sometimes we get an error saying there was insufficient read
                # converage. We need to drop that line.
                # E.g. HLA_DQB1 no call due to insufficient reads at this locus
                if line[1] == 'no':
                    continue
                if line[4] != 'NA' and not line[1].endswith('xx'):
                    mhc_alleles[line[0]].append((line[1], line[4]))
                if line[5] != 'NA' and not line[2].endswith('xx'):
                    mhc_alleles[line[0]].append((line[2], line[5]))
    # Get most probable alleles for each allele group and print to output
    with open(os.path.join(work_dir, 'mhci_alleles.list'), 'w') as mhci_file, \
            open(os.path.join(work_dir, 'mhcii_alleles.list'), 'w') as mhcii_file:
        for mhci_group in ['HLA_A', 'HLA_B', 'HLA_C']:
            mpa = most_probable_alleles(mhc_alleles[mhci_group])
            print('\n'.join([''.join(['HLA-', x]) for x in mpa]), file=mhci_file)
        drb_mpa = most_probable_alleles(mhc_alleles['HLA_DRB'])
        print('\n'.join([''.join(['HLA-', x]) for x in drb_mpa]), file=mhcii_file)
        dqa_mpa = most_probable_alleles(mhc_alleles['HLA_DQA'])
        dqb_mpa = most_probable_alleles(mhc_alleles['HLA_DQB'])
        for dqa_allele in dqa_mpa:
            for dqb_allele in dqb_mpa:
                print(''.join(['HLA-', dqa_allele, '/', dqb_allele]), file=mhcii_file)
    output_files = defaultdict()
    for allele_file in ['mhci_alleles.list', 'mhcii_alleles.list']:
        output_files[allele_file] = job.fileStore.writeGlobalFile(os.path.join(work_dir,
                                                                               allele_file))
    return output_files


def spawn_antigen_predictors(job, transgened_files, phlat_files, univ_options, mhc_options):
    """
    Based on the number of alleles obtained from node 14, this module will spawn callers to predict
    MHCI:peptide and MHCII:peptide binding on the peptide files from node 17.  Once all MHC:peptide
    predictions are made, merge them via a follow-on job.

    ARGUMENTS
    1. transgened_files: REFER RETURN VALUE of run_transgene()
    2. phlat_files: REFER RETURN VALUE of merge_phlat_calls()
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. mhc_options: Dict of dicts of parameters specific to mhci and mhcii
                    respectively
         mhc_options
              |- 'mhci'
              |     |- 'method_file': <JSid for json file containing data
              |     |                  linking alleles, peptide lengths, and
              |     |                  prediction methods>
              |     +- 'pred': String describing prediction method to use
              +- 'mhcii'
                    |- 'method_file': <JSid for json file containing data
                    |                  linking alleles and prediction methods>
                    +- 'pred': String describing prediction method to use

    RETURN VALUES
    1. tuple of (mhci_preds, mhcii_preds)
         mhci_preds: Dict of return value from running predictions on a given
                     mhc for all peptides of length 9 and 10.
             mhci_preds
                |- <MHC molecule 1>_9_mer.faa: <PromisedJobReturnValue>
                |- <MHC molecule 1>_10_mer.faa: <PromisedJobReturnValue>
                |
                ..
                +- <MHC molecule n>_10_mer.faa: <PromisedJobReturnValue>
         mhcii_preds: Dict of return value from running predictions on a given
                     mhc for all peptides of length 15.
             mhci_preds
                |- <MHC molecule 1>_15_mer.faa: <PromisedJobReturnValue>
                |
                ..
                +- <MHC molecule n>_15_mer.faa: <PromisedJobReturnValue>

    This module corresponds to node 18 on the tree
    """
    job.fileStore.logToMaster('Running spawn_anti on %s' % univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    mhci_options, mhcii_options = mhc_options
    pept_files = {
        '9_mer.faa': transgened_files['transgened_tumor_9_mer_snpeffed.faa'],
        '10_mer.faa': transgened_files['transgened_tumor_10_mer_snpeffed.faa'],
        '15_mer.faa': transgened_files['transgened_tumor_15_mer_snpeffed.faa']}
    input_files = {
        'mhci_alleles.list': phlat_files['mhci_alleles.list'],
        'mhcii_alleles.list': phlat_files['mhcii_alleles.list'],
        'mhci_restrictions.list': mhci_options['method_file'],
        'mhcii_restrictions.list': mhcii_options['method_file']}
    input_files = get_files_from_filestore(job, input_files, work_dir)
    # pept_files = get_files_from_filestore(job, pept_files, work_dir)
    mhci_alleles, mhcii_alleles = [], []
    with open(input_files['mhci_alleles.list'], 'r') as mhci_file:
        for line in mhci_file:
            mhci_alleles.append(line.strip())
    with open(input_files['mhcii_alleles.list'], 'r') as mhcii_file:
        for line in mhcii_file:
            mhcii_alleles.append(line.strip())
    # This file contains the list of allele:pept length combinations supported
    # by each prediction type.
    with open(input_files['mhci_restrictions.list'], 'r') as restfile:
        mhci_restrictions = json.load(restfile)
    with open(input_files['mhcii_restrictions.list'], 'r') as restfile:
        mhcii_restrictions = json.load(restfile)
    # Make a regexp to convert non alphanumeric characters in HLA names to _
    strip_allele_re = re.compile('[^A-Z0-9]')
    # For each mhci allele:peptfile combination, spawn a job and store the job handle in the dict.
    # Then do the same for mhcii
    mhci_preds, mhcii_preds = {}, {}
    for allele in mhci_alleles:
        stripped_allele = re.sub(strip_allele_re, '_', allele)
        for peptfile in ['9_mer.faa', '10_mer.faa']:
            peplen = peptfile.split('_')[0]
            # Ensure that the allele is among the list of accepted alleles
            try:
                if not mhci_restrictions[allele][peplen]:
                    continue
            except KeyError:
                continue
            predfile = ''.join([stripped_allele, '_', peptfile[:-4], '_mer.pred'])
            mhci_preds[predfile] = job.addChildJobFn(predict_mhci_binding, pept_files[peptfile],
                                                     allele, peplen, univ_options,
                                                     mhci_options, disk='10G').rv()
    for allele in mhcii_alleles:
        stripped_allele = re.sub(strip_allele_re, '_', allele)
        predfile = ''.join([stripped_allele, '_15_mer.pred'])
        if allele not in mhcii_restrictions[mhcii_options['pred']]:
            continue
        mhcii_preds[predfile] = job.addChildJobFn(predict_mhcii_binding, pept_files['15_mer.faa'],
                                                  allele, univ_options, mhcii_options,
                                                  disk='10G').rv()
    return mhci_preds, mhcii_preds


def predict_mhci_binding(job, peptfile, allele, peplen, univ_options,
                         mhci_options):
    """
    This module will predict MHC:peptide binding for peptides in the files created in node XX to
    ALLELE.  ALLELE represents an MHCI allele.

    This module corresponds to node 18 on the tree
    """
    job.fileStore.logToMaster('Running mhci on %s:%s:%s' % (univ_options['patient'], allele,
                                                            peplen))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'peptfile.faa': peptfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = [mhci_options['pred'],
                  allele,
                  peplen,
                  input_files['peptfile.faa']]
    with open('/'.join([work_dir, 'predictions.tsv']), 'w') as predfile:
        docker_call(tool='mhci', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=predfile, interactive=True)
    output_file = job.fileStore.writeGlobalFile(predfile.name)
    return output_file


def predict_mhcii_binding(job, peptfile, allele, univ_options, mhcii_options):
    """
    This module will predict MHC:peptide binding for peptides in the files created in node YY to
    ALLELE.  ALLELE represents an MHCII allele.

    The module returns (PREDFILE, PREDICTOR) where PREDFILE contains the predictions and PREDICTOR
    is the predictor used (Consensus, NetMHCIIpan, or Sturniolo).

    This module corresponds to node 19 on the tree
    """
    job.fileStore.logToMaster('Running mhcii on %s:%s' % (univ_options['patient'], allele))
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'peptfile.faa': peptfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = [mhcii_options['pred'],
                  allele,
                  input_files['peptfile.faa']]
    with open('/'.join([work_dir, 'predictions.tsv']), 'w') as predfile:
        docker_call(tool='mhcii', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=predfile, interactive=True)
    run_netMHCIIpan = True
    with open(predfile.name, 'r') as predfile:
        for line in predfile:
            if not line.startswith('HLA'):
                continue
            if line.strip().split('\t')[5] == 'NetMHCIIpan':
                break
            # If the predictor type is sturniolo then it needs to be processed differently
            elif line.strip().split('\t')[5] == 'Sturniolo':
                predictor = 'Sturniolo'
            else:
                predictor = 'Consensus'
            run_netMHCIIpan = False
            break
    if run_netMHCIIpan:
        NetMHCIIpan = job.addChildJobFn(predict_netmhcii_binding, peptfile, allele, univ_options,
                                        disk='10G')
        return NetMHCIIpan.rv()
    else:
        output_file = job.fileStore.writeGlobalFile(predfile.name)
        return output_file, predictor


def predict_netmhcii_binding(job, peptfile, allele, univ_options):
    """
    This module will predict MHC:peptide binding for peptides in the files created in node YY to
    ALLELE.  ALLELE represents an MHCII allele.

    This module corresponds to node 19 on the tree
    """
    job.fileStore.logToMaster('Running netmhciipan on %s' % allele)
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {
        'peptfile.faa': peptfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    # netMHCIIpan accepts differently formatted alleles so we need to modify the input alleles
    if allele.startswith('HLA-DQA') or allele.startswith('HLA-DPA'):
        allele = re.sub(r'[*:]', '', allele)
        allele = re.sub(r'/', '-', allele)
    elif allele.startswith('HLA-DRB'):
        allele = re.sub(r':', '', allele)
        allele = re.sub(r'\*', '_', allele)
        allele = allele.lstrip('HLA-')
    else:
        raise RuntimeError('Unknown allele seen')
    parameters = ['-a', allele,
                  '-xls', '1',
                  '-xlsfile', 'predictions.tsv',
                  '-f', input_files['peptfile.faa']]
    # netMHC writes a lot of useless stuff to sys.stdout so we open /dev/null and dump output there.
    with open(os.devnull, 'w') as output_catcher:
        docker_call(tool='netmhciipan:final', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=output_catcher)
    output_file = job.fileStore.writeGlobalFile('/'.join([work_dir, 'predictions.tsv']))
    return output_file, 'netMHCIIpan'


def merge_mhc_peptide_calls(job, antigen_predictions, transgened_files):
    """
    This module will merge all the calls from nodes 18 and 19, and will filter for the top X%% of
    binders of each allele.  The module will then call the rank boosting script to finish off the
    pipeline.

    This module corresponds to node 19 on the tree
    """
    job.fileStore.logToMaster('Merging MHC calls')
    work_dir = job.fileStore.getLocalTempDir()
    pept_files = {
        '10_mer.faa': transgened_files['transgened_tumor_10_mer_snpeffed.faa'],
        '10_mer.faa.map': transgened_files['transgened_tumor_10_mer_snpeffed.faa.map'],
        '15_mer.faa': transgened_files['transgened_tumor_15_mer_snpeffed.faa'],
        '15_mer.faa.map': transgened_files['transgened_tumor_15_mer_snpeffed.faa.map']}
    mhci_preds, mhcii_preds = antigen_predictions
    mhci_files = get_files_from_filestore(job, mhci_preds, work_dir)
    # First split mhcii_preds into prediction files and predictors and maintain keys so we can later
    # reference them in pairs
    mhcii_predictors = {x: y[1] for x, y in mhcii_preds.items()}
    mhcii_files = {x: y[0] for x, y in mhcii_preds.items()}
    mhcii_files = get_files_from_filestore(job, mhcii_files, work_dir)
    # Get peptide files
    pept_files = get_files_from_filestore(job, pept_files, work_dir)

    # Merge MHCI calls
    # Read 10-mer pepts into memory
    peptides = read_peptide_file(pept_files['10_mer.faa'])
    with open(pept_files['10_mer.faa.map'], 'r') as mapfile:
        pepmap = json.load(mapfile)
    # Incorporate peptide names into the merged calls
    with open('/'.join([work_dir, 'mhci_merged_files.list']), 'w') as mhci_resfile:
        for mhcifile in mhci_files.values():
            with open(mhcifile, 'r') as mf:
                for line in mf:
                    # Skip header lines
                    if not line.startswith('HLA'):
                        continue
                    line = line.strip().split('\t')
                    allele = line[0]
                    pept = line[5]
                    pred = line[7]
                    if float(pred) > 5.00:
                        continue
                    print_mhc_peptide((allele, pept, pred, pept), peptides, pepmap, mhci_resfile)
    # Merge MHCII calls
    # read 15-mer pepts into memory
    peptides = read_peptide_file(pept_files['15_mer.faa'])
    with open(pept_files['15_mer.faa.map'], 'r') as mapfile:
        pepmap = json.load(mapfile)
    # Incorporate peptide names into the merged calls
    with open('/'.join([work_dir, 'mhcii_merged_files.list']), 'w') as \
            mhcii_resfile:
        for mhciifile in mhcii_files.keys():
            core_col = None  # Variable to hold the column number with the core
            if mhcii_predictors[mhciifile] == 'Consensus':
                with open(mhcii_files[mhciifile], 'r') as mf:
                    for line in mf:
                        # Skip header lines
                        if not line.startswith('HLA'):
                            continue
                        line = line.strip().split('\t')
                        allele = line[0]
                        pept = line[4]
                        pred = line[6]
                        if core_col:
                            core = line[core_col] if core_col else 'NOCORE'
                        else:
                            methods = line[5].lstrip('Consensus(').rstrip(')')
                            methods = methods.split(',')
                            if 'NN' in methods:
                                core_col = 13
                            elif 'netMHCIIpan' in methods:
                                core_col = 17
                            elif 'Sturniolo' in methods:
                                core_col = 19
                            elif 'SMM' in methods:
                                core_col = 10
                            core = line[core_col] if core_col else 'NOCORE'
                        if float(pred) > 5.00:
                            continue
                        print_mhc_peptide((allele, pept, pred, core), peptides, pepmap,
                                          mhcii_resfile)
            elif mhcii_predictors[mhciifile] == 'Sturniolo':
                with open(mhcii_files[mhciifile], 'r') as mf:
                    for line in mf:
                        # Skip header lines
                        if not line.startswith('HLA'):
                            continue
                        line = line.strip().split('\t')
                        allele = line[0]
                        pept = line[5]
                        pred = line[6]
                        core = line[19]  #
                        if float(pred) > 5.00:
                            continue
                        print_mhc_peptide((allele, pept, pred, core), peptides, pepmap,
                                          mhcii_resfile)
            elif mhcii_predictors[mhciifile] == 'netMHCIIpan':
                with open(mhcii_files[mhciifile], 'r') as mf:
                    # Get the allele from the first line and skip the second line
                    allele = re.sub('-DQB', '/DQB', mf.readline().strip())
                    _ = mf.readline()
                    for line in mf:
                        line = line.strip().split('\t')
                        pept = line[1]
                        pred = line[5]
                        core = 'NOCORE'
                        peptide_name = line[2]
                        if float(pred) > 5.00:
                            continue
                        print(allele, pept, peptide_name, core, '0', pred, pepmap[peptide_name],
                              sep='\t', file=mhcii_resfile)
            else:
                raise RuntimeError('Shouldn\'t ever see this!!!')
    output_files = defaultdict()
    for mhc_file in [mhci_resfile.name, mhcii_resfile.name]:
        output_files[os.path.split(mhc_file)[1]] = job.fileStore.writeGlobalFile(mhc_file)
    return output_files


def boost_ranks(job, gene_expression, merged_mhc_calls, transgene_out, univ_options,
                rank_boost_options):
    """
    This is the final module in the pipeline.  It will call the rank boosting R
    script.

    This module corresponds to node 21 in the tree
    """
    job.fileStore.logToMaster('Running boost_ranks on %s' % univ_options['patient'])
    work_dir = os.path.join(job.fileStore.getLocalTempDir(), univ_options['patient'])
    os.mkdir(work_dir)
    input_files = {
        'rsem_quant.tsv': gene_expression,
        'mhci_merged_files.tsv': merged_mhc_calls['mhci_merged_files.list'],
        'mhcii_merged_files.tsv': merged_mhc_calls['mhcii_merged_files.list'],
        'mhci_peptides.faa': transgene_out['transgened_tumor_10_mer_snpeffed.faa'],
        'mhcii_peptides.faa': transgene_out['transgened_tumor_15_mer_snpeffed.faa']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    output_files = {}
    for mhc in ('mhci', 'mhcii'):
        parameters = [mhc,
                      input_files[''.join([mhc, '_merged_files.tsv'])],
                      input_files['rsem_quant.tsv'],
                      input_files[''.join([mhc, '_peptides.faa'])],
                      rank_boost_options[''.join([mhc, '_combo'])]
                      ]
        docker_call(tool='rankboost', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'])
        output_files[mhc] = {
            ''.join([mhc, '_concise_results.tsv']):
                job.fileStore.writeGlobalFile(''.join([work_dir, '/', mhc,
                                                       '_merged_files_concise_results.tsv'])),
            ''.join([mhc, '_detailed_results.tsv']):
                job.fileStore.writeGlobalFile(''.join([work_dir, '/', mhc,
                                                       '_merged_files_detailed_results.tsv']))}
    export_results(work_dir, univ_options)
    return output_files


def tool_specific_param_generator(job, config_file):
    """
    This is a generator function to parse and yield the various groups of parameters from
    CONFIG_FILE one at a time.
    Input arguments:
        config_file - a file handle to an open file stream that is reading the
                      input config file to the pipeline.
    Return (Yielded) Values:
        group_name - The name of the group that is being yielded
        group_params - The parameters for the group GROUPNAME
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Initialize the return values.  group_name == None will be used to bypass the first #-prefixed
    # group in the file
    group_params = defaultdict()
    group_name = None
    for line in config_file:
        line = line.strip()
        if line.startswith('##') or len(line) == 0:
            continue
        if line.startswith('#'):
            # Skip first #-prefixed string
            if group_name is None:
                group_name = line.lstrip('#').strip()
                continue
            else:
                yield group_name, group_params
                group_params = defaultdict(int)
                group_name = line.lstrip('#').strip()
                continue
        else:
            line = line.strip().split()
            if len(line) != 2:
                raise ParameterError('Found a problem in the config file while attempting to ' +
                                     'parse %s in group %s' % (line[0], group_name) + '.  Every ' +
                                     'parameter takes ONLY one argument.')
            # If a file is of the type file, vcf, tar or fasta, it needs to be downloaded from S3 if
            # reqd, then written to job store.
            if [x for x in ['file', 'vcf', 'tar', 'fasta', 'fai', 'idx', 'dict'] if x in line[0]]:
                group_params[line[0]] = job.addChildJobFn(get_pipeline_inputs, line[0],
                                                          line[1]).rv()
            else:
                group_params[line[0]] = line[1]
    yield group_name, group_params


def get_pipeline_inputs(job, input_flag, input_file):
    """
    Get the input file from s3 or disk, untargz if necessary and then write to file job store.
    :param job: job
    :param str input_flag: The name of the flag
    :param str input_file: The value passed in the config file
    :return: The jobstore ID for the file
    """
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.logToMaster('Obtaining file (%s) to the file job store' % os.path.basename(
            input_file))
    if input_file.startswith('http'):
        assert input_file.startswith('https://s3'), input_file + ' is not an S3 file'
        input_file = get_file_from_s3(job, input_file, write_to_jobstore=False)
    else:
        assert os.path.exists(input_file), 'Bogus Input : ' + input_file
    # If the file isn't a tarball, then it is a single file that is tar.gzipped for the
    # sake of maximum compression instead of enveloping a folder. Thus it should be
    # decompressed before writing to job store. Also, this is cool but they will by
    # default also be dumped into the cache!
    if 'tar' not in input_flag:
        input_file = untargz(input_file, work_dir)
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
    allowed_samples = {'tumor_dna_fastq_prefix', 'tumor_rna_fastq_prefix',
                       'normal_dna_fastq_prefix'}
    if set(fastqs.keys()).difference(allowed_samples) != {'patient_id'}:
        raise ParameterError('Sample with the following parameters has an error:\n' +
                             '\n'.join(fastqs.values()))
    # For each sample type, check if the prefix is an S3 link or a regular file
    # Download S3 files.
    for sample_type in ['tumor_dna', 'tumor_rna', 'normal_dna']:
        prefix = fastqs[''.join([sample_type, '_fastq_prefix'])]
        # if the file was an xml, process it before moving further.
        if prefix.endswith('.xml'):
            prefix = get_file_from_cghub(job, fastqs[''.join([sample_type, '_fastq_prefix'])],
                                         univ_options['cghub_key'], univ_options,
                                         write_to_jobstore=False)
        # If gzipped, remember that
        if prefix.endswith('.gz'):
            prefix = os.path.splitext(prefix)[0]
            fastqs['gzipped'] = True
        else:
            fastqs['gzipped'] = False
        # Is the file .fastq or .fq
        if prefix.endswith(('.fastq', '.fq')):
            prefix, extn = os.path.splitext(prefix)
            # If the file was gzipped, add that to the extension
            if fastqs['gzipped']:
                extn = ''.join([extn, '.gz'])
        else:
            raise ParameterError('Are the inputs for patient (%s) fastq/fq?' % fastqs['patient_id'])
        # Handle the R1/R2 identifiers in the prefix.  That is added by the program.
        assert prefix.endswith('1'), 'Prefix didn\'t end with 1.<fastq/fq>[.gz]: (%s.%s)' % (prefix,
                                                                                             extn)
        prefix = prefix[:-1]
        # If it is a weblink, it HAS to be from S3
        if prefix.startswith('http'):
            assert prefix.startswith('https://s3'), 'Not an S3 link'
            fastqs[sample_type] = [
                get_file_from_s3(job, ''.join([prefix, '1', extn]), univ_options['sse_key']),
                get_file_from_s3(job, ''.join([prefix, '2', extn]), univ_options['sse_key'])]
        else:
            # Relies heavily on the assumption that the pair will be in the same
            # folder
            assert os.path.exists(''.join([prefix, '1', extn])), 'Bogus input: %s' % ''.join(
                    [prefix, '1', extn])
            # Python lists maintain order hence the values are always guaranteed to be
            # [fastq1, fastq2]
            fastqs[sample_type] = [
                job.fileStore.writeGlobalFile(''.join([prefix, '1', extn])),
                job.fileStore.writeGlobalFile(''.join([prefix, '2', extn]))]
    [fastqs.pop(x) for x in allowed_samples]
    return fastqs


def get_files_from_filestore(job, files, work_dir, cache=True, docker=False):
    """
    This is adapted from John Vivian's return_input_paths from the RNA-Seq pipeline.

    Returns the paths of files from the FileStore if they are not present.
    If docker=True, return the docker path for the file.
    If the file extension is tar.gz, then tar -zxvf it.

    files is a dict with:
        keys = the name of the file to be returned in toil space
        value = the input value for the file (can be toil temp file)
    work_dir is the location where the file should be stored
    cache indiciates whether caching should be used
    """
    for name in files.keys():
        outfile = job.fileStore.readGlobalFile(files[name], '/'.join([work_dir, name]), cache=cache)
        # If the file pointed to a tarball, extract it to WORK_DIR
        if tarfile.is_tarfile(outfile) and file_xext(outfile).startswith('.tar'):
            untar_name = os.path.basename(strip_xext(outfile))
            files[untar_name] = untargz(outfile, work_dir)
            files.pop(name)
            name = os.path.basename(untar_name)
        # If the file is gzipped but NOT a tarfile, gunzip it to work_dir. However, the file is
        # already named x.gz so we need to write to a temporary file x.gz_temp then do a move
        # operation to overwrite x.gz.
        elif is_gzipfile(outfile) and file_xext(outfile) == '.gz':
            ungz_name = strip_xext(outfile)
            with gzip.open(outfile, 'rb') as gz_in, open(ungz_name, 'w') as ungz_out:
                shutil.copyfileobj(gz_in, ungz_out)
            files[os.path.basename(ungz_name)] = outfile
            files.pop(name)
            name = os.path.basename(ungz_name)
        else:
            files[name] = outfile
        # If the files will be sent to docker, we will mount work_dir to the container as /data and
        # we want the /data prefixed path to the file
        if docker:
            files[name] = docker_path(files[name])
    return files


def most_probable_alleles(allele_list):
    """
    This module accepts a list of tuples of (allele, p_value) pairs. It returns the 2 most probable
    alleles for that group.
    """
    all_alleles = defaultdict()
    # First collect all the keys.  Make a dict with allele as key and list of pvalues as value
    for allele, pvalue in allele_list:
        allele = re.split(':', allele)
        # Ensure allele has enough resolution for mhc:peptide prediciton.
        # HLA-A*02:01:04  -> ['HLA-A*02', '01', '04']   => At least 2 fields are required for
        #                                                  satisfying criteria.
        if len(allele) < 2:
            continue
        allele = ':'.join([allele[0], allele[1]])  # stitch back together
        try:
            all_alleles[allele].append(float(pvalue))
        except KeyError:
            all_alleles[allele] = [float(pvalue)]
    # If there are less than 2 alleles, report all
    if len(all_alleles.keys()) <= 2:
        return all_alleles.keys()
    # Else, get the two with most evidence.  Evidence is gauged by
    # a) How many files (of the 3) thought that Allele was present
    # b) In a tie, who has a lower avg p value
    # In the lambda function, if 2 alleles have the same number of calls, the sum of the p values is
    # a measure of the avg because avg = sum / n and n is equal in both of them.
    else:
        return sorted(all_alleles.keys(), key=lambda x: \
            (-len(all_alleles[x]), sum(all_alleles[x])))[0:2]


def read_peptide_file(in_peptfile):
    """
    This module reads an input peptide fasta file into memory in the form of a dict with
    key = fasta record name
    value = corresponding peptide sequence
    """
    peptides = defaultdict()
    with open(in_peptfile, 'r') as peptfile:
        for line in peptfile:
            if line.startswith('>'):
                pept = line.strip().lstrip('>')
                peptides[pept] = ''
            else:
                peptides[pept] = line.strip()
    return peptides


def merge_vcfs(vcf_file, merged_mut_file):
    """
    This module will accept the vcf files for mutect and radia read into memory in a dict object
    VCF_FILE and will merge the calls.  Merged calls are printed to MERGED_MUT_FILE.

    VCF_FILE is a dict with
    key : mutation caller (mutect or radia)
    value : dict with
            key: (chrom, pos, ref, alt)
            value: vcf line in list form (split by tab)
    """
    mutect_keys = set(vcf_file['mutect'].keys())
    radia_keys = set(vcf_file['radia'].keys())
    common_keys = radia_keys.intersection(mutect_keys)
    # Open as append since the header is already written
    with open(merged_mut_file, 'a') as outfile:
        for mutation in common_keys:
            print('\t'.join(vcf_file['radia'][mutation]), file=outfile)
    return None


def parse_radia_multi_alt(infile, outfile):
    """
    This function will parse the vcf to detect sites having multiple alt alleles and pick out on the
    most likely ones.
    INFILE : open file handle for the input vcf
    OUTFILE : open file handle for the output vcf

    The columns in INFILE are
     [0] CHROM
     [1] POS
     [2] ID
     [3] REF
     [4] ALT
     [5] QUAL
     [6] FILTER
     [7] INFO
     [8] FORMAT
     [9] DNA_NORMAL
    [10] DNA_TUMOR
    [11] RNA_TUMOR  -  Not always present
    """
    for line in infile:
        # Print header to putfile
        if line.startswith('#'):
            print(line.strip(), file=outfile)
            continue
        line = line.strip().split('\t')
        # If there is just 1 ALT allele, print and continue
        if len(line[4]) == 1:
            print('\t'.join(line), file=outfile)
        # If not, process
        else:
            seq_field_indeces = [9, 10]
            alleles = [line[3]] + line[4].split(',')  # all alleles, incl. REF
            # collect tumor, normal and (if present) rna AD and AFs
            # AD = Depth of reads supporting each allele
            # AF = Fraction of reads supporting each allele
            normal_AD = line[9].split(':')[5].split(',')
            normal_AF = line[9].split(':')[6].split(',')
            tumor_AD = line[10].split(':')[5].split(',')
            tumor_AF = line[10].split(':')[6].split(',')
            if len(line[11]) > 1:
                rna_AD = line[11].split(':')[5].split(',')
                rna_AF = line[11].split(':')[6].split(',')
                seq_field_indeces += [11]  # append rna since it is present
            else:
                # If rna is missing, set RNA_AD and RNA_AF to null sets for easily integrating into
                # the logic in the following code
                rna_AD = rna_AF = [0, 0, 0, 0]
            # Initialise variables to store the probable ALT alleles and the index values of the
            # same wrt AD and AF
            out_alleles = set([])
            out_AF_AD_index = {0}
            # parse AD and AF to get most probable ALT alleles
            for i in range(1, len(normal_AF)):
                # Criteria for selection = AD > 4 and AF >0.1 in either tumor or RNA, given normal
                # AF < 0.1
                if ((float(tumor_AF[i]) >= 0.1 and int(tumor_AD[i]) >= 4) \
                        or (float(rna_AF[i]) >= 0.1 and int(rna_AD[i]) >= 4)) \
                        and (float(normal_AF[i]) < 0.1):
                    out_alleles.add(alleles[i])
                    out_AF_AD_index.add(i)
            # If the number of probable alleles is greater than 0 the print to outfile with the
            # modified allele fraction representing reads corrresponding to all alleles
            if len(out_alleles) > 0:
                line[4] = ','.join(out_alleles)  # set alt alleles
                # Modify the AD and AF values in the TUMOR/NORMAL/RNA fields
                # one at a time.  Seq fields contain
                # [0] GT* - Genotype
                # [1] DP - Read depth at this position in the sample
                # [2] INDEL - Number of indels
                # [3] START - Number of reads starting at this position
                # [4] STOP - Number of reads stopping at this position
                # [5] AD* - Depth of reads supporting alleles
                # [6] AF* - Fraction of reads supporting alleles
                # [7] BQ* - Avg base quality for reads supporting alleles
                # [8] SB* - Strand Bias for reads supporting alleles
                # Fields marked with *s are teh ones that contain info for each seq field and need
                # to be modified
                for seq_field_index in seq_field_indeces:
                    # Get the details for seq_field
                    deets = line[seq_field_index].split(':')
                    # modify fields 5 thu 8 to hold only info for the probable
                    # alleles
                    for field_index in range(5, 9):
                        field = deets[field_index].split(",")
                        deets[field_index] = ",".join([x for i, x in enumerate(field)
                                                       if i in out_AF_AD_index])
                    # Modify DP to hold the new total of reads
                    deets[1] = str(sum([int(x) for x in deets[5].split(",")]))
                    # get the most likely genotypes based on AD and AF
                    GT_by_AD = set([i for i, x in enumerate(deets[5].split(",")) if int(x) >= 4])
                    GT_by_AF = set([i for i, x in enumerate(deets[6].split(",")) \
                                    if float(x) >= 0.1])
                    # Get the consensus genotype
                    GT = GT_by_AD.intersection(GT_by_AF)
                    if len(GT) == 0:
                        deets[0] = "0/0"
                    elif len(GT) == 1:
                        deets[0] = "/".join([str(x) for x in GT] + [str(x) for x in GT])
                    elif len(GT) == 2:
                        deets[0] = "/".join([str(x) for x in GT])
                    else:
                        print("ERROR : triple genotype detected", file=sys.stderr)
                        print(line, file=sys.stdout)
                    # Rejoin the details line
                    line[seq_field_index] = ":".join(deets)
                # Print the modified line to output
                print("\t".join(line), file=outfile)
            # Else do nothing
            else:
                pass


def print_mhc_peptide(neoepitope_info, peptides, pepmap, outfile):
    """
    To reduce code redundancy, this module will accept data from merge_mhc_peptide_calls for a given
    neoepitope and print it to outfile
    ARGUMENTS
    1. neoepitope_info: Tuple of (<allele>, <peptide_sequence>,
                                  <binding_prediction>)
    2. peptides: Dict of all IARS considered
           peptides
              |- 'neoepitope_1': <peptide_sequence>
              ..
              |- 'neoepitope_n': <peptide_sequence>
    3. pepmap: Info correlating neoepitope with the gene and transcript level
               mutations.
           peptides
              |- 'neoepitope_1':
              |      'ensembl_gene\thugo_gene\tcomma_sep_transcript_mutations'
              ..
              +- 'neoepitope_n':
                     'ensembl_gene\thugo_gene\tcomma_sep_transcript_mutations'

    """
    allele, pept, pred, core = neoepitope_info
    peptide_names = [x for x, y in peptides.items() if pept in y]
    # For each peptide, append the ensembl gene
    for peptide_name in peptide_names:
        print(allele, pept, peptide_name, core, '0', pred, pepmap[peptide_name], sep='\t',
              file=outfile)
    return None


def docker_path(filepath):
    """
    Given a path, returns that files path inside the docker mount directory
    (/data).
    """
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(tool, tool_parameters, work_dir, java_opts=None, outfile=None,
                dockerhub='aarjunrao', interactive=False):
    """
    Makes subprocess call of a command to a docker container. work_dir MUST BE AN ABSOLUTE PATH or
    the call will fail.  outfile is an open file descriptor to a writeable file.
    """
    # If an outifle has been provided, then ensure that it is of type file, it is writeable, and
    # that it is open.
    if outfile:
        assert isinstance(outfile, file), 'outfile was not passsed a file'
        assert outfile.mode in ['w', 'a', 'wb', 'ab'], 'outfile not writeable'
        assert not outfile.closed, 'outfile is closed'
    # If the call is interactive, set intereactive to -i
    if interactive:
        interactive = '-i'
    else:
        interactive = ''
    # If a tag is passed along with the image, use it.
    if ':' in tool:
        docker_tool = '/'.join([dockerhub, tool])
    # Else use 'latest'
    else:
        docker_tool = ''.join([dockerhub, '/', tool, ':latest'])
    # Get the docker image on the worker if needed
    call = ['docker', 'images']
    dimg_rv = subprocess.check_output(call)
    existing_images = [':'.join(x.split()[0:2]) for x in dimg_rv.splitlines()
                       if x.startswith(dockerhub)]
    if docker_tool not in existing_images:
        try:
            call = ' '.join(['docker', 'pull', docker_tool]).split()
            subprocess.check_call(call)
        except subprocess.CalledProcessError as err:
            raise RuntimeError('docker command returned a non-zero exit status ' +
                               '(%s)' % err.returncode + 'for command \"%s\"' % ' '.join(call),)
        except OSError:
            raise RuntimeError('docker not found on system. Install on all' +
                               ' nodes.')
    # If java options have been provided, it needs to be in the docker call
    if java_opts:
        base_docker_call = ' docker run -e JAVA_OPTS=-Xmx{} '.format(java_opts) + '--rm=true ' + \
            '-v {}:/data --log-driver=none '.format(work_dir) + interactive
    else:
        base_docker_call = ' docker run --rm=true -v {}:/data '.format(work_dir) + \
            '--log-driver=none ' + interactive
    call = base_docker_call.split() + [docker_tool] + tool_parameters
    try:
        subprocess.check_call(call, stdout=outfile)
    except subprocess.CalledProcessError as err:
        raise RuntimeError('docker command returned a non-zero exit status (%s)' % err.returncode +
                           'for command \"%s\"' % ' '.join(call),)
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def untargz(input_targz_file, untar_to_dir):
    """
    This module accepts a tar.gz archive and untars it.

    RETURN VALUE: path to the untar-ed directory/file

    NOTE: this module expects the multiple files to be in a directory before
          being tar-ed.
    """
    assert tarfile.is_tarfile(input_targz_file), 'Not a tar file.'
    tarball = tarfile.open(input_targz_file)
    return_value = os.path.join(untar_to_dir, tarball.getmembers()[0].name)
    tarball.extractall(path=untar_to_dir)
    tarball.close()
    return return_value


def is_gzipfile(filename):
    """
    This function attempts to ascertain the gzip status of a file based on the "magic signatures" of
    the file. This was taken from the stack overflow
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress
    """
    assert os.path.exists(filename), 'Input {} does not '.format(filename) + \
        'point to a file.'
    with file(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == '\x1f\x8b\x08':
            # bam files are bgzipped and they share the magic sequence with gzip.  Pysam will error
            # if the input is gzip but not if it is a bam.
            try:
                _ = Samfile(filename)
            except ValueError:
                return True
            else:
                return False
        else:
            return False


def generate_unique_key(master_key, url):
    """
    This module will take a master key and a url, and then make a new key specific to the url, based
    off the master.
    """
    with open(master_key, 'r') as keyfile:
        master_key = keyfile.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
        'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not ' + \
        '32 characters: {}'.format(new_key)
    return new_key


def get_file_from_s3(job, s3_url, encryption_key=None, write_to_jobstore=True):
    """
    Downloads a supplied URL that points to an unencrypted, unprotected file on Amazon S3. The file
    is downloaded and a subsequently written to the jobstore and the return value is a the path to
    the file in the jobstore.
    """
    work_dir = job.fileStore.getLocalTempDir()
    filename = '/'.join([work_dir, os.path.basename(s3_url)])
    # This is common to encrypted and unencrypted downloads
    download_call = ['curl', '-fs', '--retry', '5']
    # If an encryption key was provided, use it to create teh headers that need to be injected into
    # the curl script and append to the call
    if encryption_key:
        key = generate_unique_key(encryption_key, s3_url)
        encoded_key = base64.b64encode(key)
        encoded_key_md5 = base64.b64encode( hashlib.md5(key).digest() )
        h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
        h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
        h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
        download_call.extend(['-H', h1, '-H', h2, '-H', h3])
    # This is also common to both types of downloads
    download_call.extend([s3_url, '-o', filename])
    try:
        subprocess.check_call(download_call)
    except subprocess.CalledProcessError:
        raise RuntimeError('Curl returned a non-zero exit status processing %s. Do you' % s3_url +
                           'have premssions to access the file?')
    except OSError:
        raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(filename)
    if write_to_jobstore:
        filename = job.fileStore.writeGlobalFile(filename)
    return filename


def get_file_from_cghub(job, cghub_xml, cghub_key, univ_options, write_to_jobstore=True):
    """
    This function will download the file from cghub using the xml specified by cghub_xml

    ARGUMENTS
    1. cghub_xml: Path to an xml file for cghub.
    2. cghub_key: Credentials for a cghub download operation.
    3. write_to_jobstore: Flag indicating whether the final product should be written to jobStore.

    RETURN VALUES
    1. A path to the prefix for the fastqs that is compatible with the pipeline.
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Get from S3 if required
    if cghub_xml.startswith('http'):
        assert cghub_xml.startswith('https://s3'), 'Not an S3 link'
        cghub_xml = get_file_from_s3(job, cghub_xml, encryption_key=univ_options['sse_key'],
                                     write_to_jobstore=False)
    else:
        assert os.path.exists(cghub_xml), 'Could not find file: %s' % cghub_xml
    shutil.copy(cghub_xml, os.path.join(work_dir, 'cghub.xml'))
    cghub_xml = os.path.join(work_dir, 'cghub.xml')
    assert os.path.exists(cghub_key), 'Could not find file: %s' % cghub_key
    shutil.copy(cghub_key, os.path.join(work_dir, 'cghub.key'))
    cghub_key = os.path.join(work_dir, 'cghub.key')
    temp_fastqdir = os.path.join(work_dir, 'temp_fastqdir')
    os.mkdir(temp_fastqdir)
    base_parameters = ['-d',  docker_path(cghub_xml),
                  '-c', docker_path(cghub_key),
                  '-p', docker_path(temp_fastqdir)]
    attemptNumber = 0
    while True:
        # timeout increases by 10 mins per try
        parameters = base_parameters + ['-k', str((attemptNumber + 1) * 10)]
        try:
            docker_call('genetorrent', tool_parameters=parameters, work_dir=work_dir,
                        dockerhub=univ_options['dockerhub'])
        except RuntimeError as err:
            time.sleep(600)
            job.fileStore.logToMaster(err.message)
            attemptNumber += 1
            if attemptNumber == 3:
                raise
            else:
                continue
        else:
            break
    analysis_id = [x for x in os.listdir(temp_fastqdir)
                   if not (x.startswith('.') or x.endswith('.gto'))][0]
    files = [x for x in os.listdir(os.path.join(temp_fastqdir, analysis_id))
             if not x.startswith('.')]
    if len(files) == 2:
        prefixes = [os.path.splitext(x)[1] for x in files]
        if {'.bam', '.bai'} - set(prefixes):
            raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are ' +
                               'sure it is, email aarao@ucsc.edu with details.')
        else:
            bamfile = os.path.join(temp_fastqdir, analysis_id,
                                   [x for x in files if x.endswith('.bam')][0])
            return bam2fastq(job, bamfile, univ_options)
    elif len(files) == 1:
        if not files[0].endswith('.tar.gz'):
            raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are ' +
                               'sure it is, email aarao@ucsc.edu with details.')
        else:
            outFastqDir = os.path.join(work_dir, 'fastqs')
            os.mkdir(outFastqDir)
            fastq_file = untargz(os.path.join(temp_fastqdir, analysis_id, files[0]), outFastqDir)
            if fastq_file.endswith(('.fastq', '.fastq.gz')):
                return re.sub('_2.fastq', '_1.fastq', fastq_file)
            else:
                raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you ' +
                                   'are sure it is, email aarao@ucsc.edu with details.')
    else:
        raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are sure ' +
                           'it is, email aarao@ucsc.edu with details.')


def bam2fastq(job, bamfile, univ_options):
    """
    split an input bam to paired fastqs.

    ARGUMENTS
    1. bamfile: Path to a bam file
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                |- 'dockerhub': <dockerhub to use>
                +- 'java_Xmx': value for max heap passed to java
    """
    work_dir = os.path.split(bamfile)[0]
    base_name = os.path.split(os.path.splitext(bamfile)[0])[1]
    parameters = ['SamToFastq',
                  ''.join(['I=', docker_path(bamfile)]),
                  ''.join(['F=/data/', base_name, '_1.fastq']),
                  ''.join(['F2=/data/', base_name, '_2.fastq']),
                  ''.join(['FU=/data/', base_name, '_UP.fastq'])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=univ_options['java_Xmx'])
    first_fastq = ''.join([work_dir, '/', base_name, '_1.fastq'])
    assert os.path.exists(first_fastq)
    return first_fastq


def export_results(file_path, univ_options):
    """
    Write out a file to a given location. The location can be either a directory on the local
    machine, or a folder with a bucket on AWS.
    TODO: Azure support
    :param file_path: The path to the file that neeeds to be transferred to the new location.
    :param univ_options: A dict of the universal options passed to this script. The important dict
                         entries are ouput_folder and storage_location.
                         * storage_location: 'Local' or an 'aws:<bucket_name>'.
                         * output_folder: The folder to store the file. This must exist on the local
                           machine if storage_location is 'Local'. If the storage_location is an aws
                           bucket,  this string represents the path to the file in the bucket.  To
                           keep it in the base directory, specify 'NA'.

    :return: None
    """
    try:
        assert univ_options['output_folder'], 'Need a path to a folder to write out files'
        assert univ_options['storage_location'], 'Need to know where the files need to go. ' + \
                                                 'Local or AWS/Azure, etc.'
    except AssertionError as err:
        # This isn't a game killer.  Continue the pipeline without erroring out but do inform the
        # user about it.
        print('ERROR:', err.message, file=sys.stderr)
        return
    assert os.path.exists(file_path), "Can't copy a file that doesn't exist!"
    if univ_options['output_folder'] == 'NA':
        if univ_options['storage_location'].lower == 'local':
            print('ERROR: Cannot have NA as output folder if storage location is Local',
                  file=sys.stderr)
            return
        output_folder = ''
    else:
        output_folder = univ_options['output_folder']
    # Handle Local
    if univ_options['storage_location'].lower() == 'local':
        # Create the directory if required
        try:
            os.makedirs(univ_options['output_folder'], 755)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise
        output_file = os.path.join(output_folder, os.path.basename(file_path))
        shutil.copy(file_path, output_file)
    # Handle AWS
    elif univ_options['storage_location'].startswith('aws'):
        bucket_name = univ_options['storage_location'].split(':')[-1]
        write_to_s3(file_path, univ_options['sse_key'], bucket_name, output_folder)
    # Can't do Azure or google yet.
    else:
        print("Currently doesn't support anything but Local and aws.")


def file_xext(filepath):
    """
    Get the file extension wrt compression from the filename (is it tar or targz)
    :param str filepath: Path to the file
    :return str ext: Compression extension name
    """
    ext = os.path.splitext(filepath)[1]
    if ext == '.gz':
        xext = os.path.splitext(os.path.splitext(filepath)[0])[1]
        if xext == '.tar':
            ext = xext + ext
    elif ext == '.tar':
        pass # ext is already .tar
    else:
        ext = ''
    return ext


def strip_xext(filepath):
    """
    Strips the compression extension from the filename
    :param filepath: Path to compressed file.
    :return str filepath: Path to the file with the compression extension stripped off.
    """
    ext_size = len(file_xext(filepath).split('.')) - 1
    for i in xrange(0, ext_size):
        filepath = os.path.splitext(filepath)[0]
    return filepath


# Exception for bad parameters provided
class ParameterError(Exception):
    """
    This Error Class will be raised  in the case of a bad parameter provided.
    """
    pass


def main():
    """
    This is the main function for the UCSC Precision Immuno pipeline.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file', dest='config_file', help='Config file to be used in the' +
                        'run.', type=str, required=True, default=None)
    Job.Runner.addToilOptions(parser)
    params = parser.parse_args()
    START = Job.wrapJobFn(parse_config_file, params.config_file).encapsulate()
    Job.Runner.startToil(START, params)
    return None


if __name__ == '__main__':
    main()
