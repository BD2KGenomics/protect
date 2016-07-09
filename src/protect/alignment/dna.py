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
from __future__ import absolute_import, print_function
from math import ceil
from protect.alignment.common import index_bamfile, index_disk
from protect.common import docker_call, docker_path, get_files_from_filestore, is_gzipfile, untargz
from toil.job import PromisedRequirement

import os


# disk for bwa-related functions
def bwa_disk(dna_fastqs, bwa_index):
    return 6 * ceil(sum([f.size for f in dna_fastqs]) + 524288) + 2 * ceil(bwa_index.size + 524288)


def sam2bam_disk(samfile):
    return ceil(1.5 * samfile.size + 524288)


def reheader_disk(bamfile):
    return ceil(2 * bamfile.size + 524288)


def regroup_disk(reheader_bam):
    return ceil(4 * reheader_bam.size + 524288)


def align_dna(job, fastqs, sample_type, univ_options, bwa_options):
    """
    This is a convenience function that runs the entire dna alignment subgraph
    """
    bwa = job.wrapJobFn(run_bwa, fastqs, sample_type, univ_options, bwa_options,
                        disk=PromisedRequirement(bwa_disk, fastqs, bwa_options['tool_index']))
    sam2bam = job.wrapJobFn(bam_conversion, bwa.rv(), sample_type, univ_options,
                            disk=PromisedRequirement(sam2bam_disk, bwa.rv()))
    # reheader takes the same disk as sam2bam so we can serialize this on the same worker.
    reheader = job.wrapJobFn(fix_bam_header, sam2bam.rv(), sample_type, univ_options,
                             disk=PromisedRequirement(sam2bam_disk, bwa.rv()))
    regroup = job.wrapJobFn(add_readgroups, reheader.rv(), sample_type, univ_options,
                            disk=PromisedRequirement(regroup_disk, reheader.rv()))
    index = job.wrapJobFn(index_bamfile, regroup.rv(), sample_type, univ_options,
                          disk=PromisedRequirement(index_disk, regroup.rv()))
    job.addChild(bwa)
    bwa.addChild(sam2bam)
    sam2bam.addChild(reheader)
    reheader.addChild(regroup)
    regroup.addChild(index)
    return index.rv()


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
              |- 'tool_index': <JSid for the bwa index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_files: Dict of aligned bam + reference (nested return)
         output_files
             |- '<ST>_fix_pg_sorted.bam': <JSid>
             +- '<ST>_fix_pg_sorted.bam.bai': <JSid>

    This module corresponds to nodes 3 and 4 on the tree
    """
    job.fileStore.logToMaster('Running bwa on %s:%s' % (univ_options['patient'], sample_type))
    work_dir = os.getcwd()
    input_files = {
        'dna_1.fastq': fastqs[0],
        'dna_2.fastq': fastqs[1],
        'bwa_index.tar.gz': bwa_options['tool_index']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    # Handle gzipped file
    gz = '.gz' if is_gzipfile(input_files['dna_1.fastq']) else ''
    if gz:
        for read_file in 'dna_1.fastq', 'dna_2.fastq':
            os.symlink(read_file, read_file + gz)
            input_files[read_file + gz] = input_files[read_file] + gz
    # Untar the index
    input_files['bwa_index'] = untargz(input_files['bwa_index.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['mem',
                  '-t', str(bwa_options['n']),
                  '-v', '1',  # Don't print INFO messages to the stderr
                  '/'.join([input_files['bwa_index'], 'hg19']),
                  input_files['dna_1.fastq' + gz],
                  input_files['dna_2.fastq' + gz]]
    with open(''.join([work_dir, '/', sample_type, '_aligned.sam']), 'w') as samfile:
        docker_call(tool='bwa', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=samfile)
    # samfile.name retains the path info
    output_file = job.fileStore.writeGlobalFile(samfile.name)
    return output_file


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
    work_dir = os.getcwd()
    input_files = {
        sample_type + '_aligned.sam': samfile}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    bamfile = '/'.join([work_dir, sample_type + '_aligned.bam'])
    parameters = ['view',
                  '-bS',
                  '-o', docker_path(bamfile),
                  input_files[sample_type + '_aligned.sam']
                  ]
    docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_file = job.fileStore.writeGlobalFile(bamfile)
    # The samfile is no longer useful so delete it
    job.fileStore.deleteGlobalFile(samfile)
    return output_file


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
    work_dir = os.getcwd()
    input_files = {
        sample_type + '_aligned.bam': bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['view',
                  '-H',
                  input_files[sample_type + '_aligned.bam']]
    with open('/'.join([work_dir, sample_type + '_aligned_bam.header']), 'w') as headerfile:
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=headerfile)
    with open(headerfile.name, 'r') as headerfile, \
            open('/'.join([work_dir, sample_type + '_output_bam.header']), 'w') as outheaderfile:
        for line in headerfile:
            if line.startswith('@PG'):
                line = '\t'.join([x for x in line.strip().split('\t') if not x.startswith('CL')])
            print(line.strip(), file=outheaderfile)
    parameters = ['reheader',
                  docker_path(outheaderfile.name),
                  input_files[sample_type + '_aligned.bam']]
    with open('/'.join([work_dir, sample_type + '_aligned_fixPG.bam']), 'w') as fixpg_bamfile:
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=fixpg_bamfile)
    output_file = job.fileStore.writeGlobalFile(fixpg_bamfile.name)
    # The old bam file is now useless.
    job.fileStore.deleteGlobalFile(bamfile)
    return output_file


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
    work_dir = os.getcwd()
    input_files = {
        sample_type + '_aligned_fixpg.bam': bamfile}
    get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['AddOrReplaceReadGroups',
                  'CREATE_INDEX=false',
                  'I=/data/' + sample_type + '_aligned_fixpg.bam',
                  'O=/data/' + sample_type + '_aligned_fixpg_sorted_reheader.bam',
                  'SO=coordinate',
                  'ID=1',
                  ''.join(['LB=', univ_options['patient']]),
                  'PL=ILLUMINA',
                  'PU=12345',
                  ''.join(['SM=', sample_type.rstrip('_dna')])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=univ_options['java_Xmx'])
    output_file = job.fileStore.writeGlobalFile(
        '/'.join([work_dir, sample_type + '_aligned_fixpg_sorted_reheader.bam']))
    # Delete the old bam file
    job.fileStore.deleteGlobalFile(bamfile)
    return output_file
