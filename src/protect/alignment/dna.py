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
from protect.mutation_calling.common import sample_chromosomes
from toil.job import PromisedRequirement

import os


# disk for bwa-related functions
def bwa_disk(dna_fastqs, bwa_index):
    return int(6 * ceil(sum([f.size for f in dna_fastqs]) + 524288) +
               2.5 * ceil(bwa_index.size + 524288))


def sam2bam_disk(samfile):
    return int(ceil(1.5 * samfile.size + 524288))


def reheader_disk(bamfile):
    return int(ceil(2 * bamfile.size + 524288))


def regroup_disk(reheader_bam):
    return int(ceil(4 * reheader_bam.size + 524288))


def mkdup_disk(regroup_bam):
    return int(ceil(4 * regroup_bam.size + 524288))


# disk for fixing a GDC bam
def fix_gdc_bam_disk(bamfile):
    return int(2.5 * ceil(bamfile[0].size + 524288))


def align_dna(job, fastqs, sample_type, univ_options, bwa_options):
    """
    A wrapper for the entire dna alignment subgraph.

    :param list fastqs: The input fastqs for alignment
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict bwa_options: Options specific to bwa
    :return: Dict containing output bam and bai
             output_files:
                 |- '<sample_type>_fix_pg_sorted.bam': fsID
                 +- '<sample_type>_fix_pg_sorted.bam.bai': fsID
    :rtype: dict
    """
    bwa = job.wrapJobFn(run_bwa, fastqs, sample_type, univ_options, bwa_options,
                        disk=PromisedRequirement(bwa_disk, fastqs, bwa_options['index']),
                        cores=bwa_options['n'])
    sam2bam = job.wrapJobFn(bam_conversion, bwa.rv(), sample_type, univ_options,
                            bwa_options['samtools'],
                            disk=PromisedRequirement(sam2bam_disk, bwa.rv()))
    # reheader takes the same disk as sam2bam so we can serialize this on the same worker.
    reheader = job.wrapJobFn(fix_bam_header, sam2bam.rv(), sample_type, univ_options,
                             bwa_options['samtools'],
                             disk=PromisedRequirement(sam2bam_disk, bwa.rv()))
    regroup = job.wrapJobFn(add_readgroups, reheader.rv(), sample_type, univ_options,
                            bwa_options['picard'],
                            disk=PromisedRequirement(regroup_disk, reheader.rv()))
    mkdup = job.wrapJobFn(mark_duplicates, regroup.rv(), sample_type, univ_options,
                          bwa_options['picard'],
                          disk=PromisedRequirement(mkdup_disk, regroup.rv()))
    index = job.wrapJobFn(index_bamfile, mkdup.rv(), sample_type, univ_options,
                          bwa_options['samtools'], sample_info='fix_pg_sorted',
                          disk=PromisedRequirement(index_disk, mkdup.rv()))
    job.addChild(bwa)
    bwa.addChild(sam2bam)
    sam2bam.addChild(reheader)
    reheader.addChild(regroup)
    regroup.addChild(mkdup)
    mkdup.addChild(index)
    return index.rv()


def run_bwa(job, fastqs, sample_type, univ_options, bwa_options):
    """
    Align a pair of fastqs with bwa.

    :param list fastqs: The input fastqs for alignment
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict bwa_options: Options specific to bwa
    :return: fsID for the generated sam
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        'dna_1.fastq': fastqs[0],
        'dna_2.fastq': fastqs[1],
        'bwa_index.tar.gz': bwa_options['index']}
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
                  '/'.join([input_files['bwa_index'], univ_options['ref']]),
                  input_files['dna_1.fastq' + gz],
                  input_files['dna_2.fastq' + gz]]
    with open(''.join([work_dir, '/', sample_type, '.sam']), 'w') as samfile:
        docker_call(tool='bwa', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=samfile,
                    tool_version=bwa_options['version'])
    # samfile.name retains the path info
    output_file = job.fileStore.writeGlobalFile(samfile.name)
    job.fileStore.logToMaster('Ran bwa on %s:%s successfully'
                              % (univ_options['patient'], sample_type))
    return output_file


def bam_conversion(job, samfile, sample_type, univ_options, samtools_options):
    """
    Convert a sam to a bam.

    :param dict samfile: The input sam file
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict samtools_options: Options specific to samtools
    :return: fsID for the generated bam
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        sample_type + '.sam': samfile}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    bamfile = '/'.join([work_dir, sample_type + '.bam'])
    parameters = ['view',
                  '-bS',
                  '-o', docker_path(bamfile),
                  input_files[sample_type + '.sam']
                  ]
    docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=samtools_options['version'])
    output_file = job.fileStore.writeGlobalFile(bamfile)
    # The samfile is no longer useful so delete it
    job.fileStore.deleteGlobalFile(samfile)
    job.fileStore.logToMaster('Ran sam2bam on %s:%s successfully'
                              % (univ_options['patient'], sample_type))
    return output_file


def fix_bam_header(job, bamfile, sample_type, univ_options, samtools_options, retained_chroms=None):
    """
    Fix the bam header to remove the command line call.  Failing to do this causes Picard to reject
    the bam.

    :param dict bamfile: The input bam file
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict samtools_options: Options specific to samtools
    :param list retained_chroms: A list of chromosomes to retain
    :return: fsID for the output bam
    :rtype: toil.fileStore.FileID
    """
    if retained_chroms is None:
        retained_chroms = []

    work_dir = os.getcwd()
    input_files = {
        sample_type + '.bam': bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['view',
                  '-H',
                  input_files[sample_type + '.bam']]
    with open('/'.join([work_dir, sample_type + '_input_bam.header']), 'w') as headerfile:
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=headerfile,
                    tool_version=samtools_options['version'])
    with open(headerfile.name, 'r') as headerfile, \
            open('/'.join([work_dir, sample_type + '_output_bam.header']), 'w') as outheaderfile:
        for line in headerfile:
            if line.startswith('@PG'):
                line = '\t'.join([x for x in line.strip().split('\t') if not x.startswith('CL')])
            if retained_chroms and line.startswith('@SQ'):
                if line.strip().split()[1].lstrip('SN:') not in retained_chroms:
                    continue
            print(line.strip(), file=outheaderfile)
    parameters = ['reheader',
                  docker_path(outheaderfile.name),
                  input_files[sample_type + '.bam']]
    with open('/'.join([work_dir, sample_type + '_fixPG.bam']), 'w') as fixpg_bamfile:
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=fixpg_bamfile,
                    tool_version=samtools_options['version'])
    output_file = job.fileStore.writeGlobalFile(fixpg_bamfile.name)
    # The old bam file is now useless.
    job.fileStore.deleteGlobalFile(bamfile)
    job.fileStore.logToMaster('Ran reheader on %s:%s successfully'
                              % (univ_options['patient'], sample_type))
    return output_file


def add_readgroups(job, bamfile, sample_type, univ_options, picard_options):
    """
    Add read groups to the bam.

    :param dict bamfile: The input bam file
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict picard_options: Options specific to picard
    :return: fsID for the output bam
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        sample_type + '.bam': bamfile}
    get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['AddOrReplaceReadGroups',
                  'CREATE_INDEX=false',
                  'I=/data/' + sample_type + '.bam',
                  'O=/data/' + sample_type + '_reheader.bam',
                  'SO=coordinate',
                  'ID=1',
                  ''.join(['LB=', univ_options['patient']]),
                  'PL=ILLUMINA',
                  'PU=12345',
                  ''.join(['SM=', sample_type.rstrip('_dna')])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_xmx=univ_options['java_Xmx'],
                tool_version=picard_options['version'])
    output_file = job.fileStore.writeGlobalFile(
        '/'.join([work_dir, sample_type + '_reheader.bam']))
    # Delete the old bam file
    job.fileStore.deleteGlobalFile(bamfile)
    job.fileStore.logToMaster('Ran add_read_groups on %s:%s successfully'
                              % (univ_options['patient'], sample_type))
    return output_file


def mark_duplicates(job, bamfile, sample_type, univ_options, picard_options):
    """
    Mark duplicates within the bam.

    :param dict bamfile: The input bam file
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict picard_options: Options specific to picard
    :return: fsID for the output bam
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        sample_type + '.bam': bamfile}
    get_files_from_filestore(job, input_files, work_dir, docker=True)

    parameters = ['MarkDuplicates',
                  'I=/data/' + sample_type + '.bam',
                  'O=/data/' + sample_type + '_mkdup.bam',
                  'M=/data/' + sample_type + '_mkdup.metrics',
                  'AS=true',
                  'CREATE_INDEX=true']

    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_xmx=univ_options['java_Xmx'],
                tool_version=picard_options['version'])
    output_file = job.fileStore.writeGlobalFile(
        '/'.join([work_dir, sample_type + '_mkdup.bam']))
    # Delete the old bam file
    job.fileStore.deleteGlobalFile(bamfile)
    job.fileStore.logToMaster('Ran mark_duplicates on %s:%s successfully'
                              % (univ_options['patient'], sample_type))
    return output_file
