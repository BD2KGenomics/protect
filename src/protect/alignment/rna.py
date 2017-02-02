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
from collections import defaultdict
from math import ceil
from protect.alignment.common import index_bamfile, index_disk
from protect.common import docker_call, get_files_from_filestore, is_gzipfile, untargz, docker_path
from toil.job import PromisedRequirement

import os


# disk for star
def star_disk(rna_fastqs, star_tar):
    return int(4 * ceil(sum([f.size for f in rna_fastqs]) + 524288) +
               2 * ceil(star_tar.size + 524288) +
               5242880)


def align_rna(job, fastqs, univ_options, star_options):
    """
    This is a convenience function that runs the entire rna alignment subgraph
    """
    star = job.wrapJobFn(run_star, fastqs, univ_options, star_options,
                         cores=star_options['n'],
                         memory=PromisedRequirement(lambda x: int(1.85 * x.size),
                                                    star_options['tool_index']),
                         disk=PromisedRequirement(star_disk, fastqs, star_options['tool_index']))
    index = job.wrapJobFn(index_star, star.rv(), univ_options,
                          disk=PromisedRequirement(star_disk, fastqs, star_options['tool_index']))
    job.addChild(star)
    star.addChild(index)
    return index.rv()


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
             |- 'tool_index': <JSid for the STAR index tarball>
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
    job.fileStore.logToMaster('Running STAR on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'rna_cutadapt_1.fastq': fastqs[0],
        'rna_cutadapt_2.fastq': fastqs[1],
        'star_index.tar.gz': star_options['tool_index']}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=False)
    # Handle gzipped file
    gz = '.gz' if is_gzipfile(input_files['rna_cutadapt_1.fastq']) else ''
    if gz:
        for read_file in 'rna_cutadapt_1.fastq', 'rna_cutadapt_2.fastq':
            os.symlink(read_file, read_file + gz)
            input_files[read_file + gz] = input_files[read_file] + gz
    # Untar the index
    input_files['star_index'] = untargz(input_files['star_index.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['--runThreadN', str(star_options['n']),
                  '--genomeDir', input_files['star_index'],
                  '--outFileNamePrefix', 'rna',
                  '--readFilesIn',
                  input_files['rna_cutadapt_1.fastq' + gz],
                  input_files['rna_cutadapt_2.fastq' + gz],
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  '--quantMode', 'TranscriptomeSAM',
                  '--limitBAMsortRAM', str(job.memory)]
    if gz:
        parameters.extend(['--readFilesCommand', 'zcat'])
    if star_options['type'] == 'star':
        docker_call(tool='star:2.4.2a', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'])
    else:
        docker_call(tool='starlong:2.4.2a', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'])
    output_files = defaultdict()
    for bam_file in ['rnaAligned.toTranscriptome.out.bam',
                     'rnaAligned.sortedByCoord.out.bam']:
        output_files[bam_file] = job.fileStore.writeGlobalFile('/'.join([
            work_dir, bam_file]))
    return output_files


def index_star(job, star_bams, univ_options):
    """
    This is a wrapper functiion for index_bamfile in protect.common which is required since run_star
    returns a dict of 2 bams
    """
    index = job.wrapJobFn(index_bamfile, star_bams['rnaAligned.sortedByCoord.out.bam'], 'rna',
                          univ_options,
                          disk=PromisedRequirement(
                              index_disk, star_bams['rnaAligned.sortedByCoord.out.bam']))
    job.addChild(index)
    star_bams['rnaAligned.sortedByCoord.out.bam'] = index.rv()
    return star_bams
