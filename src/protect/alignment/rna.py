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
from protect.common import (docker_call,
                            docker_path,
                            export_results,
                            get_files_from_filestore,
                            is_gzipfile,
                            untargz)
from toil.job import PromisedRequirement

import os


# disk for star
def star_disk(rna_fastqs, star_tar):
    return int(4 * ceil(sum([f.size for f in rna_fastqs]) + 524288) +
               2 * ceil(star_tar.size + 524288) +
               5242880)


def align_rna(job, fastqs, univ_options, star_options):
    """
    A wrapper for the entire rna alignment subgraph.

    :param list fastqs: The input fastqs for alignment
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict star_options: Options specific to star
    :return: Dict containing input bam and the generated index (.bam.bai)
    :rtype: dict
    """
    star = job.wrapJobFn(run_star, fastqs, univ_options, star_options,
                         cores=star_options['n'],
                         memory=PromisedRequirement(lambda x: int(1.85 * x.size),
                                                    star_options['index']),
                         disk=PromisedRequirement(star_disk, fastqs, star_options['index']))
    index = job.wrapJobFn(index_star, star.rv(), univ_options, star_options,
                          disk=PromisedRequirement(star_disk, fastqs, star_options['index']))
    job.addChild(star)
    star.addChild(index)
    return index.rv()


def run_star(job, fastqs, univ_options, star_options):
    """
    Align a pair of fastqs with STAR.

    :param list fastqs: The input fastqs for alignment
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict star_options: Options specific to star
    :return: Dict containing output genome bam, genome bai, and transcriptome bam
                 output_files:
                    |- 'rnaAligned.toTranscriptome.out.bam': fsID
                    +- 'rnaAligned.sortedByCoord.out.bam':
                                        +- 'rna_fix_pg_sorted.bam': fsID
    :rtype: dict
    """
    assert star_options['type'] in ('star', 'starlong')
    job.fileStore.logToMaster('Running STAR on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'rna_cutadapt_1.fastq': fastqs[0],
        'rna_cutadapt_2.fastq': fastqs[1],
        'star_index.tar.gz': star_options['index']}
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
        docker_call(tool='star', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], tool_version=star_options['version'])
    else:
        docker_call(tool='starlong', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], tool_version=star_options['version'])
    output_files = defaultdict()
    for bam_file in ['rnaAligned.toTranscriptome.out.bam',
                     'rnaAligned.sortedByCoord.out.bam']:
        output_files[bam_file] = job.fileStore.writeGlobalFile('/'.join([
            work_dir, bam_file]))
    export_results(job, output_files['rnaAligned.toTranscriptome.out.bam'], 'rna_transcriptome.bam',
                   univ_options, subfolder='alignments')
    return output_files


def index_star(job, star_bams, univ_options, star_options):
    """
    A wrapper for indexing the genomic star bam generated by run_star. It is required since run_star
    returns a dict of 2 bams

    :param dict star_bams: The bams from run_star
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict star_options: Options specific to star
    :return: Dict containing input bam and the generated index (.bam.bai)
                     output_files:
                        |- 'rnaAligned.toTranscriptome.out.bam': fsID
                        +- 'rnaAligned.sortedByCoord.out.bam':
                                        |- 'rna_fix_pg_sorted.bam': fsID
                                        +- 'rna_fix_pg_sorted.bam.bai': fsID

    :rtype: dict
    """
    index = job.wrapJobFn(index_bamfile, star_bams['rnaAligned.sortedByCoord.out.bam'], 'rna',
                          univ_options, samtools_options=star_options['samtools'],
                          disk=PromisedRequirement(
                              index_disk, star_bams['rnaAligned.sortedByCoord.out.bam']))
    job.addChild(index)
    star_bams['rnaAligned.sortedByCoord.out.bam'] = index.rv()
    return star_bams
