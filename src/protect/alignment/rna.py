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

from protect.alignment.common import index_bamfile, index_disk, sort_bamfile, sort_disk
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
    return int(5 * ceil(sum([f.size for f in rna_fastqs]) + 524288) +
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
    s_and_i = job.wrapJobFn(sort_and_index_star, star.rv(), univ_options,
                            star_options).encapsulate()
    job.addChild(star)
    star.addChild(s_and_i)
    return s_and_i.rv()


def run_star(job, fastqs, univ_options, star_options):
    """
    Align a pair of fastqs with STAR.

    :param list fastqs: The input fastqs for alignment
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict star_options: Options specific to star
    :return: Dict containing output genome bam, genome bai, and transcriptome bam
                 output_files:
                    |- 'rnaAligned.toTranscriptome.out.bam': fsID
                    +- 'rnaAligned.out.bam': fsID
                    +- 'rnaChimeric.out.junction': fsID
    :rtype: dict
    """
    assert star_options['type'] in ('star', 'starlong')
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

    # Check to see if user is using a STAR-Fusion index
    star_fusion_idx = os.path.join(input_files['star_index'], 'ref_genome.fa.star.idx')
    if os.path.exists(star_fusion_idx):
        input_files['star_index'] = star_fusion_idx

    input_files = {key: docker_path(path, work_dir=work_dir) for key, path in input_files.items()}

    # Using recommended STAR-Fusion parameters:
    # https://github.com/STAR-Fusion/STAR-Fusion/wiki
    parameters = ['--runThreadN', str(star_options['n']),
                  '--genomeDir', input_files['star_index'],
                  '--twopassMode', 'Basic',
                  '--outReadsUnmapped', 'None',
                  '--chimSegmentMin', '12',
                  '--chimJunctionOverhangMin', '12',
                  '--alignSJDBoverhangMin', '10',
                  '--alignMatesGapMax', '200000',
                  '--alignIntronMax', '200000',
                  '--chimSegmentReadGapMax', 'parameter', '3',
                  '--alignSJstitchMismatchNmax', '5', '-1', '5', '5',
                  '--outFileNamePrefix', 'rna',
                  '--readFilesIn',
                  input_files['rna_cutadapt_1.fastq' + gz],
                  input_files['rna_cutadapt_2.fastq' + gz],
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outSAMtype', 'BAM', 'Unsorted',
                  '--quantMode', 'TranscriptomeSAM']
    if gz:
        parameters.extend(['--readFilesCommand', 'zcat'])

    if star_options['type'] == 'star':
        docker_call(tool='star', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], tool_version=star_options['version'])
    else:
        docker_call(tool='starlong', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], tool_version=star_options['version'])
    output_files = defaultdict()
    for output_file in ['rnaAligned.toTranscriptome.out.bam',
                     'rnaAligned.out.bam',
                     'rnaChimeric.out.junction']:
        output_files[output_file] = job.fileStore.writeGlobalFile('/'.join([work_dir, output_file]))
    export_results(job, output_files['rnaAligned.toTranscriptome.out.bam'], 'rna_transcriptome.bam',
                   univ_options, subfolder='alignments')
    export_results(job, output_files['rnaChimeric.out.junction'], 'rna_chimeric.junction',
                   univ_options, subfolder='mutations/fusions')
    job.fileStore.logToMaster('Ran STAR on %s successfully' % univ_options['patient'])
    return output_files


def sort_and_index_star(job, star_bams, univ_options, star_options):
    """
    A wrapper for sorting and indexing the genomic star bam generated by run_star. It is required
    since run_star returns a dict of 2 bams

    :param dict star_bams: The bams from run_star
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict star_options: Options specific to star
    :return: Dict containing input bam and the generated index (.bam.bai)
                     output_files:
                        |- 'rna_transcriptome.bam': fsID
                        +- 'rna_genome':
                                 |- 'rna_sorted.bam': fsID
                                 +- 'rna_sorted.bam.bai': fsID
                        +- 'rnaChimeric.out.junction': fsID
    :rtype: dict
    """
    star_options['samtools']['n'] = star_options['n']
    sort = job.wrapJobFn(sort_bamfile, star_bams['rnaAligned.out.bam'], 'rna', univ_options,
                         samtools_options=star_options['samtools'],
                         disk=PromisedRequirement(sort_disk, star_bams['rnaAligned.out.bam']))
    index = job.wrapJobFn(index_bamfile, sort.rv(), 'rna', univ_options,
                          samtools_options=star_options['samtools'], sample_info='genome_sorted',
                          disk=PromisedRequirement(index_disk, sort.rv()))
    job.addChild(sort)
    sort.addChild(index)
    return {'rna_genome': index.rv(),
            'rna_transcriptome.bam': star_bams['rnaAligned.toTranscriptome.out.bam'],
            'rnaChimeric.out.junction': star_bams['rnaChimeric.out.junction']}
