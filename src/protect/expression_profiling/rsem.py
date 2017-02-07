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
from __future__ import print_function
from math import ceil
from protect.common import (docker_call, get_files_from_filestore, untargz, docker_path,
                            export_results)
from toil.job import PromisedRequirement

import os


# disk for rsem
def rsem_disk(star_bams, rsem_index):
    star_transcriptome_bam = star_bams['rnaAligned.sortedByCoord.out.bam']['rna_fix_pg_sorted.bam']
    return int(3 * ceil(star_transcriptome_bam.size + 524288) +
               4 * ceil(rsem_index.size + 524288))


def wrap_rsem(job, star_bams, univ_options, rsem_options):
    """
    This is a convenience function that runs rsem from the star outputs.

    :param job job: job
    :param dict star_bams: dict of results from star
    :param dict univ_options: Universal Options
    :param dict rsem_options: Options specific to rsem
    :return:
    """
    rsem = job.addChildJobFn(run_rsem, star_bams['rnaAligned.toTranscriptome.out.bam'],
                             univ_options, rsem_options, cores=rsem_options['n'],
                             disk=PromisedRequirement(rsem_disk, star_bams,
                                                      rsem_options['tool_index']))

    return rsem.rv()


def run_rsem(job, rna_bam, univ_options, rsem_options):
    """
    This module will run rsem on the RNA Bam file.

    ARGUMENTS
    1. rna_bam: <JSid of rnaAligned.toTranscriptome.out.bam>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. rsem_options: Dict of parameters specific to rsem
         rsem_options
              |- 'tool_index': <JSid for the rsem index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_file: <Jsid of rsem.isoforms.results>

    This module corresponds to node 9 on the tree
    """
    job.fileStore.logToMaster('Running rsem on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'star_transcriptome.bam': rna_bam,
        'rsem_index.tar.gz': rsem_options['tool_index']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    input_files['rsem_index'] = untargz(input_files['rsem_index.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['--paired-end',
                  '-p', str(rsem_options['n']),
                  '--bam',
                  input_files['star_transcriptome.bam'],
                  '--no-bam-output',
                  '/'.join([input_files['rsem_index'], 'hg19']),
                  'rsem']
    docker_call(tool='rsem:1.2.20', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_files = {}
    for filename in ('rsem.genes.results', 'rsem.isoforms.results'):
        output_files[filename] = job.fileStore.writeGlobalFile('/'.join([work_dir, filename]))
        export_results(job, '/'.join([work_dir, filename]), univ_options, subfolder='expression')
    return output_files
