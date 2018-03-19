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

from protect.common import docker_call, docker_path, get_files_from_filestore, is_gzipfile

import os


# disk for cutadapt
def cutadapt_disk(rna_fastqs):
    return int(2.5 * ceil(sum([f.size for f in rna_fastqs]) + 524288))


def run_cutadapt(job, fastqs, univ_options, cutadapt_options):
    """
    Runs cutadapt on the input RNA fastq files.

    :param list fastqs: List of fsIDs for input an RNA-Seq fastq pair
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict cutadapt_options: Options specific to cutadapt
    :return: List of fsIDs of cutadapted fastqs
    :rtype: list[toil.fileStore.FileID]
    """
    work_dir = os.getcwd()
    input_files = {
        'rna_1.fastq': fastqs[0],
        'rna_2.fastq': fastqs[1]}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    # Handle gzipped file
    gz = '.gz' if is_gzipfile(input_files['rna_1.fastq']) else ''
    if gz:
        for read_file in 'rna_1.fastq', 'rna_2.fastq':
            os.symlink(read_file, read_file + gz)
            input_files[read_file + gz] = input_files[read_file] + gz
    input_files = {key: docker_path(path) for key, path in input_files.items()}
    parameters = ['-a', cutadapt_options['a'],  # Fwd read 3' adapter
                  '-A', cutadapt_options['A'],  # Rev read 3' adapter
                  '-m', '35',  # Minimum size of read
                  '-o', docker_path('rna_cutadapt_1.fastq.gz'),  # Output for R1
                  '-p', docker_path('rna_cutadapt_2.fastq.gz'),  # Output for R2
                  input_files['rna_1.fastq' + gz],
                  input_files['rna_2.fastq' + gz]]
    docker_call(tool='cutadapt', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=cutadapt_options['version'])
    output_files = []
    for fastq_file in ['rna_cutadapt_1.fastq.gz', 'rna_cutadapt_2.fastq.gz']:
        output_files.append(job.fileStore.writeGlobalFile('/'.join([work_dir, fastq_file])))
    job.fileStore.logToMaster('Ran cutadapt on %s successfully' % univ_options['patient'])
    return output_files
