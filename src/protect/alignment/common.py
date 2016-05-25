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
from __future__ import absolute_import

import os
from protect.common import docker_call, get_files_from_filestore, export_results


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
    work_dir = os.getcwd()
    in_bamfile = '_'.join([sample_type, 'fix_pg_sorted.bam'])
    input_files = {
        in_bamfile: bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['index',
                  input_files[in_bamfile]]
    docker_call(tool='samtools', tool_parameters=parameters,
                work_dir=work_dir, dockerhub=univ_options['dockerhub'])
    out_bai = '/'.join([work_dir, in_bamfile + '.bai'])
    output_files = {in_bamfile: bamfile,
                    in_bamfile + '.bai': job.fileStore.writeGlobalFile(out_bai)}
    export_results(job, os.path.splitext(out_bai)[0], univ_options, subfolder='alignments')
    export_results(job, out_bai, univ_options, subfolder='alignments')
    return output_files
