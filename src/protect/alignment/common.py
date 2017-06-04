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
from math import ceil
from protect.common import docker_call, docker_path, export_results, get_files_from_filestore

import os


# disk for indexing
def index_disk(bamfile):
    return int(1.1 * ceil(bamfile.size + 524288))


# disk for sorting
def sort_disk(bamfile):
    return int(2.5 * ceil(bamfile.size + 524288))


def index_bamfile(job, bamfile, sample_type, univ_options, samtools_options, sample_info=None,
                  export=True):
    """
    Index `bamfile` using samtools

    :param toil.fileStore.FileID bamfile: fsID for the bam file
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict samtools_options: Options specific to samtools
    :param str sample_info: Information regarding the sample that will beinjected into the filename
               as `sample_type`_`sample_info`.bam(.bai)
    :param bool export: Should the bam and bai be exported to the output directory?
    :return: Dict containing input bam and the generated index (.bam.bai)
             output_files:
                 |- '<sample_type>(_<sample_info>).bam': fsID
                 +- '<sample_type>(_<sample_info>).bam.bai': fsID
    :rtype: dict
    """
    job.fileStore.logToMaster('Running samtools-index on %s:%s' % (univ_options['patient'],
                                                                   sample_type))
    work_dir = os.getcwd()
    in_bamfile = sample_type
    if sample_info is not None:
        assert isinstance(sample_info, str)
        in_bamfile = '_'.join([in_bamfile, sample_info])
    in_bamfile += '.bam'
    input_files = {
        in_bamfile: bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['index',
                  input_files[in_bamfile]]
    docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=samtools_options['version'])
    out_bai = '/'.join([work_dir, in_bamfile + '.bai'])
    output_files = {in_bamfile: bamfile,
                    in_bamfile + '.bai': job.fileStore.writeGlobalFile(out_bai)}
    if export:
        export_results(job, bamfile, os.path.splitext(out_bai)[0], univ_options,
                       subfolder='alignments')
        export_results(job, output_files[in_bamfile + '.bai'], out_bai, univ_options,
                       subfolder='alignments')
    return output_files


def sort_bamfile(job, bamfile, sample_type, univ_options, samtools_options):
    """
    Sort `bamfile` using samtools

    :param toil.fileStore.FileID bamfile: fsID for the bam file
    :param str sample_type: Description of the sample to inject into the filename
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict samtools_options: Options specific to samtools
    :return: fsID for the sorted bamfile
    :rtype: toil.fileStore.FileID
    """
    job.fileStore.logToMaster('Running samtools-sort on %s:%s' % (univ_options['patient'],
                                                                  sample_type))
    work_dir = os.getcwd()
    in_bamfile = ''.join([sample_type, '.bam'])
    out_bamfile = '_'.join([sample_type, 'sorted.bam'])
    input_files = {
        in_bamfile: bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['sort',
                  '-o', docker_path(out_bamfile),
                  '-O', 'bam',
                  '-T', 'temp_sorted',
                  '-@', str(samtools_options['n']),
                  input_files[in_bamfile]]
    docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=samtools_options['version'])
    job.fileStore.deleteGlobalFile(bamfile)
    return job.fileStore.writeGlobalFile(out_bamfile)
