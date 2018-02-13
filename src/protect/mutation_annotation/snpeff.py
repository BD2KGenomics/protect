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

from protect.common import (docker_call,
                            docker_path,
                            export_results,
                            get_files_from_filestore,
                            untargz)
import os


# disk for snpeff.
def snpeff_disk(snpeff_index):
    return int(6 * ceil(snpeff_index.size + 524288))


def run_snpeff(job, merged_mutation_file, univ_options, snpeff_options):
    """
    Run snpeff on an input vcf.

    :param toil.fileStore.FileID merged_mutation_file: fsID for input vcf
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict snpeff_options: Options specific to snpeff
    :return: fsID for the snpeffed vcf
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        'merged_mutations.vcf': merged_mutation_file,
        'snpeff_index.tar.gz': snpeff_options['index']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    input_files['snpeff_index'] = untargz(input_files['snpeff_index.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['eff',
                  '-dataDir', input_files['snpeff_index'],
                  '-c', '/'.join([input_files['snpeff_index'],
                                  'snpEff_' + univ_options['ref'] + '_gencode.config']),
                  '-no-intergenic',
                  '-no-downstream',
                  '-no-upstream',
                  # '-canon',
                  '-noStats',
                  univ_options['ref'] + '_gencode',
                  input_files['merged_mutations.vcf']]
    xmx = snpeff_options['java_Xmx'] if snpeff_options['java_Xmx'] else univ_options['java_Xmx']
    with open('/'.join([work_dir, 'mutations.vcf']), 'w') as snpeff_file:
        docker_call(tool='snpeff', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], java_xmx=xmx, outfile=snpeff_file,
                    tool_version=snpeff_options['version'])
    output_file = job.fileStore.writeGlobalFile(snpeff_file.name)
    export_results(job, output_file, snpeff_file.name, univ_options, subfolder='mutations/snpeffed')
    job.fileStore.logToMaster('Ran snpeff on %s successfully' % univ_options['patient'])
    return output_file
