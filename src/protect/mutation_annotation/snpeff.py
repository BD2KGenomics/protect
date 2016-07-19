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
from protect.common import (docker_call, get_files_from_filestore, export_results, untargz,
                            docker_path)

import os


# disk for snpeff.
def snpeff_disk(snpeff_index):
    return int(6 * ceil(snpeff_index.size + 524288))


def run_snpeff(job, merged_mutation_file, univ_options, snpeff_options):
    """
    This module will run snpeff on the aggregated mutation calls.  Currently the only mutations
    called are SNPs hence SnpEff suffices. This node will be replaced in the future with another
    translator.

    ARGUMENTS
    1. merged_mutation_file: <JSid for merged vcf>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. snpeff_options: Dict of parameters specific to snpeff
         snpeff_options
                +- 'tool_index': <JSid for the snpEff index tarball>

    RETURN VALUES
    1. output_file: <JSid for the snpeffed vcf>

    This node corresponds to node 16 on the tree
    """
    job.fileStore.logToMaster('Running snpeff on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'merged_mutations.vcf': merged_mutation_file,
        'snpeff_index.tar.gz': snpeff_options['tool_index']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    input_files['snpeff_index'] = untargz(input_files['snpeff_index.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['eff',
                  '-dataDir', input_files['snpeff_index'],
                  '-c', '/'.join([input_files['snpeff_index'], 'snpEff_hg19_gencode.config']),
                  '-no-intergenic',
                  '-no-downstream',
                  '-no-upstream',
                  # '-canon',
                  '-noStats',
                  'hg19_gencode',
                  input_files['merged_mutations.vcf']]
    xmx = snpeff_options['java_Xmx'] if snpeff_options['java_Xmx'] else univ_options['java_Xmx']
    with open('/'.join([work_dir, 'mutations.vcf']), 'w') as snpeff_file:
        docker_call(tool='snpeff', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], java_opts=xmx, outfile=snpeff_file)
    export_results(job, snpeff_file.name, univ_options, subfolder='mutations/snpeffed')
    output_file = job.fileStore.writeGlobalFile(snpeff_file.name)
    return output_file
