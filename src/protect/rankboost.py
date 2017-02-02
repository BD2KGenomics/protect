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
from protect.common import docker_call, get_files_from_filestore, export_results

import os


def wrap_rankboost(job, rsem_files, merged_mhc_calls, transgene_out, univ_options,
                   rank_boost_options):
    """
    This is a convenience function that runs rankboost on pipeline outputs.

    :param job job: job
    :param dict rsem_files: dict of results from rsem
    :param dict merged_mhc_calls: dict of results from merging mhc peptide binding predictions
    :param dict transgene_out: dict of results from running transgene
    :param dict univ_options: Universal Options
    :param dict rank_boost_options: Options specific to rank boost
    :return:
    """
    rankboost = job.addChildJobFn(boost_ranks, rsem_files['rsem.isoforms.results'],
                                  merged_mhc_calls, transgene_out, univ_options, rank_boost_options)

    return rankboost.rv()


def boost_ranks(job, isoform_expression, merged_mhc_calls, transgene_out, univ_options,
                rank_boost_options):
    """
    This is the final module in the pipeline.  It will call the rank boosting R
    script.

    This module corresponds to node 21 in the tree
    """
    job.fileStore.logToMaster('Running boost_ranks on %s' % univ_options['patient'])
    work_dir = os.path.abspath(univ_options['patient'])
    os.mkdir(work_dir)
    input_files = {
        'rsem_quant.tsv': isoform_expression,
        'mhci_merged_files.tsv': merged_mhc_calls['mhci_merged_files.list'],
        'mhcii_merged_files.tsv': merged_mhc_calls['mhcii_merged_files.list'],
        'mhci_peptides.faa': transgene_out['transgened_tumor_10_mer_snpeffed.faa'],
        'mhcii_peptides.faa': transgene_out['transgened_tumor_15_mer_snpeffed.faa']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    output_files = {}
    for mhc in ('mhci', 'mhcii'):
        parameters = [mhc,
                      input_files[''.join([mhc, '_merged_files.tsv'])],
                      input_files['rsem_quant.tsv'],
                      input_files[''.join([mhc, '_peptides.faa'])],
                      rank_boost_options[''.join([mhc, '_combo'])]
                      ]
        docker_call(tool='rankboost:1.0.0', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'])
        mhc_concise = ''.join([work_dir, '/', mhc, '_merged_files_concise_results.tsv'])
        mhc_detailed = ''.join([work_dir, '/', mhc, '_merged_files_detailed_results.tsv'])
        output_files[mhc] = {}
        if os.path.exists(mhc_concise):
            output_files[os.path.basename(mhc_concise)] = job.fileStore.writeGlobalFile(mhc_concise)
            export_results(job, mhc_concise, univ_options, subfolder='rankboost')
        else:
            output_files[os.path.basename(mhc_concise)] = None
        if os.path.exists(mhc_detailed):
            output_files[os.path.basename(mhc_detailed)] = \
                job.fileStore.writeGlobalFile(mhc_detailed)
            export_results(job, mhc_detailed, univ_options, subfolder='rankboost')
        else:
            output_files[os.path.basename(mhc_detailed)] = None
    return output_files
