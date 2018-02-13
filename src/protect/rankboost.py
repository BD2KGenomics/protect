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
                   rankboost_options):
    """
    A wrapper for boost_ranks.

    :param dict rsem_files: Dict of results from rsem
    :param dict merged_mhc_calls: Dict of results from merging mhc peptide binding predictions
    :param dict transgene_out: Dict of results from running Transgene
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict rankboost_options: Options specific to rankboost
    :return: Dict of concise and detailed results for mhci and mhcii
             output_files:
                |- 'mhcii_rankboost_concise_results.tsv': fsID
                |- 'mhcii_rankboost_detailed_results.txt': fsID
                |- 'mhci_rankboost_concise_results.tsv': fsID
                +- 'mhci_rankboost_detailed_results.txt': fsID
    :rtype: dict
    """
    rankboost = job.addChildJobFn(boost_ranks, rsem_files['rsem.isoforms.results'],
                                  merged_mhc_calls, transgene_out, univ_options, rankboost_options)

    return rankboost.rv()


def boost_ranks(job, isoform_expression, merged_mhc_calls, transgene_out, univ_options,
                rankboost_options):
    """
    Boost the ranks of the predicted peptides:MHC combinations.

    :param toil.fileStore.FileID isoform_expression: fsID of rsem isoform expression file
    :param dict merged_mhc_calls: Dict of results from merging mhc peptide binding predictions
    :param dict transgene_out: Dict of results from running Transgene
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict rankboost_options: Options specific to rankboost
    :return: Dict of concise and detailed results for mhci and mhcii
             output_files:
                |- 'mhcii_rankboost_concise_results.tsv': fsID
                |- 'mhcii_rankboost_detailed_results.txt': fsID
                |- 'mhci_rankboost_concise_results.tsv': fsID
                +- 'mhci_rankboost_detailed_results.txt': fsID
    :rtype: dict
    """
    work_dir = os.getcwd()
    input_files = {
        'rsem_quant.tsv': isoform_expression,
        'mhci_merged_files.tsv': merged_mhc_calls['mhci_merged_files.list'],
        'mhcii_merged_files.tsv': merged_mhc_calls['mhcii_merged_files.list'],
        'mhci_peptides.faa': transgene_out['transgened_tumor_10_mer_snpeffed.faa'],
        'mhcii_peptides.faa': transgene_out['transgened_tumor_15_mer_snpeffed.faa']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    output_files = {}
    for mhc in ('mhci', 'mhcii'):
        import re
        ratios = re.sub("'", '', repr(rankboost_options[''.join([mhc, '_args'])]))
        parameters = ['--' + mhc,
                      '--predictions', input_files[''.join([mhc, '_merged_files.tsv'])],
                      '--expression', input_files['rsem_quant.tsv'],
                      '--peptides', input_files[''.join([mhc, '_peptides.faa'])],
                      '--ratios', ratios
                      ]
        docker_call(tool='rankboost', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], tool_version=rankboost_options['version'])
        mhc_concise = ''.join([work_dir, '/', mhc, '_rankboost_concise_results.tsv'])
        mhc_detailed = ''.join([work_dir, '/', mhc, '_rankboost_detailed_results.txt'])
        output_files[mhc] = {}
        if os.path.exists(mhc_concise):
            output_files[os.path.basename(mhc_concise)] = job.fileStore.writeGlobalFile(mhc_concise)
            export_results(job, output_files[os.path.basename(mhc_concise)], mhc_concise,
                           univ_options, subfolder='rankboost')
        else:
            output_files[os.path.basename(mhc_concise)] = None
        if os.path.exists(mhc_detailed):
            output_files[os.path.basename(mhc_detailed)] = \
                job.fileStore.writeGlobalFile(mhc_detailed)
            export_results(job, output_files[os.path.basename(mhc_detailed)], mhc_detailed,
                           univ_options, subfolder='rankboost')
        else:
            output_files[os.path.basename(mhc_detailed)] = None
    job.fileStore.logToMaster('Ran boost_ranks on %s successfully' % univ_options['patient'])
    return output_files
