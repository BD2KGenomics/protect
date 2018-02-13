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

from protect.common import docker_call, get_files_from_filestore, read_peptide_file

import os


def predict_mhci_binding(job, peptfile, allele, peplen, univ_options, mhci_options):
    """
    Predict binding for each peptide in `peptfile` to `allele` using the IEDB mhci binding
    prediction tool.

    :param toil.fileStore.FileID peptfile: The input peptide fasta
    :param str allele: Allele to predict binding against
    :param str peplen: Length of peptides to process
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict mhci_options: Options specific to mhci binding prediction
    :return: fsID for file containing the predictions
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        'peptfile.faa': peptfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    peptides = read_peptide_file(os.path.join(os.getcwd(), 'peptfile.faa'))
    if not peptides:
        return job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())
    parameters = [mhci_options['pred'],
                  allele,
                  peplen,
                  input_files['peptfile.faa']]
    with open('/'.join([work_dir, 'predictions.tsv']), 'w') as predfile:
        docker_call(tool='mhci', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=predfile, interactive=True,
                    tool_version=mhci_options['version'])
    output_file = job.fileStore.writeGlobalFile(predfile.name)
    job.fileStore.logToMaster('Ran mhci on %s:%s:%s successfully'
                              % (univ_options['patient'], allele, peplen))
    return output_file
