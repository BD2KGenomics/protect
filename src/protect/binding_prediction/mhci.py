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

import os
from protect.common import docker_call, get_files_from_filestore


def predict_mhci_binding(job, peptfile, allele, peplen, univ_options,
                         mhci_options):
    """
    This module will predict MHC:peptide binding for peptides in the files created in node XX to
    ALLELE.  ALLELE represents an MHCI allele.

    This module corresponds to node 18 on the tree
    """
    job.fileStore.logToMaster('Running mhci on %s:%s:%s' % (univ_options['patient'], allele,
                                                            peplen))
    work_dir = os.getcwd()
    input_files = {
        'peptfile.faa': peptfile}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = [mhci_options['pred'],
                  allele,
                  peplen,
                  input_files['peptfile.faa']]
    with open('/'.join([work_dir, 'predictions.tsv']), 'w') as predfile:
        docker_call(tool='mhci:2.13', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=predfile, interactive=True)
    output_file = job.fileStore.writeGlobalFile(predfile.name)
    return output_file
