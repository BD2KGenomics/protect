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


def run_fusion_caller(job, star_bam, univ_options, fusion_options):
    """
    This module will run a fusion caller on DNA bams.  This module will be
    implemented in the future.

    This module corresponds to node 10 on the tree
    """
    job.fileStore.logToMaster('Running FUSION on %s' % univ_options['patient'])
    job.fileStore.logToMaster('Fusions are currently unsupported.... Skipping.')
    fusion_file = job.fileStore.getLocalTempFile()
    output_file = job.fileStore.writeGlobalFile(fusion_file)
    return output_file
