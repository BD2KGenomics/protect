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

"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : protect/test/test_mhc_pathway_assessment.py
"""
from __future__ import print_function
from protect.test import ProtectTest

import os
import subprocess

class TestProtect(ProtectTest):
    def setUp(self):
        super(TestProtect, self).setUp()
        self.test_dir = self._createTempDir()
        self.work_dir = os.path.join(self.test_dir, 'working')
        os.makedirs(self.work_dir)

    def test_protect(self):
        protect_call = ['ProTECT',
                        '--config_file', os.path.join(self._projectRootPath(), 'src', 'protect',
                                                      'test', 'test_inputs',
                                                      'ci_parameters.yaml'),
                        '--workDir', self.work_dir,
                        os.path.join(self.test_dir, 'jobstore')]
        # Run ProTECT
        subprocess.check_call(protect_call)