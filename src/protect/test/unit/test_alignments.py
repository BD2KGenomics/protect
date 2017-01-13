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
File : protect/test/test_file_downloads.py
"""
from __future__ import print_function
import os
import subprocess

from protect.alignment.dna import align_dna
from protect.alignment.rna import align_rna
from protect.pipeline.ProTECT import _parse_config_file
from protect.test import ProtectTest

from toil.job import Job


class TestAlignments(ProtectTest):
    def setUp(self):
        super(TestAlignments, self).setUp()
        test_dir = self._createTempDir()
        self.options = Job.Runner.getDefaultOptions(self._getTestJobStorePath())
        self.options.logLevel = 'INFO'
        self.options.workDir = test_dir
        self.options.clean = 'always'

    def test_bwa(self):
        """
        Test the functionality of align_dna
        """
        univ_options = self._getTestUnivOptions()
        config_file = os.path.join(self._projectRootPath(),
                                   'src/protect/test/test_inputs/ci_parameters.yaml')
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a = Job.wrapJobFn(self._get_test_bwa_files)
        b = Job.wrapJobFn(self._get_all_tools, config_file).encapsulate()
        c = Job.wrapJobFn(self._get_tool, b.rv(), 'bwa')
        d = Job.wrapJobFn(align_dna, a.rv(), 'tumor_dna', univ_options, c.rv()).encapsulate()
        a.addChild(b)
        b.addChild(c)
        c.addChild(d)
        Job.Runner.startToil(a, self.options)

    @staticmethod
    def _get_all_tools(job, config_file):
        sample_set, univ_options, tool_options = _parse_config_file(job, config_file,
                                                                    max_cores=None)
        return tool_options

    @staticmethod
    def _get_tool(job, all_tools, tool):
        all_tools[tool]['n'] = 2
        return all_tools[tool]

    @staticmethod
    def _get_test_bwa_files(job):
        """
        Get the test rsem file and write to jobstore

        :return: FSID for the rsem file
        """
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/ci_references/'
        subprocess.check_call((base_call + 'Tum_1.fq.gz Tum_1.fq.gz').split(' '))
        subprocess.check_call((base_call + 'Tum_2.fq.gz Tum_2.fq.gz').split(' '))
        return [job.fileStore.writeGlobalFile('Tum_1.fq.gz'),
                job.fileStore.writeGlobalFile('Tum_2.fq.gz')]

    def test_star(self):
        """
        Test the functionality of align_dna
        """
        univ_options = self._getTestUnivOptions()
        config_file = os.path.join(self._projectRootPath(),
                                   'src/protect/test/test_inputs/ci_parameters.yaml')
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a = Job.wrapJobFn(self._get_test_star_files)
        b = Job.wrapJobFn(self._get_all_tools, config_file).encapsulate()
        c = Job.wrapJobFn(self._get_tool, b.rv(), 'star')
        d = Job.wrapJobFn(align_rna, a.rv(), univ_options, c.rv()).encapsulate()
        a.addChild(b)
        b.addChild(c)
        c.addChild(d)
        Job.Runner.startToil(a, self.options)

    @staticmethod
    def _get_test_star_files(job):
        """
        Get the test rsem file and write to jobstore

        :return: FSID for the rsem file
        """
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/ci_references/'
        subprocess.check_call((base_call + 'Rna_1.fq.gz Rna_1.fq.gz').split(' '))
        subprocess.check_call((base_call + 'Rna_2.fq.gz Rna_2.fq.gz').split(' '))
        return [job.fileStore.writeGlobalFile('Rna_1.fq.gz'),
                job.fileStore.writeGlobalFile('Rna_2.fq.gz')]


# noinspection PyProtectedMember
_get_all_tools = TestAlignments._get_all_tools
# noinspection PyProtectedMember
_get_tool = TestAlignments._get_tool
# noinspection PyProtectedMember
_get_test_bwa_files = TestAlignments._get_test_bwa_files
# noinspection PyProtectedMember
_get_test_star_files = TestAlignments._get_test_star_files
