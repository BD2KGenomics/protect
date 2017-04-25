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
File : protect/test/test_rankboost.py
"""
from __future__ import print_function

from protect.pipeline.ProTECT import _parse_config_file
from protect.rankboost import wrap_rankboost
from protect.test import ProtectTest
from toil.job import Job

import os
import subprocess


class TestRankboost(ProtectTest):
    def setUp(self):
        super(TestRankboost, self).setUp()
        test_dir = self._createTempDir()
        self.options = Job.Runner.getDefaultOptions(self._getTestJobStorePath())
        self.options.logLevel = 'INFO'
        self.options.workDir = test_dir
        self.options.clean = 'always'

    def test_rank_boost(self):
        """
        Test the functionality of spawn_antigen_predictors
        """
        univ_options = self._getTestUnivOptions()
        univ_options['output_folder'] = '/mnt/ephemeral/done'
        config_file = os.path.join(self._projectRootPath(),
                                   'src/protect/test/test_inputs/ci_parameters.yaml')
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a = Job.wrapJobFn(self._get_test_rsem_files)
        b = Job.wrapJobFn(self._get_test_merge_mhc_files)
        c = Job.wrapJobFn(self._get_test_transgene_files)
        d = Job.wrapJobFn(self._get_all_tools, config_file).encapsulate()
        e = Job.wrapJobFn(self._get_tool, d.rv(), 'rankboost')
        f = Job.wrapJobFn(wrap_rankboost, a.rv(), b.rv(), c.rv(), univ_options, e.rv(), disk='100M',
                          memory='100M', cores=1).encapsulate()
        a.addChild(b)
        b.addChild(c)
        c.addChild(d)
        d.addChild(e)
        e.addChild(f)
        Job.Runner.startToil(a, self.options)

    @staticmethod
    def _get_all_tools(job, config_file):
        sample_set, univ_options, tool_options = _parse_config_file(job, config_file,
                                                                    max_cores=None)
        return tool_options

    @staticmethod
    def _get_tool(job, all_tools, tool):
        return all_tools[tool]

    @staticmethod
    def _get_test_transgene_files(job):
        """
        Get the test transgene file and write to jobstore

        :return: FSID for the tansgene file
        """
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/unit_results/peptides/'
        transgened_files = {}
        filenames = []
        for length in ['9', '10', '15']:
            for tissue in ['tumor', 'normal']:
                filename = '_'.join(['transgened', tissue, length, 'mer_snpeffed.faa'])
                filenames.append(filename)
                if length != '9' and tissue == 'tumor':
                    filenames.append(filename + '.map')
        for filename in filenames:
            call = (base_call + ('%s ' % filename)*2).strip().split(' ')
            subprocess.check_call(call)
            transgened_files[filename] = job.fileStore.writeGlobalFile(filename)
        return transgened_files

    @staticmethod
    def _get_test_merge_mhc_files(job):
        """
        Get the test merge_mhc files and write to jobstore

        :return: FSID for the phlat file
        """
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/unit_results/' \
                    'binding_predictions/'
        merge_mhc_files = {}
        for filename in ['mhci_merged_files.list', 'mhcii_merged_files.list']:
            call = (base_call + ('%s ' % filename) * 2).strip().split(' ')
            subprocess.check_call(call)
            merge_mhc_files[filename] = job.fileStore.writeGlobalFile(filename)
        return merge_mhc_files

    @staticmethod
    def _get_test_rsem_files(job):
        """
        Get the test rsem file and write to jobstore

        :return: FSID for the phlat file
        """
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/unit_results/expression/'
        rsem_files = {}
        filename = 'rsem.isoforms.results'
        call = (base_call + ('%s ' % filename) * 2).strip().split(' ')
        subprocess.check_call(call)
        rsem_files[filename] = job.fileStore.writeGlobalFile(filename)
        return rsem_files


# noinspection PyProtectedMember
_get_all_tools = TestRankboost._get_all_tools
# noinspection PyProtectedMember
_get_tool = TestRankboost._get_tool
# noinspection PyProtectedMember
_get_test_transgene_files = TestRankboost._get_test_transgene_files
# noinspection PyProtectedMember
_get_test_merge_mhc_files = TestRankboost._get_test_merge_mhc_files
# noinspection PyProtectedMember
_get_test_rsem_files = TestRankboost._get_test_rsem_files