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
File : protect/test/test_spawn_antigen_predictors.py
"""
from __future__ import print_function

from protect.binding_prediction.common import spawn_antigen_predictors, merge_mhc_peptide_calls
from protect.pipeline.ProTECT import _parse_config_file
from protect.test import ProtectTest
from toil.job import Job

import os
import subprocess


class TestSpawnAntigenPredictorsAndMerge(ProtectTest):
    def setUp(self):
        super(TestSpawnAntigenPredictorsAndMerge, self).setUp()
        test_dir = self._createTempDir()
        self.options = Job.Runner.getDefaultOptions(self._getTestJobStorePath())
        self.options.logLevel = 'INFO'
        self.options.workDir = test_dir
        self.options.clean = 'always'

    def test_spawn_antigen_predictors(self):
        """
        Test the functionality of spawn_antigen_predictors
        """
        univ_options = self._getTestUnivOptions()
        univ_options['output_folder'] = '/mnt/ephemeral/done'
        config_file = os.path.join(self._projectRootPath(),
                                   'src/protect/test/test_inputs/ci_parameters.yaml')
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a = Job.wrapJobFn(self._get_test_transgene_files)
        b = Job.wrapJobFn(self._get_test_phlat_files)
        c = Job.wrapJobFn(self._get_all_tools, config_file).encapsulate()
        d = Job.wrapJobFn(self._get_tool, c.rv(), 'mhci')
        e = Job.wrapJobFn(self._get_tool, c.rv(), 'mhcii')
        f = Job.wrapJobFn(spawn_antigen_predictors, a.rv(), b.rv(), univ_options, (d.rv(), e.rv()),
                          disk='100M', memory='100M', cores=1).encapsulate()
        g = Job.wrapJobFn(merge_mhc_peptide_calls, f.rv(), a.rv(), univ_options, disk='100M',
                          memory='100M', cores=1)
        a.addChild(b)
        a.addChild(g)
        b.addChild(c)
        c.addChild(d)
        c.addChild(e)
        d.addChild(f)
        e.addChild(f)
        f.addChild(g)
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
    def _get_test_phlat_files(job):
        """
        Get the test phlat file and write to jobstore

        :return: FSID for the phlat file
        """
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/unit_results/haplotyping/'
        phlat_files = {}
        for filename in ['mhci_alleles.list', 'mhcii_alleles.list']:
            call = (base_call + ('%s ' % filename) * 2).strip().split(' ')
            subprocess.check_call(call)
            phlat_files[filename] = job.fileStore.writeGlobalFile(filename)
        return phlat_files

# noinspection PyProtectedMember
_get_all_tools = TestSpawnAntigenPredictorsAndMerge._get_all_tools
# noinspection PyProtectedMember
_get_tool = TestSpawnAntigenPredictorsAndMerge._get_tool
# noinspection PyProtectedMember
_get_test_transgene_files = TestSpawnAntigenPredictorsAndMerge._get_test_transgene_files
# noinspection PyProtectedMember
_get_test_phlat_files = TestSpawnAntigenPredictorsAndMerge._get_test_phlat_files