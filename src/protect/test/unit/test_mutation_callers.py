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
from protect.common import untargz
from protect.mutation_calling.muse import run_muse
from protect.mutation_calling.mutect import run_mutect
from protect.mutation_calling.radia import run_radia
from protect.mutation_calling.fusion import run_fusion
from protect.mutation_calling.somaticsniper import run_somaticsniper
from protect.mutation_calling.strelka import run_strelka
from protect.pipeline.ProTECT import _parse_config_file
from protect.test import ProtectTest

from toil.job import Job

import os
import subprocess
import unittest


class TestMutationCallers(ProtectTest):
    def setUp(self):
        super(TestMutationCallers, self).setUp()
        test_dir = self._createTempDir()
        self.options = Job.Runner.getDefaultOptions(self._getTestJobStorePath())
        self.options.logLevel = 'INFO'
        self.options.workDir = test_dir
        self.options.clean = 'always'

    def test_mutect(self):
        self._test_dna_based_callers('mutect', run_mutect)

    def test_muse(self):
        self._test_dna_based_callers('muse', run_muse)

    def test_somaticsniper(self):
        self._test_dna_based_callers('somaticsniper', run_somaticsniper)

    def test_strelka(self):
        self._test_dna_based_callers('strelka', run_strelka)

    def _test_dna_based_callers(self, mut_caller, mut_caller_fn):
        """
        Test the functionality of MuSE, MuTect and Somatic Sniper
        """
        univ_options = self._getTestUnivOptions()
        config_file = os.path.join(self._projectRootPath(),
                                   'src/protect/test/test_inputs/ci_parameters.yaml')
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a1 = Job.wrapJobFn(self._get_test_dna_alignments, 'tumor_dna')
        a2 = Job.wrapJobFn(self._get_test_dna_alignments, 'normal_dna')
        b1 = Job.wrapJobFn(self._get_all_tools, config_file).encapsulate()
        b2 = Job.wrapJobFn(self._get_tool, b1.rv(), mut_caller)
        c = Job.wrapJobFn(mut_caller_fn, a1.rv(), a2.rv(), univ_options, b2.rv()).encapsulate()
        a1.addChild(a2)
        a2.addChild(b1)
        b1.addChild(b2)
        b2.addChild(c)
        Job.Runner.startToil(a1, self.options)


    def test_star_fusion(self):
        """
        Test the functionality of STAR-Fusion pipeline
        """
        univ_options = self._getTestUnivOptions()
        a1 = Job.wrapJobFn(self._get_test_fusion_reads)
        a2 = Job.wrapJobFn(self._get_test_fusion_junction)
        b = Job.wrapJobFn(self._get_fusion_options).encapsulate()
        c = Job.wrapJobFn(run_fusion,
                          a1.rv(),
                          a1.rv(),
                          univ_options,
                          b.rv('star_fusion'),
                          b.rv('fusion_inspector')).encapsulate()
        a1.addChild(a2)
        a1.addChild(b)
        b.addChild(c)
        Job.Runner.startToil(a1, self.options)

    @staticmethod
    def _get_fusion_options(job):
        star_fusion_options = {}
        fusion_inspector_options = {}
        call = 's3am download S3://cgl-pipeline-inputs/protect/ci_references/ci_star_fusion_compatible_index.tar.gz ./index.tar.gz'
        subprocess.check_call(call.split(' '))
        star_fusion_options['index'] = fusion_inspector_options['index'] = job.fileStore.writeGlobalFile('index.tar.gz')
        # Trinity now sets a minimum for the number of reads, so don't run it
        # on test data
        star_fusion_options['n'] = 2
        star_fusion_options['version'] = '1.0.0'
        fusion_inspector_options['run_trinity'] = False
        fusion_inspector_options['version'] = '1.0.1'
        return {'star_fusion': star_fusion_options, 'fusion_inspector': fusion_inspector_options}


    @unittest.skip('Takes too long')
    def test_radia(self):
        self._test_mixed_callers('radia', run_radia)

    def _test_mixed_callers(self, mut_caller, mut_caller_fn):
        """
        Test the functionality of radia
        """
        univ_options = self._getTestUnivOptions()
        config_file = os.path.join(self._projectRootPath(),
                                   'src/protect/test/test_inputs/ci_parameters.yaml')
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a1 = Job.wrapJobFn(self._get_test_rna_alignments)
        a2 = Job.wrapJobFn(self._get_test_dna_alignments, 'tumor_dna')
        a3 = Job.wrapJobFn(self._get_test_dna_alignments, 'normal_dna')
        b1 = Job.wrapJobFn(self._get_all_tools, config_file).encapsulate()
        b2 = Job.wrapJobFn(self._get_tool, b1.rv(), mut_caller)
        c = Job.wrapJobFn(mut_caller_fn, a1.rv(), a2.rv(), a3.rv(), univ_options, b2.rv()
                          ).encapsulate()
        a1.addChild(a2)
        a2.addChild(a3)
        a3.addChild(b1)
        b1.addChild(b2)
        b2.addChild(c)
        Job.Runner.startToil(a1, self.options)

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
    def _get_test_dna_alignments(job, sample_type):
        """
        Get the test bam and bai from s3

        :return: FSID for the bam and bai
        """
        assert sample_type in ('tumor_dna', 'normal_dna')
        bamfile = sample_type + '_fix_pg_sorted.bam'
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/unit_results/alignments/'
        final_call = base_call + bamfile + ' ' + bamfile
        subprocess.check_call(final_call.split(' '))
        final_call = base_call + bamfile + '.bai ' + bamfile + '.bai'
        subprocess.check_call(final_call.split(' '))
        return {bamfile: job.fileStore.writeGlobalFile(bamfile),
                bamfile + '.bai': job.fileStore.writeGlobalFile(bamfile + '.bai')}

    @staticmethod
    def _get_test_rna_alignments(job):
        """
        Get the test bam and bai from s3

        :return: FSID for the bam and bai
        """
        sample_type = 'rna'
        bamfile = sample_type + '_genome_sorted.bam'
        base_call = 's3am download s3://cgl-pipeline-inputs/protect/unit_results/alignments/'
        final_call = base_call + bamfile + ' ' + bamfile
        subprocess.check_call(final_call.split(' '))
        return {'rna_genome': {
            bamfile: job.fileStore.writeGlobalFile(bamfile),
            bamfile + '.bai': job.fileStore.writeGlobalFile(bamfile + '.bai')}}


    @staticmethod
    def _get_test_fusion_reads(job):
        """
        Get the test fusion reads from s3

        :return: FSID for each paired FASTQ
        """

        base_call = 's3am download S3://cgl-pipeline-inputs/protect/ci_references/'
        samples = ['RNA_CD74_ROS1_1.fq.gz', 'RNA_CD74_ROS1_2.fq.gz']
        for sample in samples:
            call = '{base}{sample} ./{sample}'.format(base=base_call, sample=sample)
            subprocess.check_call(call.split(' '))
        r1 = job.fileStore.writeGlobalFile('RNA_CD74_ROS1_1.fq.gz')
        r2 = job.fileStore.writeGlobalFile('RNA_CD74_ROS1_2.fq.gz')
        return r1, r2

    @staticmethod
    def _get_test_fusion_junction(job):
        """
        Get the test fusion reads from s3

        :return: FSID for each paired FASTQ
        """

        base_call = 's3am download S3://cgl-pipeline-inputs/protect/ci_references/'
        sample = 'CD74_ROS1_Chimeric.out.junction'
        call = '{base}{sample} ./ChimericJunction'.format(base=base_call, sample=sample)
        subprocess.check_call(call.split(' '))
        return job.fileStore.writeGlobalFile('ChimericJunction')


_get_all_tools = TestMutationCallers._get_all_tools
_get_tool = TestMutationCallers._get_tool
_get_test_dna_alignments = TestMutationCallers._get_test_dna_alignments
_get_test_rna_alignments = TestMutationCallers._get_test_rna_alignments
_get_fusion_options = TestMutationCallers._get_fusion_options
_get_test_fusion_reads = TestMutationCallers._get_test_fusion_reads
_get_test_fusion_junction = TestMutationCallers._get_test_fusion_junction
