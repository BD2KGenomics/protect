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
from protect.common import get_file_from_s3, untargz
from protect.addons import assess_mhc_genes
from toil.job import Job

import os


class TestMHCPathwayAssessment(ProtectTest):
    def setUp(self):
        super(TestMHCPathwayAssessment, self).setUp()
        test_dir = self._createTempDir()
        self.options = Job.Runner.getDefaultOptions(self._getTestJobStorePath())
        self.options.logLevel = 'INFO'
        self.options.workDir = test_dir
        self.options.clean = 'always'

    def test_mhc_assessment(self):
        """
        Test the functionality of assess_mhc_genes
        """
        univ_options = self._getTestUnivOptions()
        test_src_folder = os.path.join(self._projectRootPath(), 'src', 'protect', 'test')
        a = Job.wrapJobFn(self._get_test_rsem_file, test_src_folder)
        b = Job.wrapJobFn(self._get_MHC_file)
        c = Job.wrapJobFn(self._get_test_haplotype_file, test_src_folder)
        d = Job.wrapJobFn(assess_mhc_genes, a.rv(), c.rv(), univ_options, b.rv())
        e = Job.wrapJobFn(self._test_output, d.rv(), univ_options)
        a.addChild(b)
        b.addChild(c)
        c.addChild(d)
        d.addChild(e)
        Job.Runner.startToil(a, self.options)

    @staticmethod
    def _get_test_rsem_file(job, test_src_folder):
        """
        Get the test rsem file and write to jobstore

        :return: FSID for the rsem file
        """
        rsem_file = os.path.join(test_src_folder, 'test_inputs/test_rsem_quant.tsv.tar.gz')
        rsem_file = untargz(rsem_file, os.getcwd())
        return job.fileStore.writeGlobalFile(rsem_file)

    @staticmethod
    def _get_MHC_file(job):
        """
        Get the MHC file and write to jobstore

        :return: FSID for the MHC file
        """
        mhc_file = get_file_from_s3(job,
                                    'S3://cgl-protect-data/hg19_references/mhc_pathway_genes.json.tar.gz',
                                    write_to_jobstore=False)
        return {
            'genes_file': job.fileStore.writeGlobalFile(mhc_file)}

    @staticmethod
    def _get_test_haplotype_file(job, test_src_folder):
        """
        Get the test haplotype file and write to jobstore

        :return: FSID for the MHC file
        """
        rna_haplotype = os.path.join(test_src_folder, 'test_inputs/test_mhc_haplotype.sum.tar.gz')
        rna_haplotype = untargz(rna_haplotype, os.getcwd())
        return job.fileStore.writeGlobalFile(rna_haplotype)

    @staticmethod
    def _test_output(job, output_file, univ_options):
        """
        Test the results of the assessment

        :param output_file: The file created by assess_mhc_genes
        """
        outfile = job.fileStore.readGlobalFile(output_file, 'mhc_report.txt')
        # Ensure that the exported file exists
        assert os.path.exists(os.path.join(univ_options['output_folder'], 'reports',
                                           'mhc_pathway_report.txt'))
        # Ensure that the 2 input genes were processed correctly
        outdict = {}
        with open(outfile) as o_f:
            for line in o_f:
                line = line.strip().split()
                if len(line) == 0 or line[0] in ['Gene', 'TAP', 'MHC', 'MHCI', 'MHCII']:
                    continue
                else:
                    outdict[line[0]] = line[3]

        lows = {x for x in outdict if outdict[x] == 'LOW'}
        fails = {x for x in outdict if outdict[x] == 'FAIL'}
        passes = {x for x in outdict if outdict[x] == 'PASS'}
        expectedlows = {'TNF', 'HLA_A', 'HLA_B', 'HLA_DRA'}
        expectedfails = {'HLA_DQA'}
        expectedpasses = {'CTSL', 'HLA_C', 'HLA_DRB', 'HLA_DQB'}

        assert expectedlows.difference(lows) == set([])
        assert expectedfails.difference(fails) == set([])
        assert expectedpasses.difference(passes) == set([])


# noinspection PyProtectedMember
_get_test_rsem_file = TestMHCPathwayAssessment._get_test_rsem_file
# noinspection PyProtectedMember
_get_MHC_file = TestMHCPathwayAssessment._get_MHC_file
# noinspection PyProtectedMember
_test_output = TestMHCPathwayAssessment._test_output
# noinspection PyProtectedMember
_get_test_haplotype_file = TestMHCPathwayAssessment._get_test_haplotype_file