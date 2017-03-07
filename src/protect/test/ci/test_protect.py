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
        self._test_ran_successfully()

    def _test_ran_successfully(self):
        """
        Ensure that the test created the correct outputs.
        :return: True or False
        """
        results_dir = '/mnt/ephemeral/done/TEST'
        results_dir_contents = [(a, sorted(b), sorted(c)) for a, b, c in os.walk(results_dir)]
        results_dir_contents.sort()
        # Obtained by running [[(a, sorted(b), sorted(c)) for a, b, c in os.walk(results_dir)]
        # + .sort() on a successful run.
        expected_contents = [('/mnt/ephemeral/done/TEST',
                              ['alignments', 'binding_predictions', 'expression', 'haplotyping',
                               'mutations', 'peptides', 'rankboost', 'reports'],
                              []),
                             ('/mnt/ephemeral/done/TEST/alignments',
                              [],
                              ['normal_dna_fix_pg_sorted.bam', 'normal_dna_fix_pg_sorted.bam.bai',
                               'rna_fix_pg_sorted.bam', 'rna_fix_pg_sorted.bam.bai',
                               'rna_transcriptome.bam', 'tumor_dna_fix_pg_sorted.bam',
                               'tumor_dna_fix_pg_sorted.bam.bai']),
                             ('/mnt/ephemeral/done/TEST/binding_predictions',
                              [],
                              ['mhci_merged_files.list', 'mhcii_merged_files.list']),
                             ('/mnt/ephemeral/done/TEST/expression',
                              [],
                              ['rsem.genes.results', 'rsem.isoforms.results']),
                             ('/mnt/ephemeral/done/TEST/haplotyping',
                              [],
                              ['mhci_alleles.list', 'mhcii_alleles.list']),
                             ('/mnt/ephemeral/done/TEST/mutations',
                              ['merged', 'muse', 'mutect', 'radia', 'snpeffed', 'somaticsniper',
                               'strelka', 'transgened'],
                              []),
                             ('/mnt/ephemeral/done/TEST/mutations/merged',
                              [],
                              ['all_merged.vcf', 'chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/muse',
                              [],
                              ['chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/mutect',
                              [],
                              ['chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/radia',
                              [],
                              ['chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/snpeffed',
                              [],
                              ['mutations.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/somaticsniper',
                              [],
                              ['chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/strelka',
                              ['indel', 'snv'],
                              []),
                             ('/mnt/ephemeral/done/TEST/mutations/strelka/indel',
                              [],
                              ['chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/strelka/snv',
                              [],
                              ['chr6.vcf']),
                             ('/mnt/ephemeral/done/TEST/mutations/transgened',
                              [],
                              ['mutations.vcf']),
                             ('/mnt/ephemeral/done/TEST/peptides',
                              [],
                              ['transgened_normal_10_mer_snpeffed.faa',
                               'transgened_normal_15_mer_snpeffed.faa',
                               'transgened_normal_9_mer_snpeffed.faa',
                               'transgened_tumor_10_mer_snpeffed.faa',
                               'transgened_tumor_10_mer_snpeffed.faa.map',
                               'transgened_tumor_15_mer_snpeffed.faa',
                               'transgened_tumor_15_mer_snpeffed.faa.map',
                               'transgened_tumor_9_mer_snpeffed.faa',
                               'transgened_tumor_9_mer_snpeffed.faa.map']),
                             ('/mnt/ephemeral/done/TEST/rankboost',
                              [],
                              ['mhci_rankboost_concise_results.tsv',
                               'mhci_rankboost_detailed_results.txt',
                               'mhcii_rankboost_concise_results.tsv',
                               'mhcii_rankboost_detailed_results.txt']),
                             ('/mnt/ephemeral/done/TEST/reports',
                              [],
                              ['mhc_pathway_report.txt'])]
        self.assertEqual(results_dir_contents, expected_contents)