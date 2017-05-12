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
File : protect/test/test_reporting.py
"""
from __future__ import print_function

from protect.test import ProtectTest

import os
import pytest
import shutil
import subprocess


class TestProtect(ProtectTest):
    def setUp(self):
        super(TestProtect, self).setUp()
        self.test_dir = self._createTempDir()
        self.work_dir = os.path.join(self.test_dir, 'working')
        os.makedirs(self.work_dir)

    @pytest.mark.all_fastq
    def test_protect_all_fastq(self):
        # All inputs are fastq
        self._test_protect('ci_all_fastq_parameters.yaml')
        self._test_ran_successfully([{'alignments': ['normal_dna', 'rna', 'tumor_dna']},
                                     'binding_predictions', 'expression', 'haplotyping',
                                     {'mutations': ['fusions', 'merged', 'mutect', 'muse', 'radia',
                                                    'somaticsniper', 'strelka', 'snpeffed',
                                                    'transgened']}, 'peptides', 'rankboost',
                                     'reports'])

    @pytest.mark.mix_bam_fastq
    def test_protect_mix_bam_fastq(self):
        # Inputs are a mixture of fastqs, bams, and necessarily, the input haplotype
        self._test_protect('ci_mix_bam_fastq_parameters.yaml')
        self._test_ran_successfully([{'alignments': ['tumor_dna']}, 'binding_predictions',
                                     'expression', {'mutations': ['merged', 'mutect', 'muse',
                                                                  'radia', 'somaticsniper',
                                                                  'strelka', 'snpeffed',
                                                                  'transgened']}, 'peptides',
                                     'rankboost', 'reports'])

    @pytest.mark.vcf_fastq
    def test_vcf_fastq(self):
        # Inputs are fastqs and vcf
        self._test_protect('ci_vcf_fastq_parameters.yaml')
        self._test_ran_successfully([{'alignments': ['rna']}, 'binding_predictions', 'expression',
                                     'haplotyping', {'mutations': ['snpeffed', 'transgened']},
                                     'peptides', 'rankboost', 'reports'])

    def _test_protect(self, filename):
        protect_call = ['ProTECT',
                        '--config_file', os.path.join(self._projectRootPath(), 'src', 'protect',
                                                      'test', 'test_inputs',
                                                      filename),
                        '--workDir', self.work_dir,
                        os.path.join(self.test_dir, 'jobstore')]
        # Run ProTECT
        subprocess.check_call(protect_call)

    def _test_ran_successfully(self, expected_dirs):
        """
        Ensure that the test created the correct outputs.

        :param list expected_dirs: a list of dirs that the test was expected to create
        :return: True or False
        """
        results_dir = '/mnt/ephemeral/done/TEST'
        results_dir_contents = [(a, sorted(b), sorted(c)) for a, b, c in os.walk(results_dir)]
        results_dir_contents.sort()
        # Obtained by running [[(a, sorted(b), sorted(c)) for a, b, c in os.walk(results_dir)]
        # + .sort() on a successful run.
        contents_per_dir = {'alignments': {'normal_dna': ['normal_dna_fix_pg_sorted.bam',
                                                          'normal_dna_fix_pg_sorted.bam.bai'],
                                           'rna': ['rna_genome_sorted.bam',
                                                   'rna_genome_sorted.bam.bai',
                                                   'rna_transcriptome.bam'],
                                           'tumor_dna': ['tumor_dna_fix_pg_sorted.bam',
                                                         'tumor_dna_fix_pg_sorted.bam.bai']},
                            'binding_predictions': ('/mnt/ephemeral/done/TEST/binding_predictions',
                                                    [],
                                                    ['mhci_merged_files.list',
                                                     'mhcii_merged_files.list']),
                            'expression': ('/mnt/ephemeral/done/TEST/expression',
                                           [],
                                           ['rsem.genes.results', 'rsem.isoforms.results']),
                            'haplotyping': ('/mnt/ephemeral/done/TEST/haplotyping',
                                            [],
                                            ['mhci_alleles.list', 'mhcii_alleles.list']),
                            'mutations': {
                                'fusions': ('/mnt/ephemeral/done/TEST/mutations/fusions',
                                             [],
                                             ['fusion-inspector-predictions.tsv',
                                              'fusion.bedpe',
                                              'fusion.final',
                                              'rna_chimeric.junction',
                                              'star-fusion-predictions.tsv']),
                                'merged': ('/mnt/ephemeral/done/TEST/mutations/merged',
                                           [],
                                           ['all_merged.vcf', 'chr6.vcf']),
                                'mutect': ('/mnt/ephemeral/done/TEST/mutations/mutect',
                                           [],
                                           ['chr6.vcf']),
                                'muse': ('/mnt/ephemeral/done/TEST/mutations/muse',
                                         [],
                                         ['chr6.vcf']),
                                'radia': ('/mnt/ephemeral/done/TEST/mutations/radia',
                                          [],
                                          ['chr6.vcf']),
                                'somaticsniper': (
                                    '/mnt/ephemeral/done/TEST/mutations/somaticsniper',
                                    [],
                                    ['chr6.vcf']),
                                'strelka': {
                                    'indel': ('/mnt/ephemeral/done/TEST/mutations/strelka/indel',
                                              [],
                                              ['chr6.vcf']),
                                    'snv': ('/mnt/ephemeral/done/TEST/mutations/strelka/snv',
                                            [],
                                            ['chr6.vcf']),
                                },
                                'snpeffed': ('/mnt/ephemeral/done/TEST/mutations/snpeffed',
                                             [],
                                             ['mutations.vcf']),
                                'transgened': ('/mnt/ephemeral/done/TEST/mutations/transgened',
                                               [],
                                               ['mutations.vcf']),
                            },
                            'peptides': ('/mnt/ephemeral/done/TEST/peptides',
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
                            'rankboost': ('/mnt/ephemeral/done/TEST/rankboost',
                                          [],
                                          ['mhci_rankboost_concise_results.tsv',
                                           'mhci_rankboost_detailed_results.txt',
                                           'mhcii_rankboost_concise_results.tsv',
                                           'mhcii_rankboost_detailed_results.txt']),
                            'reports': ('/mnt/ephemeral/done/TEST/reports',
                                        [],
                                        ['car_t_target_report.txt',
                                         'immunotherapy_resistance_report.txt',
                                         'mhc_pathway_report.txt'])
                            }
        expected_contents = {}
        for dir in expected_dirs:
            if isinstance(dir, dict):
                assert len(dir.keys()) == 1
                if dir.keys()[0] == 'mutations':
                    expected_contents['mutations'] = ('/mnt/ephemeral/done/TEST/mutations',
                                                      sorted(dir['mutations']),
                                                      [])
                    for caller in dir['mutations']:
                        if caller == 'strelka':
                            expected_contents['mutations_strelka_1'] = (
                                '/mnt/ephemeral/done/TEST/mutations/strelka',
                                ['indel', 'snv'],
                                [])
                            expected_contents['mutations_strelka_2'] = \
                                contents_per_dir['mutations']['strelka']['indel']
                            expected_contents['mutations_strelka_3'] = \
                                contents_per_dir['mutations']['strelka']['snv']
                        else:
                            expected_contents['mutations_' + caller] = \
                                contents_per_dir['mutations'][caller]
                elif dir.keys()[0] == 'alignments':
                    alignment_files = []
                    for tissue_type in dir['alignments']:
                        alignment_files.extend(contents_per_dir['alignments'][tissue_type])
                    expected_contents['alignments'] = ('/mnt/ephemeral/done/TEST/alignments',
                                                       [],
                                                       sorted(alignment_files))
                else:
                    assert False
            else:
                expected_contents[dir] = contents_per_dir[dir]

        expected_outputs = [('/mnt/ephemeral/done/TEST',
                             sorted([x for x in expected_contents.keys()
                                      if not x.startswith('mutations_')]),
                             [])]
        expected_outputs.extend([expected_contents[d] for d in sorted(expected_contents.keys())])
        self.assertEqual(results_dir_contents, expected_outputs)
        shutil.rmtree(results_dir)