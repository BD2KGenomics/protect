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

from collections import Counter

from protect.common import export_results, get_files_from_filestore, untargz
from protect.haplotyping.phlat import parse_phlat_file

import json
import os


def run_mhc_gene_assessment(job, rsem_files, rna_haplotype, univ_options, mhc_genes_options):
    """
    A wrapper for assess_mhc_genes.

    :param dict rsem_files: Results form running rsem
    :param str rna_haplotype: The job store id for the rna haplotype file
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict mhc_genes_options: Options specific to assessing the MHC genes
    :return: The results of running assess_mhc_genes
    :rtype: toil.fileStore.FileID
    """
    return job.addChildJobFn(assess_mhc_genes, rsem_files['rsem.isoforms.results'], rna_haplotype,
                             univ_options, mhc_genes_options).rv()


def assess_mhc_genes(job, isoform_expression, rna_haplotype, univ_options, mhc_genes_options):
    """
    Assess the prevalence of the various genes in the MHC pathway and return a report in the tsv
    format.

    :param toil.fileStore.FileID isoform_expression: fsID for the rsem isoform expression file
    :param toil.fileStore.FileID rna_haplotype: fsID for the RNA PHLAT file
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict mhc_genes_options: Options specific to assessing the MHC genes
    :return: The fsID for the mhc pathway report file
    :rtype: toil.fileStore.FileID
    """
    job.fileStore.logToMaster('Running mhc gene assessment on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'rsem_quant.tsv': isoform_expression,
        'rna_haplotype.sum': rna_haplotype,
        'mhc_genes.json.tar.gz': mhc_genes_options['genes_file']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    input_files['mhc_genes.json'] = untargz(input_files['mhc_genes.json.tar.gz'], work_dir)

    # Read in the MHC genes
    with open(input_files['mhc_genes.json']) as mhc_file:
        mhc_genes = json.load(mhc_file)

    # Parse the rna phlat file
    with open(input_files['rna_haplotype.sum']) as rna_mhc:
        mhc_alleles = {'HLA_A': [], 'HLA_B': [], 'HLA_C': [], 'HLA_DPA': [], 'HLA_DQA': [],
                       'HLA_DPB': [], 'HLA_DQB': [], 'HLA_DRB': []}
        mhc_alleles = parse_phlat_file(rna_mhc, mhc_alleles)

    # Process the isoform expressions
    gene_expressions = Counter()
    with open(input_files['rsem_quant.tsv']) as rsem_file:
        line = rsem_file.readline()
        line = line.strip().split()
        assert line == ['transcript_id', 'gene_id', 'length', 'effective_length', 'expected_count',
                        'TPM', 'FPKM', 'IsoPct']
        for line in rsem_file:
            line = line.strip().split()
            gene_expressions[line[1]] += float(line[5])

    with open(os.path.join(work_dir, 'mhc_pathway_report.txt'), 'w') as mpr:
        for section in mhc_genes:
            print(section.center(48, ' '), file=mpr)
            print("{:12}{:12}{:12}{:12}".format("Gene", "Threshold", "Observed", "Result"),
                  file=mpr)
            if section == 'MHCI loading':
                for mhci_allele in 'HLA_A', 'HLA_B', 'HLA_C':
                    num_alleles = len(mhc_alleles[mhci_allele])
                    result = 'FAIL' if num_alleles == 0 else 'LOW' if num_alleles == 1 else 'PASS'
                    print("{:12}{:<12}{:<12}{:12}".format(mhci_allele, 2, num_alleles, result),
                          file=mpr)
            elif section == 'MHCII loading':
                # TODO DP alleles
                for mhcii_allele in ('HLA_DQA', 'HLA_DQB', 'HLA_DRA', 'HLA_DRB'):
                    if mhcii_allele != 'HLA_DRA':
                        num_alleles = len(mhc_alleles[mhcii_allele])
                        result = ('FAIL' if num_alleles == 0 else
                                  'LOW' if num_alleles == 1 else
                                  'PASS')
                        print("{:12}{:<12}{:<12}{:12}".format(mhcii_allele, 2, num_alleles, result),
                              file=mpr)
                    else:
                        # FIXME This is hardcoded for now. We need to change this.
                        result = 'LOW' if gene_expressions['ENSG00000204287.9'] <= 69.37 else 'PASS'
                        print("{:12}{:<12}{:<12}{:12}".format(
                                    'HLA_DRA', gene_expressions['ENSG00000204287.9'], '69.37',
                                    result), file=mpr)
            for gene, ensgene, first_quart in mhc_genes[section]:
                result = 'LOW' if gene_expressions[ensgene] <= float(first_quart) else 'PASS'
                print("{:12}{:<12}{:<12}{:12}".format(gene, float(first_quart),
                                                      gene_expressions[ensgene], result), file=mpr)
            print('', file=mpr)
    output_file = job.fileStore.writeGlobalFile(mpr.name)
    export_results(job, output_file, mpr.name, univ_options, subfolder='reports')
    return output_file
