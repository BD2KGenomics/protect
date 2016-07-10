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
from __future__ import absolute_import, print_function
from collections import defaultdict
from math import ceil
from protect.common import (docker_call, get_files_from_filestore, is_gzipfile, untargz,
                            docker_path, export_results)

import os
import re
import sys


# disk for phlat
def phlat_disk(rna_fastqs):
    return ceil(sum([f.size for f in rna_fastqs]) + 524288) + 8053063680


def run_phlat(job, fastqs, sample_type, univ_options, phlat_options):
    """
    This module will run PHLAT on SAMPLE_TYPE fastqs.

    ARGUMENTS -- <ST> depicts the sample type. Substitute with 'tumor_dna',
                 'normal_dna', or 'tumor_rna'
    1. fastqs: Dict of list of input WGS/WXS fastqs
         fastqs
              +- '<ST>': [<JSid for 1.fastq> , <JSid for 2.fastq>]
    2. sample_type: string of 'tumor' or 'normal'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. phlat_options: Dict of parameters specific to phlat
         phlat_options
              |- 'tool_index': <JSid for the PHLAT index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_file: <JSid for the allele predictions for ST>

    This module corresponds to nodes 5, 6 and 7 on the tree
    """
    job.fileStore.logToMaster('Running phlat on %s:%s' % (univ_options['patient'], sample_type))
    print(phlat_options, file=sys.stderr)
    work_dir = os.getcwd()
    input_files = {
        'input_1.fastq': fastqs[0],
        'input_2.fastq': fastqs[1],
        'phlat_index.tar.gz': phlat_options['tool_index']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    # Handle gzipped files
    gz = '.gz' if is_gzipfile(input_files['input_1.fastq']) else ''
    if gz:
        for read_file in 'input_1.fastq', 'input_2.fastq':
            os.symlink(read_file, read_file + gz)
            input_files[read_file + gz] = input_files[read_file] + gz
    # Untar the index
    input_files['phlat_index'] = untargz(input_files['phlat_index.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['-1', input_files['input_1.fastq' + gz],
                  '-2', input_files['input_2.fastq' + gz],
                  '-index', input_files['phlat_index'],
                  '-b2url', '/usr/local/bin/bowtie2',
                  '-tag', sample_type,
                  '-e', '/home/phlat-1.0',  # Phlat directory home
                  '-o', '/data',  # Output directory
                  '-p', str(phlat_options['n'])]  # Number of threads
    docker_call(tool='phlat', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_file = job.fileStore.writeGlobalFile(''.join([work_dir, '/', sample_type, '_HLA.sum']))
    return output_file


def merge_phlat_calls(job, tumor_phlat, normal_phlat, rna_phlat, univ_options):
    """
    This module will merge the results form running PHLAT on the 3 input fastq
    pairs.

    ARGUMENTS
    1. tumor_phlat: <JSid for tumor DNA called alleles>
    2. normal_phlat: <JSid for normal DNA called alleles>
    3. rna_phlat: <JSid for tumor RNA called alleles>

    RETURN VALUES
    1. output_files: Dict of JSids for consensus MHCI and MHCII alleles
             output_files
                    |- 'mhci_alleles.list': <JSid>
                    +- 'mhcii_alleles.list': <JSid>

    This module corresponds to node 14 on the tree
    """
    job.fileStore.logToMaster('Merging Phlat calls')
    work_dir = os.getcwd()
    input_files = {
        'tumor_dna': tumor_phlat,
        'normal_dna': normal_phlat,
        'tumor_rna': rna_phlat}
    input_files = get_files_from_filestore(job, input_files, work_dir)
    with open(input_files['tumor_dna'], 'r') as td_file, \
            open(input_files['normal_dna'], 'r') as nd_file, \
            open(input_files['tumor_rna'], 'r') as tr_file:
        # TODO: Could this be a defautdict?
        mhc_alleles = {'HLA_A': [], 'HLA_B': [], 'HLA_C': [], 'HLA_DPA': [], 'HLA_DQA': [],
                       'HLA_DPB': [], 'HLA_DQB': [], 'HLA_DRB': []}
        for phlatfile in td_file, nd_file, tr_file:
            mhc_alleles = parse_phlat_file(phlatfile, mhc_alleles)
    # Get most probable alleles for each allele group and print to output
    with open(os.path.join(work_dir, 'mhci_alleles.list'), 'w') as mhci_file, \
            open(os.path.join(work_dir, 'mhcii_alleles.list'), 'w') as mhcii_file:
        for mhci_group in ['HLA_A', 'HLA_B', 'HLA_C']:
            mpa = most_probable_alleles(mhc_alleles[mhci_group])
            print('\n'.join([''.join(['HLA-', x]) for x in mpa]), file=mhci_file)
        drb_mpa = most_probable_alleles(mhc_alleles['HLA_DRB'])
        print('\n'.join([''.join(['HLA-', x]) for x in drb_mpa]), file=mhcii_file)
        dqa_mpa = most_probable_alleles(mhc_alleles['HLA_DQA'])
        dqb_mpa = most_probable_alleles(mhc_alleles['HLA_DQB'])
        for dqa_allele in dqa_mpa:
            for dqb_allele in dqb_mpa:
                print(''.join(['HLA-', dqa_allele, '/', dqb_allele]), file=mhcii_file)
    output_files = defaultdict()
    for allele_file in ['mhci_alleles.list', 'mhcii_alleles.list']:
        output_files[allele_file] = job.fileStore.writeGlobalFile(os.path.join(work_dir,
                                                                               allele_file))
        export_results(job, os.path.join(work_dir, allele_file), univ_options,
                       subfolder='haplotyping')
    return output_files


def parse_phlat_file(phlatfile, mhc_alleles):
    """
    Parse the input phlat file to pull out the alleles it contains
    :param phlatfile: Open file descriptor for a phlat output sum file
    :param mhc_alleles: dictionary of alleles.
    """
    for line in phlatfile:
        if line.startswith('Locus'):
            continue
        line = line.strip().split()
        if line[0].startswith('HLA_D'):
            line[0] = line[0][:-1]  # strip the last character
        # Sometimes we get an error saying there was insufficient read
        # converage. We need to drop that line.
        # E.g. HLA_DQB1 no call due to insufficient reads at this locus
        if line[1] == 'no':
            continue
        if line[4] != 'NA':
            split_field = line[1].split(':')
            if len(split_field) >= 2 and not split_field[1] == 'xx':
                mhc_alleles[line[0]].append((line[1], line[4]))
        if line[5] != 'NA':
            split_field = line[2].split(':')
            if len(split_field) >= 2 and not split_field[1] == 'xx':
                mhc_alleles[line[0]].append((line[2], line[5]))
    return mhc_alleles


def most_probable_alleles(allele_list):
    """
    This module accepts a list of tuples of (allele, p_value) pairs. It returns the 2 most probable
    alleles for that group.
    """
    all_alleles = defaultdict()
    # First collect all the keys.  Make a dict with allele as key and list of pvalues as value
    for allele, pvalue in allele_list:
        allele = re.split(':', allele)
        # Ensure allele has enough resolution for mhc:peptide prediciton.
        # HLA-A*02:01:04  -> ['HLA-A*02', '01', '04']   => At least 2 fields are required for
        #                                                  satisfying criteria.
        if len(allele) < 2:
            continue
        allele = ':'.join([allele[0], allele[1]])  # stitch back together
        try:
            all_alleles[allele].append(float(pvalue))
        except KeyError:
            all_alleles[allele] = [float(pvalue)]
    # If there are less than 2 alleles, report all
    if len(all_alleles.keys()) <= 2:
        return all_alleles.keys()
    # Else, get the two with most evidence.  Evidence is gauged by
    # a) How many files (of the 3) thought that Allele was present
    # b) In a tie, who has a lower avg p value
    # In the lambda function, if 2 alleles have the same number of calls, the sum of the p values is
    # a measure of the avg because avg = sum / n and n is equal in both of them.
    else:
        return sorted(all_alleles.keys(),
                      key=lambda x: (-len(all_alleles[x]), sum(all_alleles[x])))[0:2]
