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

from protect.binding_prediction.mhci import predict_mhci_binding
from protect.binding_prediction.mhcii import predict_mhcii_binding
from protect.common import get_files_from_filestore, untargz, docker_path, export_results

import json
import os
import re


def spawn_antigen_predictors(job, transgened_files, phlat_files, univ_options, mhc_options):
    """
    Based on the number of alleles obtained from node 14, this module will spawn callers to predict
    MHCI:peptide and MHCII:peptide binding on the peptide files from node 17.  Once all MHC:peptide
    predictions are made, merge them via a follow-on job.

    ARGUMENTS
    1. transgened_files: REFER RETURN VALUE of run_transgene()
    2. phlat_files: REFER RETURN VALUE of merge_phlat_calls()
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. mhc_options: Dict of dicts of parameters specific to mhci and mhcii
                    respectively
         mhc_options
              |- 'mhci'
              |     |- 'method_file': <JSid for json file containing data
              |     |                  linking alleles, peptide lengths, and
              |     |                  prediction methods>
              |     +- 'pred': String describing prediction method to use
              +- 'mhcii'
                    |- 'method_file': <JSid for json file containing data
                    |                  linking alleles and prediction methods>
                    +- 'pred': String describing prediction method to use

    RETURN VALUES
    1. tuple of (mhci_preds, mhcii_preds)
         mhci_preds: Dict of return value from running predictions on a given
                     mhc for all peptides of length 9 and 10.
             mhci_preds
                |- <MHC molecule 1>_9_mer.faa: <PromisedJobReturnValue>
                |- <MHC molecule 1>_10_mer.faa: <PromisedJobReturnValue>
                |
                ..
                +- <MHC molecule n>_10_mer.faa: <PromisedJobReturnValue>
         mhcii_preds: Dict of return value from running predictions on a given
                     mhc for all peptides of length 15.
             mhci_preds
                |- <MHC molecule 1>_15_mer.faa: <PromisedJobReturnValue>
                |
                ..
                +- <MHC molecule n>_15_mer.faa: <PromisedJobReturnValue>

    This module corresponds to node 18 on the tree
    """
    job.fileStore.logToMaster('Running spawn_anti on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    mhci_options, mhcii_options = mhc_options
    pept_files = {
        '9_mer.faa': transgened_files['transgened_tumor_9_mer_snpeffed.faa'],
        '10_mer.faa': transgened_files['transgened_tumor_10_mer_snpeffed.faa'],
        '15_mer.faa': transgened_files['transgened_tumor_15_mer_snpeffed.faa']}
    input_files = {
        'mhci_alleles.list': phlat_files['mhci_alleles.list'],
        'mhcii_alleles.list': phlat_files['mhcii_alleles.list'],
        'mhci_restrictions.json.tar.gz': mhci_options['method_file'],
        'mhcii_restrictions.json.tar.gz': mhcii_options['method_file']}
    input_files = get_files_from_filestore(job, input_files, work_dir)
    for key in ('mhci_restrictions.json', 'mhcii_restrictions.json'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)

    mhci_alleles, mhcii_alleles = [], []
    with open(input_files['mhci_alleles.list'], 'r') as mhci_file:
        for line in mhci_file:
            mhci_alleles.append(line.strip())
    with open(input_files['mhcii_alleles.list'], 'r') as mhcii_file:
        for line in mhcii_file:
            mhcii_alleles.append(line.strip())
    # This file contains the list of allele:pept length combinations supported
    # by each prediction type.
    with open(input_files['mhci_restrictions.json'], 'r') as restfile:
        mhci_restrictions = json.load(restfile)
    with open(input_files['mhcii_restrictions.json'], 'r') as restfile:
        mhcii_restrictions = json.load(restfile)
    # Make a regexp to convert non alphanumeric characters in HLA names to _
    strip_allele_re = re.compile('[^A-Z0-9]')
    # For each mhci allele:peptfile combination, spawn a job and store the job handle in the dict.
    # Then do the same for mhcii
    mhci_preds, mhcii_preds = {}, {}
    for allele in mhci_alleles:
        stripped_allele = re.sub(strip_allele_re, '_', allele)
        for peptfile in ['9_mer.faa', '10_mer.faa']:
            peplen = peptfile.split('_')[0]
            # Ensure that the allele is among the list of accepted alleles
            try:
                if not mhci_restrictions[allele][peplen]:
                    continue
            except KeyError:
                continue
            predfile = ''.join([stripped_allele, '_', peptfile[:-4], '_mer.pred'])
            mhci_preds[predfile] = job.addChildJobFn(predict_mhci_binding, pept_files[peptfile],
                                                     allele, peplen, univ_options,
                                                     mhci_options, disk='10G').rv()
    for allele in mhcii_alleles:
        stripped_allele = re.sub(strip_allele_re, '_', allele)
        predfile = ''.join([stripped_allele, '_15_mer.pred'])
        if allele not in mhcii_restrictions[mhcii_options['pred']]:
            continue
        mhcii_preds[predfile] = job.addChildJobFn(predict_mhcii_binding, pept_files['15_mer.faa'],
                                                  allele, univ_options, mhcii_options,
                                                  disk='10G').rv()
    return mhci_preds, mhcii_preds


def merge_mhc_peptide_calls(job, antigen_predictions, transgened_files, univ_options):
    """
    This module will merge all the calls from nodes 18 and 19, and will filter for the top X%% of
    binders of each allele.  The module will then call the rank boosting script to finish off the
    pipeline.

    This module corresponds to node 19 on the tree
    """
    job.fileStore.logToMaster('Merging MHC calls')
    work_dir = os.getcwd()
    pept_files = {
        '10_mer.faa': transgened_files['transgened_tumor_10_mer_snpeffed.faa'],
        '10_mer.faa.map': transgened_files['transgened_tumor_10_mer_snpeffed.faa.map'],
        '15_mer.faa': transgened_files['transgened_tumor_15_mer_snpeffed.faa'],
        '15_mer.faa.map': transgened_files['transgened_tumor_15_mer_snpeffed.faa.map']}
    mhci_preds, mhcii_preds = antigen_predictions
    mhci_files = get_files_from_filestore(job, mhci_preds, work_dir)
    # First split mhcii_preds into prediction files and predictors and maintain keys so we can later
    # reference them in pairs
    mhcii_predictors = {x: y[1] for x, y in mhcii_preds.items()}
    mhcii_files = {x: y[0] for x, y in mhcii_preds.items()}
    mhcii_files = get_files_from_filestore(job, mhcii_files, work_dir)
    # Get peptide files
    pept_files = get_files_from_filestore(job, pept_files, work_dir)

    # Merge MHCI calls
    # Read 10-mer pepts into memory
    peptides = read_peptide_file(pept_files['10_mer.faa'])
    with open(pept_files['10_mer.faa.map'], 'r') as mapfile:
        pepmap = json.load(mapfile)
    # Incorporate peptide names into the merged calls
    with open('/'.join([work_dir, 'mhci_merged_files.list']), 'w') as mhci_resfile:
        for mhcifile in mhci_files.values():
            with open(mhcifile, 'r') as mf:
                for line in mf:
                    # Skip header lines
                    if not line.startswith('HLA'):
                        continue
                    line = line.strip().split('\t')
                    allele = line[0]
                    pept = line[5]
                    pred = line[7]
                    if float(pred) > 5.00:
                        continue
                    print_mhc_peptide((allele, pept, pred, pept), peptides, pepmap, mhci_resfile)
    # Merge MHCII calls
    # read 15-mer pepts into memory
    peptides = read_peptide_file(pept_files['15_mer.faa'])
    with open(pept_files['15_mer.faa.map'], 'r') as mapfile:
        pepmap = json.load(mapfile)
    # Incorporate peptide names into the merged calls
    with open('/'.join([work_dir, 'mhcii_merged_files.list']), 'w') as \
            mhcii_resfile:
        for mhciifile in mhcii_files.keys():
            core_col = None  # Variable to hold the column number with the core
            if mhcii_predictors[mhciifile] == 'Consensus':
                with open(mhcii_files[mhciifile], 'r') as mf:
                    for line in mf:
                        # Skip header lines
                        if not line.startswith('HLA'):
                            continue
                        line = line.strip().split('\t')
                        allele = line[0]
                        pept = line[4]
                        pred = line[6]
                        if core_col:
                            core = line[core_col] if core_col else 'NOCORE'
                        else:
                            methods = line[5].lstrip('Consensus(').rstrip(')')
                            methods = methods.split(',')
                            if 'NN' in methods:
                                core_col = 13
                            elif 'netMHCIIpan' in methods:
                                core_col = 17
                            elif 'Sturniolo' in methods:
                                core_col = 19
                            elif 'SMM' in methods:
                                core_col = 10
                            core = line[core_col] if core_col else 'NOCORE'
                        if float(pred) > 5.00:
                            continue
                        print_mhc_peptide((allele, pept, pred, core), peptides, pepmap,
                                          mhcii_resfile)
            elif mhcii_predictors[mhciifile] == 'Sturniolo':
                with open(mhcii_files[mhciifile], 'r') as mf:
                    for line in mf:
                        # Skip header lines
                        if not line.startswith('HLA'):
                            continue
                        line = line.strip().split('\t')
                        allele = line[0]
                        pept = line[5]
                        pred = line[6]
                        core = line[19]  #
                        if float(pred) > 5.00:
                            continue
                        print_mhc_peptide((allele, pept, pred, core), peptides, pepmap,
                                          mhcii_resfile)
            elif mhcii_predictors[mhciifile] == 'netMHCIIpan':
                with open(mhcii_files[mhciifile], 'r') as mf:
                    # Get the allele from the first line and skip the second line
                    allele = re.sub('-DQB', '/DQB', mf.readline().strip())
                    _ = mf.readline()
                    for line in mf:
                        line = line.strip().split('\t')
                        pept = line[1]
                        pred = line[5]
                        core = 'NOCORE'
                        peptide_name = line[2]
                        if float(pred) > 5.00:
                            continue
                        print(allele, pept, peptide_name, core, '0', pred, pepmap[peptide_name],
                              sep='\t', file=mhcii_resfile)
            else:
                raise RuntimeError('Shouldn\'t ever see this!!!')
    output_files = defaultdict()
    for mhc_file in [mhci_resfile.name, mhcii_resfile.name]:
        output_files[os.path.split(mhc_file)[1]] = job.fileStore.writeGlobalFile(mhc_file)
        export_results(job, mhc_file, univ_options, subfolder='binding_predictions')
    return output_files


def print_mhc_peptide(neoepitope_info, peptides, pepmap, outfile):
    """
    To reduce code redundancy, this module will accept data from merge_mhc_peptide_calls for a given
    neoepitope and print it to outfile
    ARGUMENTS
    1. neoepitope_info: Tuple of (<allele>, <peptide_sequence>,
                                  <binding_prediction>)
    2. peptides: Dict of all IARS considered
           peptides
              |- 'neoepitope_1': <peptide_sequence>
              ..
              |- 'neoepitope_n': <peptide_sequence>
    3. pepmap: Info correlating neoepitope with the gene and transcript level
               mutations.
           peptides
              |- 'neoepitope_1':
              |      'ensembl_gene\thugo_gene\tcomma_sep_transcript_mutations'
              ..
              +- 'neoepitope_n':
                     'ensembl_gene\thugo_gene\tcomma_sep_transcript_mutations'

    """
    allele, pept, pred, core = neoepitope_info
    peptide_names = [x for x, y in peptides.items() if pept in y]
    # For each peptide, append the ensembl gene
    for peptide_name in peptide_names:
        print(allele, pept, peptide_name, core, '0', pred, pepmap[peptide_name], sep='\t',
              file=outfile)
    return None


def read_peptide_file(in_peptfile):
    """
    This module reads an input peptide fasta file into memory in the form of a dict with
    key = fasta record name
    value = corresponding peptide sequence
    """
    peptides = defaultdict()
    with open(in_peptfile, 'r') as peptfile:
        for line in peptfile:
            if line.startswith('>'):
                pept = line.strip().lstrip('>')
                peptides[pept] = ''
            else:
                peptides[pept] = line.strip()
    return peptides
