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
from protect.binding_prediction.mhcii import predict_mhcii_binding, predict_netmhcii_binding
from protect.common import get_files_from_filestore, export_results, read_peptide_file, untargz

import json
import os
import pandas
import re


def spawn_antigen_predictors(job, transgened_files, phlat_files, univ_options, mhc_options):
    """
    Spawn a job to predict MHCI:peptide and MHCII:peptide binding on the input peptide file for each
    allele in in the input haplotype (phlat) files.

    :param dict transgened_files: Dict of tumor and normal peptide fsIDs and the tumor .map fsIDs
    :param dict phlat_files: Dict of MHCI and MHCII haplotypes
    :param dict univ_options: Dict of universal options used by almost all tools
    :param tuple mhc_options: Options specific to mhci and mhcii binding predictions
    :return: Tuple of dicts of mhci and mhcii predictions
                (mhci_preds, mhcii_preds)
                      |          |- <allele>:
                      |          |     |- 'tumor': fsID
                      |          |     +- 'normal': (fsID, str)
                      |          |
                      |          |- ...
                      |          |
                      |          +- <allele>:
                      |                |- 'tumor': fsID
                      |                +- 'normal': (fsID, str)
                      |
                      |- <allele>:
                      |     |- 'tumor': fsID
                      |     +- 'normal': fsID
                      |
                      |- ...
                      |
                      +- <allele>:
                            |- 'tumor': fsID
                            +- 'normal': fsID
    :rtype: tuple(dict, dict)
    """
    work_dir = os.getcwd()
    mhci_options, mhcii_options = mhc_options
    pept_files = {
        'T_9_mer.faa': transgened_files['transgened_tumor_9_mer_snpeffed.faa'],
        'N_9_mer.faa': transgened_files['transgened_normal_9_mer_snpeffed.faa'],
        'T_10_mer.faa': transgened_files['transgened_tumor_10_mer_snpeffed.faa'],
        'N_10_mer.faa': transgened_files['transgened_normal_10_mer_snpeffed.faa'],
        'T_15_mer.faa': transgened_files['transgened_tumor_15_mer_snpeffed.faa'],
        'N_15_mer.faa': transgened_files['transgened_normal_15_mer_snpeffed.faa']
    }
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
    # For each mhci allele:peptfile combination, spawn a job and store the job handle in the dict.
    # Then do the same for mhcii
    mhci_preds, mhcii_preds = {}, {}
    for allele in mhci_alleles:
        for peplen in ['9', '10']:
            peptfile = peplen + '_mer.faa'
            # Ensure that the allele is among the list of accepted alleles
            try:
                if not mhci_restrictions[allele][peplen]:
                    continue
            except KeyError:
                continue
            mhci_job = job.addChildJobFn(predict_mhci_binding, pept_files['T_' + peptfile], allele,
                                         peplen, univ_options, mhci_options, disk='100M',
                                         memory='100M', cores=1)
            mhci_preds[(allele, peplen)] = mhci_job.addChildJobFn(
                predict_normal_binding,
                mhci_job.rv(),
                {x: y for x, y in pept_files.items() if peplen in x},
                allele,
                peplen,
                univ_options,
                mhci_options,
                disk='100M',
                memory='100M',
                cores=1).rv()
    for allele in mhcii_alleles:
        if allele not in mhcii_restrictions[mhcii_options['pred']]:
            continue
        mhcii_job = job.addChildJobFn(predict_mhcii_binding, pept_files['T_15_mer.faa'],  allele,
                                      univ_options, mhcii_options, disk='100M', memory='100M',
                                      cores=1)
        mhcii_preds[(allele, 15)] = mhcii_job.addFollowOnJobFn(
            predict_normal_binding,
            mhcii_job.rv(),
            {x: y for x, y in pept_files.items() if '15' in x},
            allele,
            '15',
            univ_options,
            mhcii_options,
            disk='100M',
            memory='100M',
            cores=1).rv()
    job.fileStore.logToMaster('Ran spawn_anti on %s successfully' % univ_options['patient'])
    return mhci_preds, mhcii_preds


def read_fastas(input_files):
    """
    Read the tumor and normal fastas into a joint dict.

    :param dict input_files: A dict containing filename: filepath for T_ and N_ transgened files.
    :return: The read fastas in a dictionary of tuples
    :rtype: dict
    """
    tumor_file = [y for x, y in input_files.items() if x.startswith('T')][0]
    normal_file = [y for x, y in input_files.items() if x.startswith('N')][0]
    output_files = defaultdict(list)
    output_files = _read_fasta(tumor_file, output_files)
    num_entries = len(output_files)
    output_files = _read_fasta(normal_file, output_files)
    assert len(output_files) == num_entries
    return output_files


def _read_fasta(fasta_file, output_dict):
    """
    Read the peptide fasta into an existing dict.

    :param str fasta_file: The peptide file
    :param dict output_dict: The dict to appends results to.
    :return: output_dict
    :rtype: dict
    """
    read_name = None
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                read_name = line.lstrip('>')
            else:
                assert read_name is not None, line
                output_dict[read_name].append(line.strip())
    return output_dict


def _process_consensus_mhcii(mhc_file, normal=False):
    """
    Process the results from running IEDB MHCII binding predictions using the consensus method into
    a pandas dataframe.

    :param str mhc_file: Output file containing consensus mhcii:peptide binding predictions
    :param bool normal: Is this processing the results of a normal?
    :return: Results in a tabular format
    :rtype: pandas.DataFrame
    """
    core_col = None  # Variable to hold the column number with the core
    results = pandas.DataFrame(columns=['allele', 'pept', 'tumor_pred', 'core'])
    with open(mhc_file, 'r') as mf:
        peptides = set()
        for line in mf:
            # Skip header lines
            if not line.startswith('HLA'):
                continue
            line = line.strip().split('\t')
            allele = line[0]
            pept = line[4]
            pred = line[6]
            if core_col:
                core = line[core_col]
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
            if float(pred) > 5.00 and not normal:
                continue
            results.loc[len(results)] = [allele, pept, pred, core]
    results.drop_duplicates(inplace=True)
    return results


def _process_sturniolo_mhcii(mhc_file, normal=False):
    """
    Process the results from running IEDB MHCII binding predictions using the sturniolo method into
    a pandas dataframe.

    :param str mhc_file: Output file containing sturniolo mhcii:peptide binding predictions
    :param bool normal: Is this processing the results of a normal?
    :return: Results in a tabular format
    :rtype: pandas.DataFrame
    """
    results = pandas.DataFrame(columns=['allele', 'pept', 'tumor_pred', 'core'])
    with open(mhc_file, 'r') as mf:
        peptides = set()
        for line in mf:
            # Skip header lines
            if not line.startswith('HLA'):
                continue
            line = line.strip().split('\t')
            allele = line[0]
            pept = line[4]
            pred = line[6]
            core = line[19]
            if float(pred) > 5.00 and not normal:
                continue
            results.loc[len(results)] = [allele, pept, pred, core]
    results.drop_duplicates(inplace=True)
    return results


def _process_net_mhcii(mhc_file, normal=False):
    """
    Process the results from running NetMHCIIpan binding predictions into a pandas dataframe.

    :param str mhc_file: Output file containing netmhciipan mhcii:peptide binding predictions
    :param bool normal: Is this processing the results of a normal?
    :return: Results in a tabular format
    :rtype: pandas.DataFrame
    """
    results = pandas.DataFrame(columns=['allele', 'pept', 'tumor_pred', 'core', 'peptide_name'])
    with open(mhc_file, 'r') as mf:
        peptides = set()
        # Get the allele from the first line and skip the second line
        allele = re.sub('-DQB', '/DQB', mf.readline().strip())
        _ = mf.readline()
        for line in mf:
            line = line.strip().split('\t')
            pept = line[1]
            pred = line[5]
            core = 'NOCORE'
            peptide_name = line[2]
            if float(pred) > 5.00 and not normal:
                continue
            results.loc[len(results)] = [allele, pept, pred, core, peptide_name]
    results.drop_duplicates(inplace=True)
    return results


def _process_mhci(mhc_file, normal=False):
    """
    Process the results from running IEDB MHCI binding predictions into a pandas dataframe.

    :param str mhc_file: Output file containing netmhciipan mhci:peptide binding predictions
    :param bool normal: Is this processing the results of a normal?
    :return: Results in a tabular format
    :rtype: pandas.DataFrame
    """
    results = pandas.DataFrame(columns=['allele', 'pept', 'tumor_pred', 'core'])
    with open(mhc_file, 'r') as mf:
        peptides = set()
        for line in mf:
            # Skip header lines
            if not line.startswith('HLA'):
                continue
            line = line.strip().split('\t')
            allele = line[0]
            pept = line[5]
            pred = line[7]
            if float(pred) > 5.00 and not normal:
                continue
            results.loc[len(results)] = [allele, pept, pred, pept]
    results.drop_duplicates(inplace=True)
    return results


def pept_diff(p1, p2):
    """
    Return the number of differences betweeen 2 peptides

    :param str p1: Peptide 1
    :param str p2: Peptide 2
    :return: The number of differences between the pepetides
    :rtype: int

    >>> pept_diff('ABCDE', 'ABCDF')
    1
    >>> pept_diff('ABCDE', 'ABDFE')
    2
    >>> pept_diff('ABCDE', 'EDCBA')
    4
    >>> pept_diff('ABCDE', 'ABCDE')
    0
    """
    return sum([p1[i] != p2[i] for i in range(len(p1))])


def _get_normal_peptides(mhc_df, iars, peplen):
    """
    Get the corresponding normal peptides for the tumor peptides that have already been subjected to
    mhc:peptide binding prediction.

    :param pandas.DataFrame mhc_df: The dataframe of mhc:peptide binding results
    :param dict iars: The dict of lists of tumor and normal peptide iar sequences
    :param str peplen: Length of the peptides to consider.
    :return: normal peptides and the updated results containing the normal peptides
    :rtype: tuple(pandas.DataFrame, list)
    """
    peplen = int(peplen)
    normal_peptides = []
    for pred in mhc_df.itertuples():
        containing_iars = [i for i, sl in iars.items() if pred.pept in sl[0]]
        assert len(containing_iars) != 0, "No IARS contained the peptide"
        if len(containing_iars) > 1:
            # If there are multiple IARs, they all or none of them have to have a corresponding
            # normal.
            assert len(set([len(y) for x, y in iars.items() if x in containing_iars])) == 1
        if len(iars[containing_iars[0]]) == 1:
            # This is a fusion and has no corresponding normal
            normal_peptides.append('N' * peplen)
        else:
            tum, norm = iars[containing_iars[0]]
            pos = tum.find(pred.pept)
            temp_normal_pept = norm[pos:pos + peplen]
            if pept_diff(pred.pept, temp_normal_pept) == 1:
                normal_peptides.append(norm[pos:pos + peplen])
            else:
                if len(tum) == len(norm):
                    # Too (2+) many single nucleotide changes to warrant having a normal counterpart
                    normal_peptides.append('N' * peplen)
                else:
                    # There is an indel in play. The difference cannot be in the last AA as that
                    # would have come out properly in the first case. There is a possibility that
                    # the indel was in the first AA causing a shift. We can handle that by looking
                    # at the suffix.
                    pos = norm.find(pred.pept[1:])
                    if pos != -1:
                        # The suffix was found,
                        normal_peptides.append(norm[pos-1:pos + peplen])
                    else:
                        # The indel was too large to warrant having a normal counterpart
                        normal_peptides.append('N' * peplen)
    mhc_df['normal_pept'] = normal_peptides
    return mhc_df, normal_peptides


def predict_normal_binding(job, binding_result, transgened_files, allele, peplen, univ_options,
                           mhc_options):
    """
    Predict the binding score for the normal counterparts of the peptides in mhc_dict and then
    return the results in a properly formatted structure.

    :param str binding_result: The results from running predict_mhci_binding or
           predict_mhcii_binding on a single allele
    :param dict transgened_files: A dictionary containing the jobstore IDs for "T_<peplen>_mer.faa"
           and "N_<peplen>_mer.faa"
    :param str allele: The allele to get binding for
    :param str peplen: The peptide length
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict mhc_options: Options specific to mhci or mhcii binding predictions
    :return: A fully filled out mhc_dict with normal information
             output_dict:
                |- 'tumor': fsID
                +- 'normal': fsID or (fsID, str)     -- Depending on MHCI or MHCII
    :rtype: dict
    """
    work_dir = os.getcwd()
    results = pandas.DataFrame(columns=['allele', 'pept', 'tumor_pred', 'core'])
    input_files = get_files_from_filestore(job, transgened_files, work_dir)
    iars = read_fastas(input_files)
    if peplen == '15':  # MHCII
        mhc_file = job.fileStore.readGlobalFile(binding_result[0],
                                                os.path.join(work_dir, 'mhci_results'))
        predictor = binding_result[1]
        core_col = None  # Variable to hold the column number with the core
        if predictor is None:
            return {'tumor': None,
                    'normal': None,
                    'predictor': None}
        elif predictor == 'Consensus':
            results = _process_consensus_mhcii(mhc_file)
            results, peptides = _get_normal_peptides(results, iars, peplen)
            with open('peptides.faa', 'w') as pfile:
                for pept in peptides:
                    print('>', pept, '\n', pept, sep='', file=pfile)
            peptfile = job.fileStore.writeGlobalFile(pfile.name)
            with open('results.json', 'w') as rj:
                json.dump(results.to_json(), rj)
            return {'tumor': job.fileStore.writeGlobalFile(rj.name),
                    'normal': job.addChildJobFn(predict_mhcii_binding, peptfile, allele,
                                                univ_options, mhc_options, disk='100M',
                                                memory='100M', cores=1).rv(),
                    'predictor': 'Consensus'}
        elif predictor == 'Sturniolo':
            results = _process_sturniolo_mhcii(mhc_file)
            results, peptides = _get_normal_peptides(results, iars, peplen)
            with open('peptides.faa', 'w') as pfile:
                for pept in peptides:
                    print('>', pept, '\n', pept, sep='', file=pfile)
            peptfile = job.fileStore.writeGlobalFile(pfile.name)
            with open('results.json', 'w') as rj:
                json.dump(results.to_json(), rj)
            return {'tumor': job.fileStore.writeGlobalFile(rj.name),
                    'normal': job.addChildJobFn(predict_mhcii_binding, peptfile, allele,
                                                univ_options, mhc_options, disk='100M',
                                                memory='100M', cores=1).rv(),
                    'predictor': 'Sturniolo'}
        elif predictor == 'netMHCIIpan':
            results = _process_net_mhcii(mhc_file)
            results, peptides = _get_normal_peptides(results, iars, peplen)
            with open('peptides.faa', 'w') as pfile:
                for pept in peptides:
                    print('>', pept, '\n', pept, sep='', file=pfile)
            peptfile = job.fileStore.writeGlobalFile(pfile.name)
            with open('results.json', 'w') as rj:
                json.dump(results.to_json(), rj)
            return {'tumor': job.fileStore.writeGlobalFile(rj.name),
                    'normal': job.addChildJobFn(predict_netmhcii_binding, peptfile, allele,
                                                univ_options, mhc_options['netmhciipan'],
                                                disk='100M', memory='100M',
                                                cores=1).rv(),
                    'predictor': 'netMHCIIpan'}
        else:
            raise RuntimeError('Shouldn\'t ever see this!!!')

    else:  # MHCI
        mhc_file = job.fileStore.readGlobalFile(binding_result,
                                                os.path.join(work_dir, 'mhci_results'))
        results = _process_mhci(mhc_file)
        results, peptides = _get_normal_peptides(results, iars, peplen)
        with open('peptides.faa', 'w') as pfile:
            for pept in peptides:
                print('>', pept, '\n', pept, sep='', file=pfile)
        peptfile = job.fileStore.writeGlobalFile(pfile.name)
        with open('results.json', 'w') as rj:
            json.dump(results.to_json(), rj)
        job.fileStore.logToMaster('Ran predict_normal_binding on %s for allele %s and length %s '
                                  'successfully' % (univ_options['patient'], allele, peplen))
        return {'tumor': job.fileStore.writeGlobalFile(rj.name),
                'normal': job.addChildJobFn(predict_mhci_binding, peptfile, allele, peplen,
                                            univ_options, mhc_options, disk='100M', memory='100M',
                                            cores=1).rv()}


def merge_mhc_peptide_calls(job, antigen_predictions, transgened_files, univ_options):
    """
    Merge all the calls generated by spawn_antigen_predictors.

    :param dict antigen_predictions: The return value from running :meth:`spawn_antigen_predictors`
    :param dict transgened_files: The transgened peptide files
    :param dict univ_options: Universal options for ProTECT
    :return: merged binding predictions
             output_files:
                 |- 'mhcii_merged_files.list': fsID
                 +- 'mhci_merged_files.list': fsID
    :rtype: dict
    """
    job.fileStore.logToMaster('Merging MHC calls')
    work_dir = os.getcwd()
    pept_files = {
        '10_mer.faa': transgened_files['transgened_tumor_10_mer_snpeffed.faa'],
        '10_mer.faa.map': transgened_files['transgened_tumor_10_mer_snpeffed.faa.map'],
        '15_mer.faa': transgened_files['transgened_tumor_15_mer_snpeffed.faa'],
        '15_mer.faa.map': transgened_files['transgened_tumor_15_mer_snpeffed.faa.map']}
    pept_files = get_files_from_filestore(job, pept_files, work_dir)
    mhci_preds, mhcii_preds = antigen_predictions

    mhci_called = mhcii_called = False
    # Merge MHCI calls
    # Read 10-mer pepts into memory
    peptides = read_peptide_file(pept_files['10_mer.faa'])
    with open(pept_files['10_mer.faa.map'], 'r') as mapfile:
        pepmap = json.load(mapfile)
    with open('/'.join([work_dir, 'mhci_merged_files.list']), 'w') as mhci_resfile:
        for key in mhci_preds:
            tumor_file = job.fileStore.readGlobalFile(mhci_preds[key]['tumor'])
            with open(tumor_file) as t_f:
                tumor_df = pandas.read_json(eval(t_f.read()))
            if tumor_df.empty:
                continue
            mhci_called = True
            # TODO: There must be a better way of doing this
            normal_df = _process_mhci(job.fileStore.readGlobalFile(mhci_preds[key]['normal']),
                                      normal=True)
            normal_dict = normal_df.set_index('pept')['tumor_pred']
            normal_preds = [normal_dict[x] for x in list(tumor_df['normal_pept'])]
            tumor_df['normal_pred'] = normal_preds
            for pred in tumor_df.itertuples():
                print_mhc_peptide(pred, peptides, pepmap, mhci_resfile)
    # Merge MHCII calls
    # read 15-mer pepts into memory
    peptides = read_peptide_file(pept_files['15_mer.faa'])
    with open(pept_files['15_mer.faa.map'], 'r') as mapfile:
        pepmap = json.load(mapfile)
    # Incorporate peptide names into the merged calls
    with open('/'.join([work_dir, 'mhcii_merged_files.list']), 'w') as \
            mhcii_resfile:
        for key in mhcii_preds:
            if mhcii_preds[key]['predictor'] is None:
                continue
            mhcii_called = True
            tumor_file = job.fileStore.readGlobalFile(mhcii_preds[key]['tumor'])
            with open(tumor_file) as t_f:
                tumor_df = pandas.read_json(eval(t_f.read()))
            if tumor_df.empty:
                continue
            # TODO: There must be a better way of doing this
            if mhcii_preds[key]['predictor'] == 'Consensus':
                normal_df = _process_consensus_mhcii(
                    job.fileStore.readGlobalFile(mhcii_preds[key]['normal'][0]),
                    normal=True)
            elif mhcii_preds[key]['predictor'] == 'Sturniolo':
                normal_df = _process_sturniolo_mhcii(
                    job.fileStore.readGlobalFile(mhcii_preds[key]['normal'][0]),
                    normal=True)
            elif mhcii_preds[key]['predictor'] == 'netMHCIIpan':
                normal_df = _process_net_mhcii(
                    job.fileStore.readGlobalFile(mhcii_preds[key]['normal'][0]),
                    normal=True)
            else:
                assert False
            normal_dict = normal_df.set_index('pept')['tumor_pred']
            normal_preds = [normal_dict[x] for x in list(tumor_df['normal_pept'])]
            tumor_df['normal_pred'] = normal_preds
            for pred in tumor_df.itertuples():
                print_mhc_peptide(pred, peptides, pepmap, mhcii_resfile,
                                  netmhc=mhcii_preds[key]['predictor'] == 'netMHCIIpan')
    if not(mhci_called or mhcii_called):
        raise RuntimeError('No peptides available for ranking')
    output_files = defaultdict()
    for mhc_file in [mhci_resfile.name, mhcii_resfile.name]:
        output_files[os.path.split(mhc_file)[1]] = job.fileStore.writeGlobalFile(mhc_file)
        export_results(job, output_files[os.path.split(mhc_file)[1]], mhc_file, univ_options,
                       subfolder='binding_predictions')

    return output_files


def print_mhc_peptide(neoepitope_info, peptides, pepmap, outfile, netmhc=False):
    """
    Accept data about one neoepitope from merge_mhc_peptide_calls and print it to outfile.  This is
    a generic module to reduce code redundancy.

    :param pandas.core.frame neoepitope_info: object containing with allele, pept, pred, core,
           normal_pept, normal_pred
    :param dict peptides: Dict of pepname: pep sequence for all IARS considered
    :param dict pepmap: Dict containing teh contents from the peptide map file.
    :param file outfile: An open file descriptor to the output file
    :param bool netmhc: Does this record correspond to a netmhcIIpan record? These are processed
           differently.
    """
    if netmhc:
        peptide_names = [neoepitope_info.peptide_name]
    else:
        peptide_names = [x for x, y in peptides.items() if neoepitope_info.pept in y]
    # Convert named tuple to dict so it can be modified
    neoepitope_info = neoepitope_info._asdict()
    # Handle fusion peptides (They are characterized by having all N's as the normal partner)
    if neoepitope_info['normal_pept'] == 'N' * len(neoepitope_info['pept']):
        neoepitope_info['normal_pept'] = neoepitope_info['normal_pred'] = 'NA'
    # For each peptide, append the ensembl gene
    for peptide_name in peptide_names:
        print('{ni[allele]}\t'
              '{ni[pept]}\t'
              '{ni[normal_pept]}\t'
              '{pname}\t'
              '{ni[core]}\t'
              '0\t'
              '{ni[tumor_pred]}\t'
              '{ni[normal_pred]}\t'
              '{pmap}'.format(ni=neoepitope_info, pname=peptide_name,
                                                  pmap=pepmap[peptide_name]), file=outfile)
    return None
