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
from protect.addons.common import TCGAToGTEx
from protect.common import export_results, get_files_from_filestore, untargz
from protect.haplotyping.phlat import parse_phlat_file

import os
import pandas as pd


def run_car_t_validity_assessment(job, rsem_files, univ_options, reports_options):
    """
    A wrapper for assess_car_t_validity.

    :param dict rsem_files: Results from running rsem
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict reports_options: Options specific to reporting modules
    :return: The results of running assess_car_t_validity
    :rtype: toil.fileStore.FileID
    """
    return job.addChildJobFn(assess_car_t_validity, rsem_files['rsem.genes.results'],
                             univ_options, reports_options).rv()


def assess_car_t_validity(job, gene_expression, univ_options, reports_options):
    """
    This function creates a report on the available clinical trials and scientific literature
    available for the overexpressed genes in the specified tumor type.
    It also gives a list of clinical trials available for other types of cancer with the same
    overexpressed gene.

    :param toil.fileStore.FileID gene_expression: The resm gene expression
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict reports_options: Options specific to reporting modules
    :return: The results of running assess_car_t_validity
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()

    tumor_type = univ_options['tumor_type']

    input_files = {
        'rsem_quant.tsv': gene_expression,
        'car_t_targets.tsv.tar.gz': reports_options['car_t_targets_file']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    input_files['car_t_targets.tsv'] = untargz(input_files['car_t_targets.tsv.tar.gz'],
                                                 work_dir)

    target_data = pd.read_table(input_files['car_t_targets.tsv'], index_col=0)
    patient_df = pd.read_csv('rsem_quant.tsv', sep=' ', delimiter='\t', header='infer', index_col=0)
    patient_df.index = (patient_df.index).str.replace('\\..*$', '')

    overexpressed = []
    # Check if the tumor has a corresponding normal
    try:
        tissue_of_origin = TCGAToGTEx[tumor_type]
    except KeyError:
        tissue_of_origin = 'NA'
    # Write the report
    with open('car_t_target_report.txt', 'w') as car_t_report:
        #print(target_data.index, file=car_t_report)
        if tissue_of_origin in target_data.index:
            print('Available clinical trials for ' + str.lower(tissue_of_origin) +
                  ' cancer with GTEX and TCGA median values', file=car_t_report)

            print(('\t{:10}{:<10}{:<10}{:<10}{:<40}{:<12}\n'.format('Gene', 'GTEX',
                                                                    'TCGA N', 'Observed',
                                                                    'DOI for gene papers',
                                                                    'Clinical Trials')),
                                                                     file=car_t_report)
            collected_values = []
            # Get the gene name, GTEX, TCGA, and observed values
            for index, row in target_data.iterrows():
                if index == tissue_of_origin:
                    gene = row['ENSG']
                    gtex = '{0:.2f}'.format(float(row['GTEX']))
                    tcga = '{0:.2f}'.format(float(row['TCGA']))
                    observed = '{0:.2f}'.format(
                        float(patient_df.loc[gene, 'TPM'])) if gene in patient_df.index else 'NA'
                    doi = row['DOI']
                    target = str.upper(row['TARGET'])
                    clinical_trial = row['Clinical trials']
                    collection = [target, gtex, tcga, observed, doi, clinical_trial]
                    collected_values.append(collection)
                    if observed != 'NA':
                        if float(gtex) <= float(observed) or float(tcga) <= float(observed):
                            overexpressed.append(gene)

            collected_values = sorted(collected_values, key=lambda col: float(col[3]), reverse=True)
            for entry in collected_values:
                print(('\t{:10}{:<10}{:<10}{:<10}{:<40}{:<12}'.format(entry[0],
                                                                      entry[1], entry[2],
                                                                      str(entry[3]), entry[4],
                                                                      entry[5])), file=car_t_report)

            print('\nBased on the genes overexpressed in this cancer type, here\'s a list of clinical '
                 'trials for other types of cancer', file=car_t_report)
            if len(overexpressed) != 0:
                # Check if there are other clinical trials for other cancer types
                print(('\t{:10}{:<10}{:<10}{:<10}{:<40}{:<17}{:<20}\n'.format('Gene', 'GTEX',
                                                                              'TCGA N', 'Observed',
                                                                              'DOI for gene papers',
                                                                              'Clinical Trials',
                                                                              'Cancer')),
                                                                              file=car_t_report)
                other_trials = []
                for index, row in target_data.iterrows():
                    if row['ENSG'] in overexpressed and index != tissue_of_origin:
                        gene = row['ENSG']
                        gtex = '{0:.2f}'.format(float(row['GTEX']))
                        tcga = '{0:.2f}'.format(float(row['TCGA']))
                        doi = row['DOI']
                        target = str.upper(row['TARGET'])
                        observed = '{0:.2f}'.format(
                            float(patient_df.loc[gene, 'TPM'])) if gene in patient_df.index else 'NA'
                        collected_values = [target, gtex, tcga, observed, doi, row['Clinical trials'],
                                            index]
                        other_trials.append(collected_values)

                other_trials = sorted(other_trials, key=lambda col: col[0])
                for entry in other_trials:
                    print(('\t{:10}{:<10}{:<10}{:<10}{:<40}{:<17}{:<20}'.format(entry[0], entry[1],
                                                                                entry[2], entry[3],
                                                                                entry[4], entry[5],
                                                                                entry[6])),
                                                                                file=car_t_report)
            else:
                print("Data not available", file=car_t_report)

        else:
            print('Data not available for ' + tumor_type, file=car_t_report)

    output_file = job.fileStore.writeGlobalFile(car_t_report.name)
    export_results(job, output_file, car_t_report.name, univ_options, subfolder='reports')
    job.fileStore.logToMaster('Ran car t validity assessment on %s successfully'
                              % univ_options['patient'])
    return output_file
