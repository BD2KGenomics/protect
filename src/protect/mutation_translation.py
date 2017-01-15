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

from protect.common import (docker_call,
                            get_files_from_filestore,
                            export_results,
                            untargz,
                            docker_path)

import os


def transgene_disk(rna_bamfiles, tdna_bam=None):
    return int(ceil(rna_bamfiles['rna_genome']['rna_genome_sorted.bam'].size) +
               (ceil(tdna_bam['tumor_dna_fix_pg_sorted.bam'].size) if tdna_bam is not None else 0) +
               104857600)


def run_transgene(job, snpeffed_file, rna_bam, univ_options, transgene_options, tumor_dna_bam=None, fusion_calls=None):
    """
    Run transgene on an input snpeffed vcf file and return the peptides for MHC prediction.


    :param toil.fileStore.FileID snpeffed_file: fsID for snpeffed vcf
    :param dict rna_bam: The dict of bams returned by running star
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict transgene_options: Options specific to Transgene
    :param dict tumor_dna_bam: The dict of bams returned by running bwa
    :return: A dictionary of 9 files (9-, 10-, and 15-mer peptides each for Tumor and Normal and the
             corresponding .map files for the 3 Tumor fastas)
             output_files:
                 |- 'transgened_normal_10_mer_snpeffed.faa': fsID
                 |- 'transgened_normal_15_mer_snpeffed.faa': fsID
                 |- 'transgened_normal_9_mer_snpeffed.faa': fsID
                 |- 'transgened_tumor_10_mer_snpeffed.faa': fsID
                 |- 'transgened_tumor_10_mer_snpeffed.faa.map': fsID
                 |- 'transgened_tumor_15_mer_snpeffed.faa': fsID
                 |- 'transgened_tumor_15_mer_snpeffed.faa.map': fsID
                 |- 'transgened_tumor_9_mer_snpeffed.faa': fsID
                 +- 'transgened_tumor_9_mer_snpeffed.faa.map': fsID
    :rtype: dict
    """
    job.fileStore.logToMaster('Running transgene on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'snpeffed_muts.vcf': snpeffed_file,
        'rna.bam': rna_bam['rna_genome']['rna_genome_sorted.bam'],
        'rna.bam.bai': rna_bam['rna_genome']['rna_genome_sorted.bam.bai'],
        'pepts.fa.tar.gz': transgene_options['gencode_peptide_fasta']}
    if tumor_dna_bam is not None:
        input_files.update({
            'tumor_dna.bam': tumor_dna_bam['tumor_dna_fix_pg_sorted.bam'],
            'tumor_dna.bam.bai': tumor_dna_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        })
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    input_files['pepts.fa'] = untargz(input_files['pepts.fa.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['--peptides', input_files['pepts.fa'],
                  '--snpeff', input_files['snpeffed_muts.vcf'],
                  '--rna_file', input_files['rna.bam'],
                  '--prefix', 'transgened',
                  '--pep_lens', '9,10,15',
                  '--cores', str(transgene_options['n'])]

    if tumor_dna_bam is not None:
        parameters.extend(['--dna_file', input_files['tumor_dna.bam']])

    if fusion_calls:
        fusion_files = {'fusion_calls': fusion_calls,
                        'transcripts.fa.tar.gz': transgene_options['gencode_transcript_fasta'],
                        'annotation.gtf.tar.gz': transgene_options['gencode_annotation_gtf'],
                        'genome.fa.tar.gz': transgene_options['genome_fasta']}

        fusion_files = get_files_from_filestore(job, fusion_files, work_dir, docker=False)
        fusion_files['transcripts.fa'] = untargz(fusion_files['transcripts.fa.tar.gz'], work_dir)
        fusion_files['genome.fa'] = untargz(fusion_files['genome.fa.tar.gz'], work_dir)
        fusion_files['annotation.gtf'] = untargz(fusion_files['annotation.gtf.tar.gz'], work_dir)
        fusion_files = {key: docker_path(path) for key, path in fusion_files.items()}
        parameters += ['--transcripts', fusion_files['transcripts.fa'],
                       '--fusions', fusion_files['fusion_calls'],
                       '--genome', fusion_files['genome.fa'],
                       '--annotation', fusion_files['annotation.gtf']]

    docker_call(tool='transgene',
                tool_parameters=parameters,
                work_dir=work_dir,
                dockerhub=univ_options['dockerhub'],
                tool_version=transgene_options['version'])

    output_files = defaultdict()
    for peplen in ['9', '10', '15']:
        for tissue_type in ['tumor', 'normal']:
            pepfile = '_'.join(['transgened', tissue_type,  peplen, 'mer_snpeffed.faa'])
            output_files[pepfile] = job.fileStore.writeGlobalFile(os.path.join(work_dir, pepfile))
            export_results(job, output_files[pepfile], pepfile, univ_options, subfolder='peptides')
        mapfile = '_'.join(['transgened_tumor', peplen, 'mer_snpeffed.faa.map'])
        output_files[mapfile] = job.fileStore.writeGlobalFile(os.path.join(work_dir, mapfile))
        export_results(job, output_files[mapfile], mapfile, univ_options, subfolder='peptides')
    os.rename('transgened_transgened.vcf', 'mutations.vcf')
    export_results(job, job.fileStore.writeGlobalFile('mutations.vcf'), 'mutations.vcf',
                   univ_options, subfolder='mutations/transgened')
    return output_files
