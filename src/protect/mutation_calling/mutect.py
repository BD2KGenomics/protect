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
from collections import defaultdict
from math import ceil
from protect.common import get_files_from_filestore, docker_path, docker_call, export_results, \
    untargz, gunzip
from protect.mutation_calling.common import sample_chromosomes, merge_perchrom_vcfs

import os
import sys


# disk for mutect.
from toil.job import PromisedRequirement


def mutect_disk(tumor_bam, normal_bam, fasta, dbsnp, cosmic):
    return int(ceil(tumor_bam.size) +
               ceil(normal_bam.size) +
               4 * ceil(fasta.size) +
               10 * ceil(dbsnp.size) +
               2 * ceil(cosmic.size))


def run_mutect_with_merge(job, tumor_bam, normal_bam, univ_options, mutect_options):
    """
    This is a convenience function that runs the entire mutect sub-graph.
    """
    spawn = job.wrapJobFn(run_mutect, tumor_bam, normal_bam, univ_options,
                          mutect_options).encapsulate()
    merge = job.wrapJobFn(merge_perchrom_vcfs, spawn.rv())
    job.addChild(spawn)
    spawn.addChild(merge)
    return merge.rv()


def run_mutect(job, tumor_bam, normal_bam, univ_options, mutect_options):
    """
    This module will spawn a mutect job for each chromosome on the DNA bams.

    ARGUMENTS
    1. tumor_bam: Dict of input tumor WGS/WSQ bam + bai
         tumor_bam
              |- 'tumor_fix_pg_sorted.bam': <JSid>
              +- 'tumor_fix_pg_sorted.bam.bai': <JSid>
    2. normal_bam: Dict of input normal WGS/WSQ bam + bai
         normal_bam
              |- 'normal_fix_pg_sorted.bam': <JSid>
              +- 'normal_fix_pg_sorted.bam.bai': <JSid>
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    4. mutect_options: Dict of parameters specific to mutect
         mutect_options
              |- 'dbsnp_vcf': <JSid for dnsnp vcf file>
              |- 'dbsnp_idx': <JSid for dnsnp vcf index file>
              |- 'cosmic_vcf': <JSid for cosmic vcf file>
              |- 'cosmic_idx': <JSid for cosmic vcf index file>
              |- 'genome_fasta': <JSid for genome fasta file>
              +- 'genome_dict': <JSid for genome fasta dict file>
              +- 'genome_fai': <JSid for genome fasta index file>

    RETURN VALUES
    1. perchrom_mutect: Dict of results of mutect per chromosome
         perchrom_mutect
              |- 'chr1'
              |   +- 'mutect_chr1.vcf': <JSid>
              |   +- 'mutect_chr1.out': <JSid>
              |- 'chr2'
              |   |- 'mutect_chr2.vcf': <JSid>
              |   +- 'mutect_chr2.out': <JSid>
             etc...

    This module corresponds to node 11 on the tree
    """
    # Get a list of chromosomes to handle
    chromosomes = sample_chromosomes(job, mutect_options['genome_fai'])
    perchrom_mutect = defaultdict()
    for chrom in chromosomes:
        perchrom_mutect[chrom] = job.addChildJobFn(
            run_mutect_perchrom, tumor_bam, normal_bam, univ_options, mutect_options, chrom,
            memory='6G', disk=PromisedRequirement(mutect_disk,
                                                  tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                                  normal_bam['normal_dna_fix_pg_sorted.bam'],
                                                  mutect_options['genome_fasta'],
                                                  mutect_options['dbsnp_vcf'],
                                                  mutect_options['cosmic_vcf'])).rv()
    return perchrom_mutect


def run_mutect_perchrom(job, tumor_bam, normal_bam, univ_options, mutect_options, chrom):
    """
    This module will run mutect on the DNA bams

    ARGUMENTS
    1. tumor_bam: REFER ARGUMENTS of spawn_mutect()
    2. normal_bam: REFER ARGUMENTS of spawn_mutect()
    3. univ_options: REFER ARGUMENTS of spawn_mutect()
    4. mutect_options: REFER ARGUMENTS of spawn_mutect()
    5. chrom: String containing chromosome name with chr appended

    RETURN VALUES
    1. output_files: Dict of results of mutect for chromosome
            output_files
              |- 'mutect_CHROM.vcf': <JSid>
              +- 'mutect_CHROM.out': <JSid>

    This module corresponds to node 12 on the tree
    """
    job.fileStore.logToMaster('Running mutect on %s:%s' % (univ_options['patient'], chrom))
    work_dir = os.getcwd()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'normal.bam': normal_bam['normal_dna_fix_pg_sorted.bam'],
        'normal.bam.bai': normal_bam['normal_dna_fix_pg_sorted.bam.bai'],
        'genome.fa.tar.gz': mutect_options['genome_fasta'],
        'genome.fa.fai.tar.gz': mutect_options['genome_fai'],
        'genome.dict.tar.gz': mutect_options['genome_dict'],
        'cosmic.vcf.tar.gz': mutect_options['cosmic_vcf'],
        'cosmic.vcf.idx.tar.gz': mutect_options['cosmic_idx'],
        'dbsnp.vcf.gz': mutect_options['dbsnp_vcf'],
        'dbsnp.vcf.idx.tar.gz': mutect_options['dbsnp_idx']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    # dbsnp.vcf should be bgzipped, but all others should be tar.gz'd
    input_files['dbsnp.vcf'] = gunzip(input_files['dbsnp.vcf.gz'])
    for key in ('genome.fa', 'genome.fa.fai', 'genome.dict', 'cosmic.vcf', 'cosmic.vcf.idx',
                'dbsnp.vcf.idx'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    mutout = ''.join([work_dir, '/', chrom, '.out'])
    mutvcf = ''.join([work_dir, '/', chrom, '.vcf'])
    parameters = ['-R', input_files['genome.fa'],
                  '--cosmic', input_files['cosmic.vcf'],
                  '--dbsnp', input_files['dbsnp.vcf'],
                  '--input_file:normal', input_files['normal.bam'],
                  '--input_file:tumor', input_files['tumor.bam'],
                  # '--tumor_lod', str(10),
                  # '--initial_tumor_lod', str(4.0),
                  '-L', chrom,
                  '--out', docker_path(mutout),
                  '--vcf', docker_path(mutvcf)
                  ]
    java_xmx = mutect_options['java_Xmx'] if mutect_options['java_Xmx'] \
        else univ_options['java_Xmx']
    docker_call(tool='mutect', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=java_xmx,
                tool_version=mutect_options['version'])
    output_file = job.fileStore.writeGlobalFile(mutvcf)
    export_results(job, output_file, mutvcf, univ_options, subfolder='mutations/mutect')
    return output_file


def process_mutect_vcf(job, mutect_vcf, work_dir, univ_options):
    """
    This will process all the mutect vcfs to return only passing calls
    :param job: job
    :param str mutect_vcf: Job Store ID corresponding to a mutect vcf for 1 chromosome
    :param univ_options: Universal options
    :returns dict: Dict with chromosomes as keys and path to the corresponding parsed mutect vcfs as
                   values
    """
    mutect_vcf = job.fileStore.readGlobalFile(mutect_vcf)

    with open(mutect_vcf, 'r') as infile, open(mutect_vcf + 'mutect_parsed.tmp', 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith('#'):
                print(line, file=outfile)
                continue
            line = line.split('\t')
            if line[6] != 'REJECT':
                print('\t'.join(line), file=outfile)
    return outfile.name
