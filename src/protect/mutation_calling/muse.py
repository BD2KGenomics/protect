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

import shutil
from collections import defaultdict

import sys

import time

from protect.common import get_files_from_filestore, docker_path, docker_call, untargz, \
    export_results
from protect.mutation_calling.common import sample_chromosomes, merge_perchrom_vcfs

import os


def run_muse_with_merge(job, tumor_bam, normal_bam, univ_options, muse_options):
    """
    This is a convenience function that runs the entire muse sub-graph.
    """
    spawn = job.wrapJobFn(run_muse, tumor_bam, normal_bam, univ_options,
                          muse_options).encapsulate()
    merge = job.wrapJobFn(merge_perchrom_vcfs, spawn.rv())
    spawn.addChild(merge)
    return merge.rv()


def run_muse(job, tumor_bam, normal_bam, univ_options, muse_options):
    """
    This module will spawn a muse job for each chromosome on the DNA bams.

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
    4. muse_options: Dict of parameters specific to muse
         muse_options
              |- 'dbsnp_vcf': <JSid for dnsnp vcf file>
              |- 'dbsnp_idx': <JSid for dnsnp vcf index file>
              |- 'cosmic_vcf': <JSid for cosmic vcf file>
              |- 'cosmic_idx': <JSid for cosmic vcf index file>
              |- 'genome_fasta': <JSid for genome fasta file>
              +- 'genome_dict': <JSid for genome fasta dict file>
              +- 'genome_fai': <JSid for genome fasta index file>

    RETURN VALUES
    1. perchrom_muse: Dict of results of muse per chromosome
         perchrom_muse
              |- 'chr1'
              |   +- <JSid for muse_chr1.vcf>
              |- 'chr2'
              |   +- <JSid for muse_chr2.vcf>
             etc...

    This module corresponds to node 11 on the tree
    """
    # Get a list of chromosomes to handle
    chromosomes = sample_chromosomes(job, muse_options['genome_fai'])
    perchrom_muse = defaultdict()
    for chrom in chromosomes:
        call = job.addChildJobFn(run_muse_perchrom, tumor_bam, normal_bam, univ_options,
                                 muse_options, chrom, disk='60G', memory='6G')
        sump = call.addChildJobFn(run_muse_sump_perchrom, call.rv(), univ_options, muse_options,
                                  chrom, disk='60G', memory='6G')
        perchrom_muse[chrom] = sump.rv()
    return perchrom_muse


def run_muse_perchrom(job, tumor_bam, normal_bam, univ_options, muse_options, chrom):
    """
    This module will run muse on the DNA bams

    ARGUMENTS
    1. tumor_bam: REFER ARGUMENTS of spawn_muse()
    2. normal_bam: REFER ARGUMENTS of spawn_muse()
    3. univ_options: REFER ARGUMENTS of spawn_muse()
    4. muse_options: REFER ARGUMENTS of spawn_muse()
    5. chrom: String containing chromosome name with chr appended

    RETURN VALUES
    1. output_files: <JSid for CHROM.MuSe.txt>

    This module corresponds to node 12 on the tree
    """
    job.fileStore.logToMaster('Running muse on %s:%s' % (univ_options['patient'], chrom))
    work_dir = os.getcwd()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'normal.bam': normal_bam['normal_dna_fix_pg_sorted.bam'],
        'normal.bam.bai': normal_bam['normal_dna_fix_pg_sorted.bam.bai'],
        'genome.fa.tar.gz': muse_options['genome_fasta'],
        'genome.fa.fai.tar.gz': muse_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    output_prefix = os.path.join(work_dir, chrom)

    parameters = ['call',
                  '-f', input_files['genome.fa'],
                  '-r', chrom,
                  '-O', docker_path(output_prefix),
                  input_files['tumor.bam'],
                  input_files['normal.bam']]
    docker_call(tool='muse', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    outfile = job.fileStore.writeGlobalFile(''.join([output_prefix, '.MuSE.txt']))
    return outfile


def run_muse_sump_perchrom(job, muse_output, univ_options, muse_options, chrom):
    """
    This module will run muse sump on the muse output
    """
    job.fileStore.logToMaster('Running muse sump on %s:%s' % (univ_options['patient'], chrom))
    work_dir = os.getcwd()
    input_files = {
        'MuSE.txt': muse_output,
        'dbsnp_coding.vcf.gz': muse_options['dbsnp_vcf'],
        'dbsnp_coding.vcf.gz.tbi.tmp': muse_options['dbsnp_tbi']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    tbi = os.path.splitext(input_files['dbsnp_coding.vcf.gz.tbi.tmp'])[0]
    print({x: os.stat(x) for x in os.listdir(work_dir)}, file=sys.stderr)
    time.sleep(2)
    shutil.copy(input_files['dbsnp_coding.vcf.gz.tbi.tmp'], tbi)
    os.chmod(tbi, 0777)
    open(tbi, 'a').close()
    input_files = {key: docker_path(path) for key, path in input_files.items()}
    print({x: os.stat(x) for x in os.listdir(work_dir)}, file=sys.stderr)
    output_file = ''.join([work_dir, '/', chrom, '.vcf'])

    parameters = ['sump',
                  '-I', input_files['MuSE.txt'],
                  '-O', docker_path(output_file),
                  '-D', input_files['dbsnp_coding.vcf.gz'],
                  '-E']

    docker_call(tool='muse', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    export_results(job, output_file, univ_options, subfolder='mutations/muse')
    outfile = job.fileStore.writeGlobalFile(output_file)
    return outfile


def process_muse_vcf(job, muse_vcf, work_dir, univ_options):
    """
    This does nothing for now except to download and return the muse vcfs
    :param job: job
    :param str muse_vcf: Job Store ID corresponding to a muse vcf for 1 chromosome
    :param univ_options: Universal options
    :returns dict: Dict with chromosomes as keys and path to the corresponding muse vcfs as values
    """
    muse_vcf = job.fileStore.readGlobalFile(muse_vcf)

    with open(muse_vcf, 'r') as infile, open(muse_vcf + 'muse_parsed.tmp', 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith('#'):
                print(line, file=outfile)
                continue
            line = line.split('\t')
            if line[6] == 'PASS':
                print('\t'.join(line), file=outfile)
    return outfile.name
