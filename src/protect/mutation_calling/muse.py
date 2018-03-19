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

from protect.common import (docker_call,
                            docker_path,
                            export_results,
                            get_files_from_filestore,
                            untargz)
from protect.mutation_calling.common import merge_perchrom_vcfs, sample_chromosomes
from toil.job import PromisedRequirement

import os
import shutil
import time


# disk for muse and muse_sump.
def muse_disk(tumor_bam, normal_bam, fasta):
    return int(ceil(tumor_bam.size) +
               ceil(normal_bam.size) +
               5 * ceil(fasta.size))


def muse_sump_disk(dbsnp):
    return int(1.5 * ceil(dbsnp.size + 524288))


def run_muse_with_merge(job, tumor_bam, normal_bam, univ_options, muse_options):
    """
    A wrapper for the the entire MuSE sub-graph.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict muse_options: Options specific to MuSE
    :return: fsID to the merged MuSE calls
    :rtype: toil.fileStore.FileID
    """
    spawn = job.wrapJobFn(run_muse, tumor_bam, normal_bam, univ_options, muse_options,
                          disk='100M').encapsulate()
    merge = job.wrapJobFn(merge_perchrom_vcfs, spawn.rv(), disk='100M')
    job.addChild(spawn)
    spawn.addChild(merge)
    return merge.rv()


def run_muse(job, tumor_bam, normal_bam, univ_options, muse_options):
    """
    Spawn a MuSE job for each chromosome on the DNA bams.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict muse_options: Options specific to MuSE
    :return: Dict of results from running MuSE on every chromosome
             perchrom_muse:
                 |- 'chr1': fsID
                 |- 'chr2' fsID
                 |
                 |-...
                 |
                 +- 'chrM': fsID
    :rtype: dict
    """
    # Get a list of chromosomes to handle
    if muse_options['chromosomes']:
        chromosomes = muse_options['chromosomes']
    else:
        chromosomes = sample_chromosomes(job, muse_options['genome_fai'])
    perchrom_muse = defaultdict()
    for chrom in chromosomes:
        call = job.addChildJobFn(run_muse_perchrom, tumor_bam, normal_bam, univ_options,
                                 muse_options, chrom, disk=PromisedRequirement(
                                     muse_disk,
                                     tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                     normal_bam['normal_dna_fix_pg_sorted.bam'],
                                     muse_options['genome_fasta']),
                                 memory='6G')
        sump = call.addChildJobFn(run_muse_sump_perchrom, call.rv(), univ_options, muse_options,
                                  chrom,
                                  disk=PromisedRequirement(muse_sump_disk,
                                                           muse_options['dbsnp_vcf']),
                                  memory='6G')
        perchrom_muse[chrom] = sump.rv()
    return perchrom_muse


def run_muse_perchrom(job, tumor_bam, normal_bam, univ_options, muse_options, chrom):
    """
    Run MuSE call on a single chromosome in the input bams.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict muse_options: Options specific to MuSE
    :param str chrom: Chromosome to process
    :return: fsID for the chromsome vcf
    :rtype: toil.fileStore.FileID
    """
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
                dockerhub=univ_options['dockerhub'], tool_version=muse_options['version'])
    outfile = job.fileStore.writeGlobalFile(''.join([output_prefix, '.MuSE.txt']))
    job.fileStore.logToMaster('Ran MuSE on %s:%s successfully' % (univ_options['patient'], chrom))
    return outfile


def run_muse_sump_perchrom(job, muse_output, univ_options, muse_options, chrom):
    """
    Run MuSE sump on the MuSE call generated vcf.

    :param toil.fileStore.FileID muse_output: vcf generated by MuSE call
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict muse_options: Options specific to MuSE
    :param str chrom: Chromosome to process
    :return: fsID for the chromsome vcf
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    input_files = {
        'MuSE.txt': muse_output,
        'dbsnp_coding.vcf.gz': muse_options['dbsnp_vcf'],
        'dbsnp_coding.vcf.gz.tbi.tmp': muse_options['dbsnp_tbi']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    tbi = os.path.splitext(input_files['dbsnp_coding.vcf.gz.tbi.tmp'])[0]
    time.sleep(2)
    shutil.copy(input_files['dbsnp_coding.vcf.gz.tbi.tmp'], tbi)
    os.chmod(tbi, 0777)
    open(tbi, 'a').close()
    input_files = {key: docker_path(path) for key, path in input_files.items()}
    output_file = ''.join([work_dir, '/', chrom, '.vcf'])

    parameters = ['sump',
                  '-I', input_files['MuSE.txt'],
                  '-O', docker_path(output_file),
                  '-D', input_files['dbsnp_coding.vcf.gz'],
                  '-E']

    docker_call(tool='muse', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=muse_options['version'])
    outfile = job.fileStore.writeGlobalFile(output_file)
    export_results(job, outfile, output_file, univ_options, subfolder='mutations/muse')
    job.fileStore.logToMaster('Ran MuSE sump on %s:%s successfully'
                              % (univ_options['patient'], chrom))
    return outfile


def process_muse_vcf(job, muse_vcf, work_dir, univ_options):
    """
    Process the MuSE vcf for accepted calls.

    :param toil.fileStore.FileID muse_vcf: fsID for a MuSE generated chromosome vcf
    :param str work_dir: Working directory
    :param dict univ_options: Dict of universal options used by almost all tools
    :return: Path to the processed vcf
    :rtype: str
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
