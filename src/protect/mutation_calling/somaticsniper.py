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
from protect.common import (get_files_from_filestore,
                            docker_path,
                            docker_call,
                            untargz,
                            export_results)
from protect.mutation_calling.common import (chromosomes_from_fai,
                                             sample_chromosomes)
from toil.job import PromisedRequirement

import os


# disk for somatic sniper, and for filtering
def sniper_disk(tumor_bam, normal_bam, fasta):
    return int(ceil(tumor_bam.size) +
               ceil(normal_bam.size) +
               4 * ceil(fasta.size))


def pileup_disk(tumor_bam, fasta):
    return int(ceil(tumor_bam.size) +
               2 * ceil(fasta.size))


def sniper_filter_disk(tumor_bam, fasta):
    return int(ceil(tumor_bam.size) +
               2 * ceil(fasta.size))


# Somatic Sniper is a different tool from the other callers because it is run in full
def run_somaticsniper_with_merge(job, tumor_bam, normal_bam, univ_options, somaticsniper_options):
    """
    This is a convenience function that runs the entire somaticsniper sub-graph.
    """
    spawn = job.wrapJobFn(run_somaticsniper, tumor_bam, normal_bam, univ_options,
                          somaticsniper_options, split=False).encapsulate()
    job.addChild(spawn)
    return spawn.rv()


def run_somaticsniper(job, tumor_bam, normal_bam, univ_options, somaticsniper_options, split=True):
    """
    This module will spawn a somaticsniper job for each chromosome on the DNA bams.

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
    4. somaticsniper_options: Dict of parameters specific to somaticsniper
         somaticsniper_options
              |- 'dbsnp_vcf': <JSid for dnsnp vcf file>
              |- 'dbsnp_idx': <JSid for dnsnp vcf index file>
              |- 'cosmic_vcf': <JSid for cosmic vcf file>
              |- 'cosmic_idx': <JSid for cosmic vcf index file>
              |- 'genome_fasta': <JSid for genome fasta file>
              +- 'genome_dict': <JSid for genome fasta dict file>
              +- 'genome_fai': <JSid for genome fasta index file>

    RETURN VALUES
    1. perchrom_somaticsniper: Dict of results of somaticsniper per chromosome
         perchrom_somaticsniper
              |- 'chr1'
              |   +- <JSid for somaticsniper_chr1.vcf>
              |- 'chr2'
              |   +- <JSid for somaticsniper_chr2.vcf>
             etc...

    This module corresponds to node 11 on the tree
    """
    # Get a list of chromosomes to handle
    chromosomes = sample_chromosomes(job, somaticsniper_options['genome_fai'])
    perchrom_somaticsniper = defaultdict()
    snipe = job.wrapJobFn(run_somaticsniper_full, tumor_bam, normal_bam, univ_options,
                          somaticsniper_options,
                          disk=PromisedRequirement(sniper_disk,
                                                   tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                                   normal_bam['normal_dna_fix_pg_sorted.bam'],
                                                   somaticsniper_options['genome_fasta']),
                          memory='6G')
    pileup = job.wrapJobFn(run_pileup, tumor_bam, univ_options, somaticsniper_options,
                           disk=PromisedRequirement(pileup_disk,
                                                    tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                                    somaticsniper_options['genome_fasta']),
                           memory='6G')
    filtersnipes = job.wrapJobFn(filter_somaticsniper, tumor_bam, snipe.rv(), pileup.rv(),
                                 univ_options, somaticsniper_options,
                                 disk=PromisedRequirement(sniper_filter_disk,
                                                          tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                                          somaticsniper_options['genome_fasta']),
                                 memory='6G')

    job.addChild(snipe)
    job.addChild(pileup)
    snipe.addChild(filtersnipes)
    pileup.addChild(filtersnipes)
    if split:
        unmerge_snipes = job.wrapJobFn(unmerge, filtersnipes.rv(), somaticsniper_options,
                                       univ_options)
        filtersnipes.addChild(unmerge_snipes)
        return unmerge_snipes.rv()
    else:
        return filtersnipes.rv()


def run_somaticsniper_full(job, tumor_bam, normal_bam, univ_options, somaticsniper_options):
    """
    This module will run somaticsniper on the DNA bams.

    ARGUMENTS
    :param dict tumor_bam: REFER ARGUMENTS of spawn_somaticsniper()
    :param dict normal_bam: REFER ARGUMENTS of spawn_somaticsniper()
    :param dict univ_options: REFER ARGUMENTS of spawn_somaticsniper()
    :param dict somaticsniper_options: REFER ARGUMENTS of spawn_somaticsniper()

    RETURN VALUES
    :returns: dict of output vcfs for each chromosome
    :rtype: dict
    """
    job.fileStore.logToMaster('Running somaticsniper on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'normal.bam': normal_bam['normal_dna_fix_pg_sorted.bam'],
        'normal.bam.bai': normal_bam['normal_dna_fix_pg_sorted.bam.bai'],
        'genome.fa.tar.gz': somaticsniper_options['genome_fasta'],
        'genome.fa.fai.tar.gz': somaticsniper_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    output_file = os.path.join(work_dir, 'somatic-sniper_full.vcf')
    parameters = ['-f', input_files['genome.fa'],
                  '-F', 'vcf',
                  '-G',
                  '-L',
                  '-q', '1',
                  '-Q', '15',
                  input_files['tumor.bam'],
                  input_files['normal.bam'],
                  docker_path(output_file)]
    docker_call(tool='somaticsniper', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    outfile = job.fileStore.writeGlobalFile(output_file)
    return outfile


def filter_somaticsniper(job, tumor_bam, somaticsniper_output, tumor_pileup, univ_options,
                         somaticsniper_options):
    """
    This module will filter the somaticsniper output for a single chromosome

    :param toil.Job job: Job
    :param dict tumor_bam: Tumor bam file and it's bai
    :param str somaticsniper_output: jsID from somatic sniper
    :param str tumor_pileup: jsID for pileup file for this chromsome
    :param dict univ_options: Universal options
    :param dict somaticsniper_options: Options specific to Somatic Sniper
    :returns: filtered chromsome vcf
    :rtype: str
    """
    job.fileStore.logToMaster('Filtering somaticsniper for %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'input.vcf': somaticsniper_output,
        'pileup.txt': tumor_pileup,
        'genome.fa.tar.gz': somaticsniper_options['genome_fasta'],
        'genome.fa.fai.tar.gz': somaticsniper_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    # Run snpfilter.pl
    parameters = ['snpfilter.pl',
                  '--snp-file', input_files['input.vcf'],
                  '--indel-file', input_files['pileup.txt']]
    # Creates /data/input.vcf.SNPfilter
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])

    # Run prepare_for_readcount.pl
    parameters = ['prepare_for_readcount.pl',
                  '--snp-file', input_files['input.vcf'] + '.SNPfilter']
    # Creates /data/input.vcf.SNPfilter.pos
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])

    # Run  bam-readcount
    parameters = ['-b', '15',
                  '-f', input_files['genome.fa'],
                  '-l', input_files['input.vcf'] + '.SNPfilter.pos',
                  '-w', '1',
                  input_files['tumor.bam']]
    # Creates the read counts file
    with open(os.path.join(work_dir, 'readcounts.txt'), 'w') as readcounts_file:
        docker_call(tool='bam-readcount', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=readcounts_file)

    # Run fpfilter.pl
    parameters = ['fpfilter.pl',
                  '--snp-file', input_files['input.vcf'] + '.SNPfilter',
                  '--readcount-file', docker_path(readcounts_file.name)]

    # Creates input.vcf.SNPfilter.fp_pass and input.vcf.SNPfilter.fp_fail
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])

    # Run highconfidence.pl
    parameters = ['highconfidence.pl',
                  '--snp-file', input_files['input.vcf'] + '.SNPfilter.fp_pass']

    # Creates input.vcf.SNPfilter.fp_pass.hc
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])

    outfile = job.fileStore.writeGlobalFile(os.path.join(os.getcwd(),
                                                         'input.vcf.SNPfilter.fp_pass.hc'))
    return outfile


def process_somaticsniper_vcf(job, somaticsniper_vcf, work_dir, univ_options):
    """
    This does nothing for now except to download and return the somaticsniper vcfs
    :param job: job
    :param str somaticsniper_vcf: Job Store ID corresponding to a somaticsniper vcf for 1 chromosome
    :param univ_options: Universal options
    :returns dict: Dict with chromosomes as keys and path to the corresponding somaticsniper vcfs as
    values
    """
    somaticsniper_vcf = job.fileStore.readGlobalFile(somaticsniper_vcf)

    with open(somaticsniper_vcf, 'r') as infile, open(somaticsniper_vcf +
                                                      'somaticsniper_parsed.tmp', 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith('#'):
                print(line, file=outfile)
                continue
            print(line, file=outfile)
    return outfile.name


def unmerge(job, input_vcf, somaticsniper_options, univ_options):
    """
    Un-merges a vcf file into a file per chromosome.

    :param str input_vcf: Input vcf
    :param dict somaticsniper_options: Options specific to Somatic Sniper
    :param dict univ_options: Universal options
    :returns: dict of jsIDs, onr for each chromosomal vcf
    :rtype: dict
    """
    work_dir = os.getcwd()
    input_files = {
        'input.vcf': input_vcf,
        'genome.fa.fai.tar.gz': somaticsniper_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    input_files['genome.fa.fai'] = untargz(input_files['genome.fa.fai.tar.gz'], work_dir)

    chromosomes = chromosomes_from_fai(input_files['genome.fa.fai'])

    read_chromosomes = defaultdict()
    with open(input_files['input.vcf'], 'r') as in_vcf:
        header = []
        for line in in_vcf:
            if line.startswith('#'):
                header.append(line)
                continue
            line = line.strip()
            chrom = line.split()[0]
            if chrom in read_chromosomes:
                print(line, file=read_chromosomes[chrom])
            else:
                read_chromosomes[chrom] = open(os.path.join(os.getcwd(), chrom + '.vcf'), 'w')
                print(''.join(header), file=read_chromosomes[chrom], end='')
                print(line, file=read_chromosomes[chrom])
    # Process chromosomes that had no mutations
    for chrom in set(chromosomes).difference(set(read_chromosomes.keys())):
        read_chromosomes[chrom] = open(os.path.join(os.getcwd(), chrom + '.vcf'), 'w')
        print(''.join(header), file=read_chromosomes[chrom], end='')
    outdict = {}
    for chrom, chromvcf in read_chromosomes.items():
        chromvcf.close()
        export_results(job, chromvcf.name, univ_options, subfolder='mutations/somaticsniper')
        outdict[chrom] = job.fileStore.writeGlobalFile(chromvcf.name)
    return outdict


def run_pileup(job, tumor_bam, univ_options, somaticsniper_options):
    """
    Runs a samtools pileup on the tumor bam.

    :param toil.Job job: job
    :param dict tumor_bam: Tumor bam file
    :param dict univ_options: Universal Options
    :returns: jsID for the chromsome pileup file
    :rtype: str
    """
    job.fileStore.logToMaster(
        'Running samtools pileup on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'genome.fa.tar.gz': somaticsniper_options['genome_fasta'],
        'genome.fa.fai.tar.gz': somaticsniper_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = ['pileup',
                  '-cvi',
                  '-f', docker_path(input_files['genome.fa']),
                  docker_path(input_files['tumor.bam'])]

    with open(os.path.join(work_dir, 'pileup.txt'), 'w') as pileup_file:
        docker_call(tool='samtools:0.1.8', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=pileup_file)
    outfile = job.fileStore.writeGlobalFile(pileup_file.name)
    return outfile
