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

from protect.common import (docker_path,
                            docker_call,
                            get_files_from_filestore,
                            untargz)
from protect.mutation_calling.common import (sample_chromosomes,
                                             unmerge)
from toil.job import PromisedRequirement

import os


# disk for somatic sniper, and for filtering
def sniper_disk(tumor_bam, normal_bam, fasta):
    return int(ceil(tumor_bam.size) +
               ceil(normal_bam.size) +
               5 * ceil(fasta.size))


def pileup_disk(tumor_bam, fasta):
    return int(ceil(tumor_bam.size) +
               5 * ceil(fasta.size))


def sniper_filter_disk(tumor_bam, fasta):
    return int(ceil(tumor_bam.size) +
               5 * ceil(fasta.size))


# Somatic Sniper is a different tool from the other callers because it is run in full.  Certain
# redundant modules in this file are kept to maintain a similar naming schema to other callers.
def run_somaticsniper_with_merge(job, tumor_bam, normal_bam, univ_options, somaticsniper_options):
    """
    A wrapper for the the entire SomaticSniper sub-graph.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict somaticsniper_options: Options specific to SomaticSniper
    :return: fsID to the merged SomaticSniper calls
    :rtype: toil.fileStore.FileID
    """
    spawn = job.wrapJobFn(run_somaticsniper, tumor_bam, normal_bam, univ_options,
                          somaticsniper_options, split=False).encapsulate()
    job.addChild(spawn)
    return spawn.rv()


def run_somaticsniper(job, tumor_bam, normal_bam, univ_options, somaticsniper_options, split=True):
    """
    Run the SomaticSniper subgraph on the DNA bams.  Optionally split the results into
    per-chromosome vcfs.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict somaticsniper_options: Options specific to SomaticSniper
    :param bool split: Should the results be split into perchrom vcfs?
    :return: Either the fsID to the genome-level vcf or a dict of results from running SomaticSniper
             on every chromosome
             perchrom_somaticsniper:
                 |- 'chr1': fsID
                 |- 'chr2' fsID
                 |
                 |-...
                 |
                 +- 'chrM': fsID
    :rtype: toil.fileStore.FileID|dict
    """
    # Get a list of chromosomes to handle
    if somaticsniper_options['chromosomes']:
        chromosomes = somaticsniper_options['chromosomes']
    else:
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
        unmerge_snipes = job.wrapJobFn(unmerge, filtersnipes.rv(), 'somaticsniper', chromosomes,
                                       somaticsniper_options, univ_options)
        filtersnipes.addChild(unmerge_snipes)
        return unmerge_snipes.rv()
    else:
        return filtersnipes.rv()


def run_somaticsniper_full(job, tumor_bam, normal_bam, univ_options, somaticsniper_options):
    """
    Run SomaticSniper on the DNA bams.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict somaticsniper_options: Options specific to SomaticSniper
    :return: fsID to the genome-level vcf
    :rtype: toil.fileStore.FileID
    """
    job.fileStore.logToMaster('Running SomaticSniper on %s' % univ_options['patient'])
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
                dockerhub=univ_options['dockerhub'], tool_version=somaticsniper_options['version'])
    outfile = job.fileStore.writeGlobalFile(output_file)
    return outfile


def filter_somaticsniper(job, tumor_bam, somaticsniper_output, tumor_pileup, univ_options,
                         somaticsniper_options):
    """
    Filter SomaticSniper calls.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param toil.fileStore.FileID somaticsniper_output: SomaticSniper output vcf
    :param toil.fileStore.FileID tumor_pileup: Pileup generated for the tumor bam
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict somaticsniper_options: Options specific to SomaticSniper
    :returns: fsID for the filtered genome-level vcf
    :rtype: toil.fileStore.FileID
    """
    job.fileStore.logToMaster('Filtering SomaticSniper for %s' % univ_options['patient'])
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
                dockerhub=univ_options['dockerhub'], tool_version=somaticsniper_options['version'])

    # Run prepare_for_readcount.pl
    parameters = ['prepare_for_readcount.pl',
                  '--snp-file', input_files['input.vcf'] + '.SNPfilter']
    # Creates /data/input.vcf.SNPfilter.pos
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=somaticsniper_options['version'])

    # Run  bam-readcount
    parameters = ['-b', '15',
                  '-f', input_files['genome.fa'],
                  '-l', input_files['input.vcf'] + '.SNPfilter.pos',
                  '-w', '1',
                  input_files['tumor.bam']]
    # Creates the read counts file
    with open(os.path.join(work_dir, 'readcounts.txt'), 'w') as readcounts_file:
        docker_call(tool='bam-readcount', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=readcounts_file,
                    tool_version=somaticsniper_options['bam_readcount']['version'])

    # Run fpfilter.pl
    parameters = ['fpfilter.pl',
                  '--snp-file', input_files['input.vcf'] + '.SNPfilter',
                  '--readcount-file', docker_path(readcounts_file.name)]

    # Creates input.vcf.SNPfilter.fp_pass and input.vcf.SNPfilter.fp_fail
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=somaticsniper_options['version'])

    # Run highconfidence.pl
    parameters = ['highconfidence.pl',
                  '--snp-file', input_files['input.vcf'] + '.SNPfilter.fp_pass']

    # Creates input.vcf.SNPfilter.fp_pass.hc
    docker_call(tool='somaticsniper-addons', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], tool_version=somaticsniper_options['version'])

    outfile = job.fileStore.writeGlobalFile(os.path.join(os.getcwd(),
                                                         'input.vcf.SNPfilter.fp_pass.hc'))
    return outfile


def process_somaticsniper_vcf(job, somaticsniper_vcf, work_dir, univ_options):
    """
    Process the SomaticSniper vcf for accepted calls.  Since the calls are pre-filtered, we just
    obtain the file from the file store and return it.

    :param toil.fileStore.FileID somaticsniper_vcf: fsID for a SomaticSniper generated chromosome
           vcf
    :param str work_dir: Working directory
    :param dict univ_options: Dict of universal options used by almost all tools
    :return: Path to the processed vcf
    :rtype: str
    """
    return job.fileStore.readGlobalFile(somaticsniper_vcf)


def run_pileup(job, tumor_bam, univ_options, somaticsniper_options):
    """
    Runs a samtools pileup on the tumor bam.

    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict somaticsniper_options: Options specific to SomaticSniper
    :return: fsID for the pileup file
    :rtype: toil.fileStore.FileID
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
        docker_call(tool='samtools', tool_parameters=parameters, work_dir=work_dir,
                    dockerhub=univ_options['dockerhub'], outfile=pileup_file,
                    tool_version=somaticsniper_options['samtools']['version'])
    outfile = job.fileStore.writeGlobalFile(pileup_file.name)
    return outfile
