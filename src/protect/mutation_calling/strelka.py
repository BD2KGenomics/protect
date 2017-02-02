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
from math import ceil
from protect.common import (docker_call,
                            docker_path,
                            get_files_from_filestore,
                            untargz)
from protect.mutation_calling.common import (sample_chromosomes,
                                             unmerge)
from toil.job import PromisedRequirement

import os


# disk for strelka.
def strelka_disk(tumor_bam, normal_bam, fasta):
    return int(ceil(tumor_bam.size) +
               ceil(normal_bam.size) +
               4 * ceil(fasta.size))


# Strelka is a different tool from the other callers because it is run in full
def run_strelka_with_merge(job, tumor_bam, normal_bam, univ_options, strelka_options):
    """
    This is a convenience function that runs the entire strelka sub-graph.
    """
    spawn = job.wrapJobFn(run_strelka, tumor_bam, normal_bam, univ_options,
                          strelka_options, split=False).encapsulate()
    job.addChild(spawn)
    return spawn.rv()


def run_strelka(job, tumor_bam, normal_bam, univ_options, strelka_options, split=True):
    """
    This module will spawn a strelka job for each chromosome on the DNA bams.

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
    4. strelka_options: Dict of parameters specific to strelka
         strelka_options
              |- 'dbsnp_vcf': <JSid for dnsnp vcf file>
              |- 'dbsnp_idx': <JSid for dnsnp vcf index file>
              |- 'cosmic_vcf': <JSid for cosmic vcf file>
              |- 'cosmic_idx': <JSid for cosmic vcf index file>
              |- 'genome_fasta': <JSid for genome fasta file>
              +- 'genome_dict': <JSid for genome fasta dict file>
              +- 'genome_fai': <JSid for genome fasta index file>

    RETURN VALUES
    1. perchrom_strelka: Dict of results of strelka per chromosome
         perchrom_strelka
              |- 'chr1'
              |   +- 'strelka_chr1.vcf': <JSid>
              |   +- 'strelka_chr1.out': <JSid>
              |- 'chr2'
              |   |- 'strelka_chr2.vcf': <JSid>
              |   +- 'strelka_chr2.out': <JSid>
             etc...

    This module corresponds to node 11 on the tree
    """
    chromosomes = sample_chromosomes(job, strelka_options['genome_fai'])
    num_cores = min(len(chromosomes), univ_options['max_cores'])
    strelka = job.wrapJobFn(run_strelka_full, tumor_bam, normal_bam, univ_options,
                            strelka_options,
                            disk=PromisedRequirement(strelka_disk,
                                                     tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                                     normal_bam['normal_dna_fix_pg_sorted.bam'],
                                                     strelka_options['genome_fasta']),
                            memory='6G',
                            cores=num_cores)
    job.addChild(strelka)
    if split:
        unmerge_strelka = job.wrapJobFn(wrap_unmerge, strelka.rv(), strelka_options, univ_options
                                        ).encapsulate()
        strelka.addChild(unmerge_strelka)
        return unmerge_strelka.rv()
    else:
        return strelka.rv()


def run_strelka_full(job, tumor_bam, normal_bam, univ_options, strelka_options):
    """
    This module will run strelka on the DNA bams.

    ARGUMENTS
    :param dict tumor_bam: REFER ARGUMENTS of spawn_strelka()
    :param dict normal_bam: REFER ARGUMENTS of spawn_strelka()
    :param dict univ_options: REFER ARGUMENTS of spawn_strelka()
    :param dict strelka_options: REFER ARGUMENTS of spawn_strelka()

    RETURN VALUES
    :returns: dict of output vcfs for each chromosome
    :rtype: dict
    """
    job.fileStore.logToMaster('Running strelka on %s' % univ_options['patient'])
    work_dir = os.getcwd()
    input_files = {
        'tumor.bam': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
        'tumor.bam.bai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
        'normal.bam': normal_bam['normal_dna_fix_pg_sorted.bam'],
        'normal.bam.bai': normal_bam['normal_dna_fix_pg_sorted.bam.bai'],
        'genome.fa.tar.gz': strelka_options['genome_fasta'],
        'genome.fa.fai.tar.gz': strelka_options['genome_fai'],
        'config.ini.tar.gz': strelka_options['strelka_config']
    }
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai', 'config.ini'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    parameters = [input_files['config.ini'],
                  input_files['tumor.bam'],
                  input_files['normal.bam'],
                  input_files['genome.fa'],
                  str(job.cores)
                  ]
    docker_call(tool='strelka:1.0.15', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'])
    output_dict = {}
    for mutation_type in ['snvs', 'indels']:
        output_dict[mutation_type] = job.fileStore.writeGlobalFile(os.path.join(
            work_dir, 'strelka_out', 'results', 'passed.somatic.' + mutation_type + '.vcf'))
    return output_dict


def process_strelka_vcf(job, strelka_vcf, work_dir, univ_options):
    """
    This will process all the strelka vcfs to return only passing calls
    :param job: job
    :param str strelka_vcf: Job Store ID corresponding to a strelka vcf for 1 chromosome
    :param univ_options: Universal options
    :returns dict: Dict with chromosomes as keys and path to the corresponding parsed strelka vcfs as
                   values
    """
    return job.fileStore.readGlobalFile(strelka_vcf)


def wrap_unmerge(job, strelka_out, strelka_options, univ_options):
    """
    This is a convenience function that wraps the unmerge for strelkas snvs and indels

    :param toil.Job job: job
    :param dict strelka_out:
    :param dict strelka_options:
    :param dict univ_options:
    """
    return {'snvs': job.addChildJobFn(unmerge, strelka_out['snvs'], 'strelka', strelka_options,
                                      univ_options).rv(),
            'indels': job.addChildJobFn(unmerge, strelka_out['indels'], 'strelka', strelka_options,
                                        univ_options).rv()}
