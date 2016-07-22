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
from protect.common import export_results, get_files_from_filestore, untargz

import itertools
import os


def sample_chromosomes(job, genome_fai_file):
    """
    Get a list of chromosomes in the input data

    :param job: job
    :param string genome_fai_file: Job store file ID for the genome fai file
    :returns list: Chromosomes in the sample
    """
    work_dir = os.getcwd()
    genome_fai = untargz(job.fileStore.readGlobalFile(genome_fai_file), work_dir)
    return chromosomes_from_fai(genome_fai)


def chromosomes_from_fai(genome_fai):
    chromosomes = []
    with open(genome_fai) as fai_file:
        for line in fai_file:
            line = line.strip().split()
            chromosomes.append(line[0])
    return chromosomes


def run_mutation_aggregator(job, mutation_results, univ_options):
    """
    This module will aggregate all the mutations called in the previous steps

    :param job: job
    :param dict mutation_results: Dict of dicts of the various mutation callers in a per chromosome
                                  format
    :param dict univ_options: Universal Options
    :param int numthreads: number of threads to use in the merge step
    :returns: job store file id fo rthe merged mutations file

    This module corresponds to node 15 on the tree
    """
    job.fileStore.logToMaster('Aggregating mutations for %s' % univ_options['patient'])
    # Setup an input data structure for the merge function

    out = {}
    for chrom in mutation_results['mutect'].keys():
        out[chrom] = job.addChildJobFn(merge_perchrom_mutations, chrom, mutation_results,
                                       univ_options).rv()
    merged_snvs = job.addFollowOnJobFn(merge_perchrom_vcfs, out, 'merged', univ_options)
    return merged_snvs.rv()


def merge_perchrom_mutations(job, chrom, mutations, univ_options):
    """
    This module will accept job store ids for vcf files for all snvs, and will merge the calls for a
    single provided chromosome.

    :param job: job
    :param str chrom: Chromosome to process
    :param dict mutations: dict of dicts of the various mutation caller names as keys, and a dict of
                           per chromosome job store ids for vcfs as value
    :param dict univ_options: Universal Options
    :returns dict of merged vcf
    """
    work_dir = os.getcwd()
    from protect.mutation_calling.muse import process_muse_vcf
    from protect.mutation_calling.mutect import process_mutect_vcf
    from protect.mutation_calling.radia import process_radia_vcf
    from protect.mutation_calling.somaticsniper import process_somaticsniper_vcf
    from protect.mutation_calling.strelka import process_strelka_vcf
    mutations.pop('indels')
    mutations.pop('fusions')
    mutations['strelka'] = mutations['strelka']['snvs']
    vcf_processor = {'mutect': process_mutect_vcf,
                     'muse': process_muse_vcf,
                     'radia': process_radia_vcf,
                     'somaticsniper': process_somaticsniper_vcf,
                     'strelka': process_strelka_vcf,
                     }
    #                 'fusions': lambda x: None,
    #                 'indels': lambda x: None}
    # For now, let's just say 2 out of n need to call it.
    # num_preds = len(mutations)
    # majority = int((num_preds + 0.5) / 2)
    majority = 2
    # Get input files
    perchrom_mutations = {caller: vcf_processor[caller](job, mutations[caller][chrom],
                                                        work_dir, univ_options)
                          for caller in mutations.keys()}

    # Read in each file to a dict
    vcf_lists = {caller: read_vcf(vcf_file) for caller, vcf_file in perchrom_mutations.items()}
    all_positions = list(set(itertools.chain(*vcf_lists.values())))
    with open(''.join([work_dir, '/', chrom, '.vcf']), 'w') as outfile:
        print('##fileformat=VCFv4.0', file=outfile)
        print('##INFO=<ID=callers,Number=.,Type=String,Description=List of supporting callers.',
              file=outfile)
        print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=outfile)
        for position in sorted(all_positions):
            hits = {caller: position in vcf_lists[caller] for caller in perchrom_mutations.keys()}
            if sum(hits.values()) >= majority:
                print(position[0], position[1], '.', position[2], position[3], '.', 'PASS',
                      'callers=' + ','.join([caller for caller, hit in hits.items() if hit]),
                      sep='\t', file=outfile)
    export_results(job, outfile.name, univ_options, subfolder='mutations/merged')
    outfile = job.fileStore.writeGlobalFile(outfile.name)
    return outfile


def read_vcf(vcf_file):
    """
    Reads a vcf file to a list of lists

    :param vcf_file:
    :return:
    """
    vcf_dict = []
    with open(vcf_file, 'r') as invcf:
        for line in invcf:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            vcf_dict.append((line[0], line[1], line[3], line[4]))
    return vcf_dict


def chrom_sorted(in_chroms):
    """
    This module will sort a list of chromosomes in the order 1..22, X, Y, M.
    :param list in_chroms: input chromsomes
    :return: sorted chromosomes
    :rtype: list
    """
    chr_prefix = False
    if in_chroms[0].startswith('chr'):
        in_chroms = [x.lstrip('chr') for x in in_chroms]
        chr_prefix = True
    assert in_chroms[0] in [str(x) for x in range(1,23)] + ['X', 'Y', 'M']
    in_chroms = sorted(in_chroms, key=lambda x: int(x) if x not in ('X', 'Y', 'M') else x)
    try:
        m_index = in_chroms.index('M')
    except ValueError:
        pass
    else:
        in_chroms.pop(m_index)
        in_chroms.append('M')
    # At this point it should be nicely sorted
    if chr_prefix:
        in_chroms = [''.join(['chr', x]) for x in in_chroms]
    return in_chroms


def merge_perchrom_vcfs(job, perchrom_vcfs, tool_name, univ_options):
    """
    This module will merge per-chromosome vcf files into a single genome level vcf.

    :param dict perchrom_vcfs: Dictionary with chromosome name as key and jobstore ID of
                               corresponding vcf as value
    :param str tool_name: Name of the tool that generated the vcfs

    :returns: Job Store File ID for the merged vcf
    """
    job.fileStore.logToMaster('Running merge_perchrom_vcfs  for %s' % tool_name)
    work_dir = os.getcwd()
    input_files = {''.join([chrom, '.vcf']): jsid for chrom, jsid in perchrom_vcfs.items()}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    first = True
    with open(''.join([work_dir, '/', 'all_merged.vcf']), 'w') as outvcf:
        for chromvcfname in chrom_sorted([x.rstrip('.vcf') for x in input_files.keys()]):
            with open(input_files[chromvcfname + '.vcf'], 'r') as infile:
                for line in infile:
                    line = line.strip()
                    if line.startswith('#'):
                        if first:
                            print(line, file=outvcf)
                        continue
                    first = False
                    print(line, file=outvcf)
    export_results(job, outvcf.name, univ_options, subfolder='mutations/' + tool_name)
    output_file = job.fileStore.writeGlobalFile(outvcf.name)
    return output_file


def unmerge(job, input_vcf, tool_name, tool_options, univ_options):
    """
    Un-merges a vcf file into a file per chromosome.

    :param str input_vcf: Input vcf
    :param str tool_name: The name of the mutation caller
    :param dict tool_options: Options specific to Somatic Sniper
    :param dict univ_options: Universal options
    :returns: dict of jsIDs, onr for each chromosomal vcf
    :rtype: dict
    """
    work_dir = os.getcwd()
    input_files = {
        'input.vcf': input_vcf,
        'genome.fa.fai.tar.gz': tool_options['genome_fai']}
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
        export_results(job, chromvcf.name, univ_options, subfolder='mutations/' + 'tool_name')
        outdict[chrom] = job.fileStore.writeGlobalFile(chromvcf.name)
    return outdict