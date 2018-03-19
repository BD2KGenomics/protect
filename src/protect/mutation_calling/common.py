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

from protect.common import chrom_sorted, export_results, get_files_from_filestore, untargz

import itertools
import os


def sample_chromosomes(job, genome_fai_file):
    """
    Get a list of chromosomes in the input data.

    :param toil.fileStore.FileID genome_fai_file: Job store file ID for the genome fai file
    :return: Chromosomes in the sample
    :rtype: list[str]
    """
    work_dir = os.getcwd()
    genome_fai = untargz(job.fileStore.readGlobalFile(genome_fai_file), work_dir)
    return chromosomes_from_fai(genome_fai)


def chromosomes_from_fai(genome_fai):
    """
    Read a fasta index (fai) file and parse the input chromosomes.

    :param str genome_fai: Path to the fai file.
    :return: list of input chromosomes
    :rtype: list[str]
    """
    chromosomes = []
    with open(genome_fai) as fai_file:
        for line in fai_file:
            line = line.strip().split()
            chromosomes.append(line[0])
    return chromosomes


def run_mutation_aggregator(job, mutation_results, univ_options):
    """
    Aggregate all the called mutations.

    :param dict mutation_results: Dict of dicts of the various mutation callers in a per chromosome
           format
    :param dict univ_options: Dict of universal options used by almost all tools
    :returns: fsID for the merged mutations file
    :rtype: toil.fileStore.FileID
    """
    # Setup an input data structure for the merge function
    out = {}
    for chrom in mutation_results['mutect'].keys():
        out[chrom] = job.addChildJobFn(merge_perchrom_mutations, chrom, mutation_results,
                                       univ_options).rv()
    merged_snvs = job.addFollowOnJobFn(merge_perchrom_vcfs, out, 'merged', univ_options)
    job.fileStore.logToMaster('Aggregated mutations for %s successfully' % univ_options['patient'])
    return merged_snvs.rv()


def merge_perchrom_mutations(job, chrom, mutations, univ_options):
    """
    Merge the mutation calls for a single chromosome.

    :param str chrom: Chromosome to process
    :param dict mutations: dict of dicts of the various mutation caller names as keys, and a dict of
           per chromosome job store ids for vcfs as value
    :param dict univ_options: Dict of universal options used by almost all tools
    :returns fsID for vcf contaning merged calls for the given chromosome
    :rtype: toil.fileStore.FileID
    """
    work_dir = os.getcwd()
    from protect.mutation_calling.muse import process_muse_vcf
    from protect.mutation_calling.mutect import process_mutect_vcf
    from protect.mutation_calling.radia import process_radia_vcf
    from protect.mutation_calling.somaticsniper import process_somaticsniper_vcf
    from protect.mutation_calling.strelka import process_strelka_vcf
    mutations.pop('indels')
    mutations['strelka_indels'] = mutations['strelka']['indels']
    mutations['strelka_snvs'] = mutations['strelka']['snvs']
    vcf_processor = {'snvs': {'mutect': process_mutect_vcf,
                              'muse': process_muse_vcf,
                              'radia': process_radia_vcf,
                              'somaticsniper': process_somaticsniper_vcf,
                              'strelka_snvs': process_strelka_vcf
                              },
                     'indels': {'strelka_indels': process_strelka_vcf
                                }
                     }
    #                 'fusions': lambda x: None,
    #                 'indels': lambda x: None}
    # For now, let's just say 2 out of n need to call it.
    # num_preds = len(mutations)
    # majority = int((num_preds + 0.5) / 2)
    majority = {'snvs': 2,
                'indels': 1}

    accepted_hits = defaultdict(dict)

    for mut_type in vcf_processor.keys():
        # Get input files
        perchrom_mutations = {caller: vcf_processor[mut_type][caller](job, mutations[caller][chrom],
                                                                      work_dir, univ_options)
                              for caller in vcf_processor[mut_type]}
        # Process the strelka key
        perchrom_mutations['strelka'] = perchrom_mutations['strelka_' + mut_type]
        perchrom_mutations.pop('strelka_' + mut_type)
        # Read in each file to a dict
        vcf_lists = {caller: read_vcf(vcf_file) for caller, vcf_file in perchrom_mutations.items()}
        all_positions = list(set(itertools.chain(*vcf_lists.values())))
        for position in sorted(all_positions):
            hits = {caller: position in vcf_lists[caller] for caller in perchrom_mutations.keys()}
            if sum(hits.values()) >= majority[mut_type]:
                callers = ','.join([caller for caller, hit in hits.items() if hit])
                assert position[1] not in accepted_hits[position[0]]
                accepted_hits[position[0]][position[1]] = (position[2], position[3], callers)

    with open(''.join([work_dir, '/', chrom, '.vcf']), 'w') as outfile:
        print('##fileformat=VCFv4.0', file=outfile)
        print('##INFO=<ID=callers,Number=.,Type=String,Description=List of supporting callers.',
              file=outfile)
        print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=outfile)
        for chrom in chrom_sorted(accepted_hits.keys()):
            for position in sorted(accepted_hits[chrom]):
                    print(chrom, position, '.', accepted_hits[chrom][position][0],
                          accepted_hits[chrom][position][1], '.', 'PASS',
                          'callers=' + accepted_hits[chrom][position][2], sep='\t', file=outfile)
    fsid = job.fileStore.writeGlobalFile(outfile.name)
    export_results(job, fsid, outfile.name, univ_options, subfolder='mutations/merged')
    return fsid


def read_vcf(vcf_file):
    """
    Read a vcf file to a dict of lists.

    :param str vcf_file: Path to a vcf file.
    :return: dict of lists of vcf records
    :rtype: dict
    """
    vcf_dict = []
    with open(vcf_file, 'r') as invcf:
        for line in invcf:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            vcf_dict.append((line[0], line[1], line[3], line[4]))
    return vcf_dict


def merge_perchrom_vcfs(job, perchrom_vcfs, tool_name, univ_options):
    """
    Merge per-chromosome vcf files into a single genome level vcf.

    :param dict perchrom_vcfs: Dictionary with chromosome name as key and fsID of the corresponding
           vcf as value
    :param str tool_name: Name of the tool that generated the vcfs
    :returns: fsID for the merged vcf
    :rtype: toil.fileStore.FileID
    """
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
    output_file = job.fileStore.writeGlobalFile(outvcf.name)
    export_results(job, output_file, outvcf.name, univ_options, subfolder='mutations/' + tool_name)
    job.fileStore.logToMaster('Ran merge_perchrom_vcfs for %s successfully' % tool_name)
    return output_file


def unmerge(job, input_vcf, tool_name, chromosomes, tool_options, univ_options):
    """
    Un-merge a vcf file into per-chromosome vcfs.

    :param str input_vcf: Input vcf
    :param str tool_name: The name of the mutation caller
    :param list chromosomes: List of chromosomes to retain
    :param dict tool_options: Options specific to the mutation caller
    :param dict univ_options: Dict of universal options used by almost all tools
    :return: dict of fsIDs, one for each chromosomal vcf
    :rtype: dict
    """
    work_dir = os.getcwd()
    input_files = {
        'input.vcf': input_vcf,
        'genome.fa.fai.tar.gz': tool_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    input_files['genome.fa.fai'] = untargz(input_files['genome.fa.fai.tar.gz'], work_dir)

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
    chroms = set(chromosomes).intersection(set(read_chromosomes.keys()))
    for chrom, chromvcf in read_chromosomes.items():
        chromvcf.close()
        if chrom not in chroms:
            continue
        outdict[chrom] = job.fileStore.writeGlobalFile(chromvcf.name)
        export_results(job, outdict[chrom], chromvcf.name, univ_options,
                       subfolder='mutations/' + tool_name)
    return outdict
