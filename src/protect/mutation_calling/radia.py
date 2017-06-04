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
                            export_results,
                            get_files_from_filestore,
                            untargz)
from protect.mutation_calling.common import sample_chromosomes, merge_perchrom_vcfs
from toil.job import PromisedRequirement

import os
import sys


# disk for radia and filterradia.
def radia_disk(tumor_bam, normal_bam, rna_bam, fasta):
    return int(ceil(tumor_bam.size) +
               ceil(normal_bam.size) +
               ceil(rna_bam.size) +
               5 * ceil(fasta.size))


def run_radia_with_merge(job, rna_bam, tumor_bam, normal_bam, univ_options, radia_options):
    """
    A wrapper for the the entire RADIA sub-graph.

    :param dict rna_bam: Dict dicts of bam and bai for tumor RNA-Seq obtained by running STAR within
           ProTECT.
    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict radia_options: Options specific to RADIA
    :return: fsID to the merged RADIA calls
    :rtype: toil.fileStore.FileID
    """
    spawn = job.wrapJobFn(run_radia, rna_bam['rna_genome'], tumor_bam,
                          normal_bam, univ_options, radia_options, disk='100M',
                          memory='100M').encapsulate()
    merge = job.wrapJobFn(merge_perchrom_vcfs, spawn.rv(), univ_options, disk='100M', memory='100M')
    job.addChild(spawn)
    spawn.addChild(merge)
    return merge.rv()


def run_radia(job, rna_bam, tumor_bam, normal_bam, univ_options, radia_options):
    """
    Spawn a RADIA job for each chromosome on the input bam trios.

    :param dict rna_bam: Dict of bam and bai for tumor DNA-Seq.  It can be one of two formats
           rna_bam:   # Just the genomic bam and bai
                |- 'rna_genome_sorted.bam': fsID
                +- 'rna_genome_sorted.bam.bai': fsID
           OR
           rna_bam:   # The output from run_star
               |- 'rna_transcriptome.bam': fsID
               |- 'rna_genome':     # Only this part will be used
                       |- 'rna_genome_sorted.bam': fsID
                       +- 'rna_genome_sorted.bam.bai': fsID
    :param dict tumor_bam: Dict of bam and bai for tumor DNA-Seq
    :param dict normal_bam: Dict of bam and bai for normal DNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict radia_options: Options specific to RADIA
    :return: Dict of results from running RADIA on every chromosome
             perchrom_radia:
                 |- 'chr1': fsID
                 |- 'chr2' fsID
                 |
                 |-...
                 |
                 +- 'chrM': fsID
    :rtype: dict
    """
    job.fileStore.logToMaster('Running spawn_radia on %s' % univ_options['patient'])
    if 'rna_genome' in rna_bam.keys():
        rna_bam = rna_bam['rna_genome']
    elif set(rna_bam.keys()) == {'rna_genome_sorted.bam', 'rna_genome_sorted.bam.bai'}:
        pass
    else:
        raise RuntimeError('An improperly formatted dict was passed to rna_bam.')
        
    bams = {'tumor_rna': rna_bam['rna_genome_sorted.bam'],
            'tumor_rnai': rna_bam['rna_genome_sorted.bam.bai'],
            'tumor_dna': tumor_bam['tumor_dna_fix_pg_sorted.bam'],
            'tumor_dnai': tumor_bam['tumor_dna_fix_pg_sorted.bam.bai'],
            'normal_dna': normal_bam['normal_dna_fix_pg_sorted.bam'],
            'normal_dnai': normal_bam['normal_dna_fix_pg_sorted.bam.bai']}
    # Get a list of chromosomes to process
    if radia_options['chromosomes']:
        chromosomes = radia_options['chromosomes']
    else:
        chromosomes = sample_chromosomes(job, radia_options['genome_fai'])
    perchrom_radia = defaultdict()
    for chrom in chromosomes:
        radia = job.addChildJobFn(run_radia_perchrom, bams, univ_options, radia_options, chrom,
                                  memory='6G',
                                  disk=PromisedRequirement(
                                      radia_disk, tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                      normal_bam['normal_dna_fix_pg_sorted.bam'],
                                      rna_bam['rna_genome_sorted.bam'],
                                      radia_options['genome_fasta']))
        filter_radia = radia.addChildJobFn(run_filter_radia, bams, radia.rv(), univ_options,
                                           radia_options, chrom, memory='6G',
                                           disk=PromisedRequirement(
                                               radia_disk, tumor_bam['tumor_dna_fix_pg_sorted.bam'],
                                               normal_bam['normal_dna_fix_pg_sorted.bam'],
                                               rna_bam['rna_genome_sorted.bam'],
                                               radia_options['genome_fasta']))
        perchrom_radia[chrom] = filter_radia.rv()
    return perchrom_radia


def run_radia_perchrom(job, bams, univ_options, radia_options, chrom):
    """
    Run RADIA call on a single chromosome in the input bams.

    :param dict bams: Dict of bam and bai for tumor DNA-Seq, normal DNA-Seq and tumor RNA-Seq
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict radia_options: Options specific to RADIA
    :param str chrom: Chromosome to process
    :return: fsID for the chromsome vcf
    :rtype: toil.fileStore.FileID
    """
    job.fileStore.logToMaster('Running radia on %s:%s' % (univ_options['patient'], chrom))
    work_dir = os.getcwd()
    input_files = {
        'rna.bam': bams['tumor_rna'],
        'rna.bam.bai': bams['tumor_rnai'],
        'tumor.bam': bams['tumor_dna'],
        'tumor.bam.bai': bams['tumor_dnai'],
        'normal.bam': bams['normal_dna'],
        'normal.bam.bai': bams['normal_dnai'],
        'genome.fa.tar.gz': radia_options['genome_fasta'],
        'genome.fa.fai.tar.gz': radia_options['genome_fai']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    input_files = {key: docker_path(path) for key, path in input_files.items()}

    radia_output = ''.join([work_dir, '/radia_', chrom, '.vcf'])
    radia_log = ''.join([work_dir, '/radia_', chrom, '_radia.log'])
    parameters = [univ_options['patient'],  # shortID
                  chrom,
                  '-n', input_files['normal.bam'],
                  '-t', input_files['tumor.bam'],
                  '-r', input_files['rna.bam'],
                  ''.join(['--rnaTumorFasta=', input_files['genome.fa']]),
                  '-f', input_files['genome.fa'],
                  '-o', docker_path(radia_output),
                  '-i', univ_options['ref'],
                  '-m', input_files['genome.fa'],
                  '-d', 'aarjunrao@soe.ucsc.edu',
                  '-q', 'Illumina',
                  '--disease', 'CANCER',
                  '-l', 'INFO',
                  '-g', docker_path(radia_log)]
    docker_call(tool='radia', tool_parameters=parameters,
                work_dir=work_dir, dockerhub=univ_options['dockerhub'],
                tool_version=radia_options['version'])
    output_file = job.fileStore.writeGlobalFile(radia_output)
    return output_file


def run_filter_radia(job, bams, radia_file, univ_options, radia_options, chrom):
    """
    Run filterradia on the RADIA output.

    :param dict bams: Dict of bam and bai for tumor DNA-Seq, normal DNA-Seq and tumor RNA-Seq
    :param toil.fileStore.FileID radia_file: The vcf from runnning RADIA
    :param dict univ_options: Dict of universal options used by almost all tools
    :param dict radia_options: Options specific to RADIA
    :param str chrom: Chromosome to process
    :return: fsID for the filtered chromsome vcf
    :rtype: toil.fileStore.FileID
    """
    job.fileStore.logToMaster('Running filter-radia on %s:%s' % (univ_options['patient'], chrom))
    work_dir = os.getcwd()
    input_files = {
        'rna.bam': bams['tumor_rna'],
        'rna.bam.bai': bams['tumor_rnai'],
        'tumor.bam': bams['tumor_dna'],
        'tumor.bam.bai': bams['tumor_dnai'],
        'normal.bam': bams['normal_dna'],
        'normal.bam.bai': bams['normal_dnai'],
        'radia.vcf': radia_file,
        'genome.fa.tar.gz': radia_options['genome_fasta'],
        'genome.fa.fai.tar.gz': radia_options['genome_fai'],
        'cosmic_beds': radia_options['cosmic_beds'],
        'dbsnp_beds': radia_options['dbsnp_beds'],
        'retrogene_beds': radia_options['retrogene_beds'],
        'pseudogene_beds': radia_options['pseudogene_beds'],
        'gencode_beds': radia_options['gencode_beds']
    }
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    for key in ('genome.fa', 'genome.fa.fai'):
        input_files[key] = untargz(input_files[key + '.tar.gz'], work_dir)
    for key in ('cosmic_beds', 'dbsnp_beds', 'retrogene_beds', 'pseudogene_beds', 'gencode_beds'):
        input_files[key] = untargz(input_files[key], work_dir)

    input_files = {key: docker_path(path) for key, path in input_files.items()}

    filterradia_log = ''.join([work_dir, '/radia_filtered_', chrom, '_radia.log'])
    parameters = [univ_options['patient'],  # shortID
                  chrom.lstrip('chr'),
                  input_files['radia.vcf'],
                  '/data',
                  '/home/radia/scripts',
                  '-d', input_files['dbsnp_beds'],
                  '-r', input_files['retrogene_beds'],
                  '-p', input_files['pseudogene_beds'],
                  '-c', input_files['cosmic_beds'],
                  '-t', input_files['gencode_beds'],
                  '--noSnpEff',
                  '--noBlacklist',
                  '--noTargets',
                  '--noRnaBlacklist',
                  '-f', input_files['genome.fa'],
                  '--log=INFO',
                  '-g', docker_path(filterradia_log)]
    docker_call(tool='filterradia',
                tool_parameters=parameters, work_dir=work_dir, dockerhub=univ_options['dockerhub'],
                tool_version=radia_options['version'])
    output_file = ''.join([work_dir, '/', chrom, '.vcf'])
    os.rename(''.join([work_dir, '/', univ_options['patient'], '_', chrom, '.vcf']), output_file)
    output_fsid = job.fileStore.writeGlobalFile(output_file)
    export_results(job, output_fsid, output_file, univ_options, subfolder='mutations/radia')
    return output_fsid


def process_radia_vcf(job, radia_vcf, work_dir, univ_options):
    """
    Process the RADIA vcf to for passing calls and additionally sites having multiple alt alleles
    to pick out on the most likely ones.

    :param toil.fileStore.FileID radia_vcf: fsID for a RADIA generated chromosome vcf
    :param str work_dir: Working directory
    :param dict univ_options: Dict of universal options used by almost all tools
    :return: Path to the processed vcf
    :rtype: str
    """
    radia_vcf = job.fileStore.readGlobalFile(radia_vcf)
    with open(radia_vcf, 'r') as infile, open(radia_vcf + 'radia_parsed.tmp', 'w') as outfile:
        # The columns in INFILE are
        # [0] CHROM
        # [1] POS
        # [2] ID
        # [3] REF
        # [4] ALT
        # [5] QUAL
        # [6] FILTER
        # [7] INFO
        # [8] FORMAT
        # [9] DNA_NORMAL
        # [10] DNA_TUMOR
        # [11] RNA_TUMOR  -  Not always present
        for line in infile:
            # Print header to outfile
            if line.startswith('#'):
                print(line.strip(), file=outfile)
                continue
            line = line.strip().split('\t')
            # If the call was not PASSing, or if the call was germline: skip
            if line[6] != 'PASS' or 'MT=GERM' in line[7]:
                continue
            # If there is just 1 ALT allele, print and continue
            if len(line[4]) == 1:
                print('\t'.join(line), file=outfile)
            # If not, process
            else:
                seq_field_indeces = [9, 10]
                alleles = [line[3]] + line[4].split(',')  # all alleles, incl. REF
                # collect tumor, normal and (if present) rna AD and AFs
                # AD = Depth of reads supporting each allele
                # AF = Fraction of reads supporting each allele
                # normal_ad = line[9].split(':')[5].split(',')
                normal_af = line[9].split(':')[6].split(',')
                tumor_ad = line[10].split(':')[5].split(',')
                tumor_af = line[10].split(':')[6].split(',')
                if len(line[11]) > 1:
                    rna_ad = line[11].split(':')[5].split(',')
                    rna_af = line[11].split(':')[6].split(',')
                    seq_field_indeces += [11]  # append rna since it is present
                else:
                    # If rna is missing, set RNA_AD and RNA_AF to null sets for easily
                    # integrating into the logic in the following code
                    rna_ad = rna_af = [0, 0, 0, 0]
                # Initialise variables to store the probable ALT alleles and the index values of
                # the same wrt AD and AF
                out_alleles = set([])
                out_af_ad_index = {0}
                # parse AD and AF to get most probable ALT alleles
                for i in range(1, len(normal_af)):
                    # Criteria for selection = AD > 4 and AF >0.1 in either tumor or RNA, given
                    # normal AF < 0.1
                    if ((float(tumor_af[i]) >= 0.1 and int(tumor_ad[i]) >= 4) or
                            (float(rna_af[i]) >= 0.1 and int(rna_ad[i]) >= 4)) and \
                            (float(normal_af[i]) < 0.1):
                        out_alleles.add(alleles[i])
                        out_af_ad_index.add(i)
                # If the number of probable alleles is greater than 0 the print to outfile with
                # the modified allele fraction representing reads corrresponding to all alleles
                if len(out_alleles) > 0:
                    line[4] = ','.join(out_alleles)  # set alt alleles
                    # Modify the AD and AF values in the TUMOR/NORMAL/RNA fields
                    # one at a time.  Seq fields contain
                    # [0] GT* - Genotype
                    # [1] DP - Read depth at this position in the sample
                    # [2] INDEL - Number of indels
                    # [3] START - Number of reads starting at this position
                    # [4] STOP - Number of reads stopping at this position
                    # [5] AD* - Depth of reads supporting alleles
                    # [6] AF* - Fraction of reads supporting alleles
                    # [7] BQ* - Avg base quality for reads supporting alleles
                    # [8] SB* - Strand Bias for reads supporting alleles
                    # Fields marked with *s are teh ones that contain info for each seq field
                    # and need to be modified
                    for seq_field_index in seq_field_indeces:
                        # Get the details for seq_field
                        deets = line[seq_field_index].split(':')
                        # modify fields 5 thu 8 to hold only info for the probable
                        # alleles
                        for field_index in range(5, 9):
                            field = deets[field_index].split(",")
                            deets[field_index] = ",".join([x for i, x in enumerate(field)
                                                           if i in out_af_ad_index])
                        # Modify DP to hold the new total of reads
                        deets[1] = str(sum([int(x) for x in deets[5].split(",")]))
                        # get the most likely genotypes based on AD and AF
                        gt_by_ad = set([i for i, x in enumerate(deets[5].split(","))
                                        if int(x) >= 4])
                        gt_by_af = set([i for i, x in enumerate(deets[6].split(","))
                                        if float(x) >= 0.1])
                        # Get the consensus genotype
                        genotype = gt_by_ad.intersection(gt_by_af)
                        if len(genotype) == 0:
                            deets[0] = "0/0"
                        elif len(genotype) == 1:
                            deets[0] = "/".join([str(x) for x in genotype] +
                                                [str(x) for x in genotype])
                        elif len(genotype) == 2:
                            deets[0] = "/".join([str(x) for x in genotype])
                        else:
                            print("ERROR : triple genotype detected", file=sys.stderr)
                            print(line, file=sys.stdout)
                        # Rejoin the details line
                        line[seq_field_index] = ":".join(deets)
                    # Print the modified line to output
                    print("\t".join(line), file=outfile)
                # Else do nothing
                else:
                    pass
    return outfile.name
