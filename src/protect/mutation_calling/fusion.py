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

from bd2k.util.expando import Expando
from toil.job import PromisedRequirement

import collections
import csv
import logging
import os
import re
import string

from protect.common import docker_call, export_results, get_files_from_filestore, untargz


def fusion_disk(rna_fastqs, star_tar):
    return int(4 * ceil(sum([f.size for f in rna_fastqs]) + 524288) +
               2 * ceil(star_tar.size + 524288) +
               5242880)


def wrap_fusion(job,
                fastqs,
                star_output,
                univ_options,
                star_fusion_options,
                fusion_inspector_options):
    """
    A wrapper for run_fusion using the results from cutadapt and star as input.

    :param tuple fastqs: RNA-Seq FASTQ Filestore IDs
    :param dict star_output: Dictionary containing STAR output files
    :param dict univ_options: universal arguments used by almost all tools
    :param dict star_fusion_options: STAR-Fusion specific parameters
    :param dict fusion_inspector_options: FusionInspector specific parameters
    :return: Transgene BEDPE file
    :rtype: toil.fileStore.FileID

    """
    # Give user option to skip fusion calling
    if not star_fusion_options['run']:
        job.fileStore.logToMaster('Skipping STAR-Fusion on %s' % univ_options['patient'])
        return

    fusion = job.wrapJobFn(run_fusion, fastqs, star_output['rnaChimeric.out.junction'],
                               univ_options, star_fusion_options, fusion_inspector_options,
                               cores=star_fusion_options['n'],
                               memory=PromisedRequirement(lambda x: int(1.85 * x.size),
                                                          star_fusion_options['index']),
                               disk=PromisedRequirement(fusion_disk,
                                                        fastqs,
                                                        star_fusion_options['index'])).encapsulate()
    job.addChild(fusion)
    return fusion.rv()


def run_fusion(job,
               fastqs,
               junction_file,
               univ_options,
               star_fusion_options,
               fusion_inspector_options):
    """
    Runs STAR-Fusion and filters fusion calls using FusionInspector

    :param tuple fastqs: RNA-Seq FASTQ Filestore IDs
    :param toil.fileStore.FileID junction_file: Chimeric junction file
    :param dict univ_options: universal arguments used by almost all tools
    :param dict star_fusion_options: STAR-Fusion specific parameters
    :return: Transgene BEDPE file
    :rtype: toil.fileStore.FileID
    """
    work_dir = job.fileStore.getLocalTempDir()

    input_files = {'rna_1.fq.gz': fastqs[0],
                   'rna_2.fq.gz': fastqs[1],
                   'tool_index.tar.gz': star_fusion_options['index']}

    parameters = []

    # If there isn't a junction file, then we can run STAR-Fusion from the fastq files
    if junction_file:
        input_files['STAR.junction'] = junction_file
        parameters.extend(['--chimeric_junction', '/data/STAR.junction'])

    else:
        parameters.extend(['--left_fq', '/data/rna_1.fq.gz', '--right_fq', '/data/rna_2.fq.gz'])

    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    input_files['tool_index'] = os.path.basename(untargz(input_files['tool_index.tar.gz'],
                                                         work_dir))

    cores = star_fusion_options['n']
    parameters.extend(['--output_dir', '/data/fusion-output',
                       '--genome_lib_dir', input_files['tool_index'],
                       '--CPU', str(cores)])

    docker_call(tool='star-fusion',
                tool_parameters=parameters,
                work_dir=work_dir,
                dockerhub=univ_options['dockerhub'],
                tool_version=star_fusion_options['version'])

    star_output = 'fusion-output/star-fusion.fusion_candidates.final.abridged'
    fusion_path = os.path.join(work_dir, star_output)

    # Export the STAR-Fusion predictions
    export_results(job,
                   job.fileStore.writeGlobalFile(fusion_path),
                   'star-fusion-predictions.tsv',
                   univ_options, subfolder='mutations/fusions')

    # Check for fusion prediction
    with open(fusion_path, 'r') as f:
        # Skip header
        f.next()
        try:
            f.next()
        except StopIteration:
            logging.warning('%s: Did not find any fusions!' % univ_options['patient'])
            return

    parameters = ['--fusions', '/data/%s' % star_output,
                  '--genome_lib', input_files['tool_index'],
                  '--left_fq', '/data/rna_1.fq.gz',
                  '--right_fq', '/data/rna_2.fq.gz',
                  '--out_dir', '/data/FusionInspector',
                  '--out_prefix', 'FusionInspector',
                  '--CPU', str(cores)]

    if fusion_inspector_options['run_trinity']:
        parameters.append('--include_Trinity')

    docker_call(tool='fusion-inspector',
                tool_parameters=parameters,
                work_dir=work_dir,
                dockerhub=univ_options['dockerhub'],
                tool_version=fusion_inspector_options['version'])

    found_fusion = False
    inspector_output = 'FusionInspector/FusionInspector.fusion_predictions.final.abridged.FFPM'
    fusion_path = os.path.join(work_dir, inspector_output)
    output_path = os.path.join(work_dir, 'fusion.final')

    # Export the FusionInpsector predictions
    export_results(job,
                   job.fileStore.writeGlobalFile(fusion_path),
                   'fusion-inspector-predictions.tsv',
                   univ_options, subfolder='mutations/fusions')

    # Remove fusions without a large anchor sequence and at least 0.1
    # fusion fragments per million reads
    if os.path.exists(fusion_path):
        with open(fusion_path, 'r') as f, open(output_path, 'w') as g:
            g.write(f.next())
            for line in f:
                fields = line.strip().split()

                # Check for a large anchor support
                ldas = fields[10]

                assert ldas in {'YES', 'NO'}, 'FusionInpsector file is malformed!'

                j_ffpm, s_ffpm = fields[-2:]

                # Fusions without a larger anchor support or low read support
                # are suspicious and should not be consider for further analysis
                if ldas == 'YES' and sum([float(j_ffpm), float(s_ffpm)]) > 0.1:
                    found_fusion = True
                    g.write(line)

    if found_fusion:
        fusion_bed_f = 'FusionInspector/FusionInspector.bed'
        fusion_bed_path = os.path.join(work_dir, fusion_bed_f)
        transcript_f = 'FusionInspector/FusionInspector.gmap_trinity_GG.fusions.fasta'
        transcript_path = os.path.join(work_dir, transcript_f)
        transcript_gff_f = 'FusionInspector/FusionInspector.gmap_trinity_GG.fusions.gff3'
        transcript_gff_path = os.path.join(work_dir, transcript_gff_f)

        transcripts = None
        transcript_annotation = None

        if os.path.exists(transcript_path):
            transcripts = job.fileStore.writeGlobalFile(transcript_path)
            export_results(job,
                           transcripts,
                           transcript_path,
                           univ_options,
                           subfolder='mutations/fusions')

        if os.path.exists(transcript_gff_path):
            transcript_annotation = job.fileStore.writeGlobalFile(transcript_gff_path)
            export_results(job,
                           transcript_annotation,
                           transcript_gff_path,
                           univ_options,
                           subfolder='mutations/fusions')

        fusion_annotation = job.fileStore.writeGlobalFile(fusion_bed_path)
        filtered_fusions = job.fileStore.writeGlobalFile(output_path)

        export_results(job,
                       filtered_fusions,
                       output_path,
                       univ_options,
                       subfolder='mutations/fusions')
        job.fileStore.logToMaster('Ran STAR-Fusion on %s successfully' % univ_options['patient'])
        return job.addChildJobFn(reformat_star_fusion_output,
                                 fusion_annotation,
                                 filtered_fusions,
                                 transcripts,
                                 transcript_annotation,
                                 univ_options).rv()

    else:
        job.fileStore.logToMaster('No fusions detected for %s' % univ_options['patient'])


def parse_star_fusion(infile):
    """
    Parses STAR-Fusion format and returns an Expando object with basic features

    :param str infile: path to STAR-Fusion prediction file
    :return: Fusion prediction attributes
    :rtype: bd2k.util.expando.Expando
    """
    reader = csv.reader(infile, delimiter='\t')
    header = reader.next()
    header = {key: index for index, key in enumerate(header)}

    features = ['LeftGene', 'LeftLocalBreakpoint', 'LeftBreakpoint',
                'RightGene', 'RightLocalBreakpoint', 'RightBreakpoint',
                'LargeAnchorSupport', 'JunctionReadCount', 'SpanningFragCount']

    for line in reader:
        yield Expando(dict((feature, line[header[feature]]) for feature in features))


def get_transcripts(transcript_file):
    """
    Parses FusionInspector transcript file and returns dictionary of sequences

    :param str transcript_file: path to transcript FASTA
    :return: de novo assembled transcripts
    :rtype: dict
    """
    with open(transcript_file, 'r') as fa:
        transcripts = {}
        regex_s = r"(?P<ID>TRINITY.*)\s(?P<fusion>.*--.*):(?P<left_start>\d+)-(?P<right_start>\d+)"
        regex = re.compile(regex_s)
        while True:
            # Usually the transcript is on one line
            try:
                info = fa.next()
                seq = fa.next()

                assert info.startswith('>')

                m = regex.search(info)
                if m:
                    transcripts[m.group('ID')] = seq.strip()

            except StopIteration:
                break

            except AssertionError:
                print("WARNING: Malformed fusion transcript file")
    return transcripts


def split_fusion_transcript(annotation_path, transcripts):
    """
    Finds the breakpoint in the fusion transcript and splits the 5' donor from the 3' acceptor

    :param str annotation_path: Path to transcript annotation file
    :param dict transcripts: Dictionary of fusion transcripts
    :return: 5' donor sequences and 3' acceptor sequences
    :rtype: tuple
    """
    annotation = collections.defaultdict(dict)

    forward = 'ACGTN'
    reverse = 'TGCAN'
    trans = string.maketrans(forward, reverse)

    # Pull in assembled transcript annotation
    five_pr_splits = collections.defaultdict(dict)
    three_pr_splits = collections.defaultdict(dict)

    regex = re.compile(r'ID=(?P<ID>.*);Name=(?P<Name>.*);Target=(?P<Target>.*)\s(?P<start>\d+)\s(?P<stop>\d+)')

    with open(annotation_path, 'r') as gff:
        for line in gff:
            print(line)
            if line.startswith('#'):
                _, eyd, fusion = line.strip().split()
                fusion, start_stop = fusion.split(':')
                left_break, right_break = start_stop.split('-')

                annotation[fusion][eyd] = {}
                annotation[fusion][eyd]['left_break'] = left_break
                annotation[fusion][eyd]['right_break'] = right_break

            else:
                line = line.strip().split('\t')
                fusion = line[0]
                strand = line[6]
                block_start = line[3]
                block_stop = line[4]
                attr = line[8]
                m = regex.search(attr)
                if m:
                    transcript_id = m.group('Name')

                    rb = any([block_start == annotation[fusion][transcript_id]['right_break'],
                              block_stop == annotation[fusion][transcript_id]['right_break']])

                    lb = any([block_start == annotation[fusion][transcript_id]['left_break'],
                              block_stop == annotation[fusion][transcript_id]['left_break']])

                    if strand == '-' and rb:
                        transcript_split = int(m.group('stop')) + 1   # Off by one
                        # Take the reverse complement to orient transcripts from 5' to 3'
                        five_seq = transcripts[transcript_id][transcript_split:]
                        five_pr_splits[fusion][transcript_id] = five_seq.translate(trans)[::-1]
                        three_seq = transcripts[transcript_id][:transcript_split]
                        three_pr_splits[fusion][transcript_id] = three_seq.translate(trans)[::-1]

                    elif strand == '+' and lb:
                        transcript_split = int(m.group('stop'))
                        s1 = transcripts[transcript_id][:transcript_split]
                        five_pr_splits[fusion][transcript_id] = s1
                        s2 = transcripts[transcript_id][transcript_split:]
                        three_pr_splits[fusion][transcript_id] = s2

    return five_pr_splits, three_pr_splits


def get_gene_ids(fusion_bed):
    """
    Parses FusionInspector bed file to ascertain the ENSEMBL gene ids

    :param str fusion_bed: path to fusion annotation
    :return: dict
    """
    with open(fusion_bed, 'r') as f:
        gene_to_id = {}
        regex = re.compile(r'(?P<gene>ENSG\d*)')
        for line in f:
            line = line.split('\t')
            transcript, gene_bit, name = line[3].split(';')
            m = regex.search(gene_bit)
            if m:
                gene_to_id[name] = m.group('gene')
    return gene_to_id


def reformat_star_fusion_output(job,
                                fusion_annot,
                                fusion_file,
                                transcript_file,
                                transcript_gff_file,
                                univ_options):
    """
    Writes STAR-Fusion results in Transgene BEDPE format

    :param toil.fileStore.FileID fusion_annot: Fusion annotation
    :param toil.fileStore.FileID fusion_file: STAR-fusion prediction file
    :param toil.fileStore.FileID transcript_file: Fusion transcript FASTA file
    :param toil.fileStore.FileID transcript_gff_file: Fusion transcript GFF file
    :param dict univ_options: universal arguments used by almost all tools
    :return: Transgene BEDPE file
    :rtype: toil.fileStore.FileID
    """
    input_files = {'results.tsv': fusion_file,
                   'fusion.bed': fusion_annot}

    if transcript_file and transcript_gff_file:
        input_files['transcripts.fa'] = transcript_file
        input_files['transcripts.gff'] = transcript_gff_file

    work_dir = job.fileStore.getLocalTempDir()
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    # Pull in assembled transcript file
    hugo_to_gene_ids = get_gene_ids(input_files['fusion.bed'])

    if transcript_file and transcript_gff_file:
        transcripts = get_transcripts(input_files['transcripts.fa'])
        five_pr_splits, three_pr_splits = split_fusion_transcript(input_files['transcripts.gff'],
                                                                  transcripts)

    else:
        five_pr_splits = collections.defaultdict(dict)
        three_pr_splits = collections.defaultdict(dict)

    # Pull in assembled transcript annotation

    # Header for BEDPE file
    header = ['# chr1', 'start1', 'end1',
              'chr2', 'start2', 'end2',
              'name', 'score',
              'strand1', 'strand2',
              'junctionSeq1', 'junctionSeq2',
              'hugo1', 'hugo2']

    output_path = os.path.join(work_dir, 'fusion_results.bedpe')
    with open(input_files['results.tsv'], 'r') as in_f, open(output_path, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(header)
        for record in parse_star_fusion(in_f):

            left_chr, left_break, left_strand = record.LeftBreakpoint.split(':')

            right_chr, right_break, right_strand = record.RightBreakpoint.split(':')

            fusion = ''.join([record.LeftGene, '--', record.RightGene])
            name = '-'.join([hugo_to_gene_ids[record.LeftGene], hugo_to_gene_ids[record.RightGene]])
            score = 'Junction:%s-Spanning:%s' % (record.JunctionReadCount, record.SpanningFragCount)

            # Add empty sequences in case Trinity doesn't output one
            if len(five_pr_splits[fusion].keys()) == 0:
                five_pr_splits[fusion]['N/A'] = '.'

            if len(three_pr_splits[fusion].keys()) == 0:
                three_pr_splits[fusion]['N/A'] = '.'

            for transcript_id in five_pr_splits[fusion].keys():
                five_prime_seq = five_pr_splits[fusion][transcript_id]
                three_prime_seq = three_pr_splits[fusion][transcript_id]

                writer.writerow([left_chr,
                                 '.',  # Donor start position is not necessary
                                 left_break,
                                 right_chr,
                                 right_break,
                                 '.',  # Acceptor end position is not necessary
                                 name,
                                 score,
                                 left_strand,
                                 right_strand,
                                 five_prime_seq,
                                 three_prime_seq,
                                 record.LeftGene,
                                 record.RightGene])

    bedpe_id = job.fileStore.writeGlobalFile(output_path)
    export_results(job, bedpe_id, 'fusion.bedpe', univ_options, subfolder='mutations/fusions')
    job.fileStore.logToMaster('Reformatted STAR-Fusion output for %s successfully'
                              % univ_options['patient'])
    return bedpe_id
