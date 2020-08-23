#!/usr/bin/env python2.7
# Copyright 2016 UCSC Computational Genomics Lab
# Original contributor: Arjun Arkal Rao
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

from pipelineWrapper import PipelineWrapperBuilder
import logging
import os
import yaml


logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


desc = """UCSC Precision Immuno pipeline"""

config = ("""patients:
    {sample_name}:
        tumor_dna_fastq_1 : {tumor_dna}
        tumor_dna_fastq_2 : {tumor_dna2}
        normal_dna_fastq_1 : {normal_dna}
        tumor_type: {tumor_type}
        ssec_encrypted: {ssec_encrypted}
        filter_for_OxoG: {filter_for_OxoG}
        normal_dna_fastq_2 : {normal_dna2}
        tumor_rna_fastq_1 : {tumor_rna}
        tumor_rna_fastq_2 : {tumor_rna2}

Universal_Options:
    dockerhub: {dockerhub}
    java_Xmx: {javaxmx}
    reference_build: {reference_build}
    sse_key: {sse_key}
    sse_key_is_master: {sse_key_is_master}
    gdc_download_token: {gdc_download_token}
    storage_location: Local
    output_folder: {work_mount}/output
    mail_to: {mail_to}

# These options are for each specific tool. You probably don't need to change any of this!
alignment:
    cutadapt:
        version : {cutadapt_ver}
        a : {cutadapt_a}
        A : {cutadapt_A}

    star:
        version: {star_ver}
        type : {star_type}
        index : {star_index}
    bwa:
        version: {bwa_ver}
        index : {bwa_index}
    post:
        samtools:
            version: {samtools_alignment_ver}
        picard:
            version: {picard_ver}

expression_estimation:
    rsem:
        version: {rsem_ver}
        index : {rsem_index}

mutation_calling:
    indexes:
        chromosomes: {chromosomes}
        genome_fasta : {genome_fasta}
        genome_fai : {genome_fai}
        genome_dict : {genome_dict}
        cosmic_vcf : {cosmic_vcf}
        cosmic_idx : {cosmic_idx}
        dbsnp_vcf : {dbsnp_vcf}
        dbsnp_idx : {dbsnp_idx}
        dbsnp_tbi : {dbsnp_tbi}
    mutect:
        version: {mutect_ver}
        java_Xmx : {mutect_javaxmx}
    muse:
        version: {muse_ver}
    radia:
        version: {radia_ver}
        cosmic_beds: {cosmic_beds}
        dbsnp_beds: {dbsnp_beds}
        retrogene_beds: {retrogene_beds}
        pseudogene_beds: {pseudogene_beds}
        gencode_beds: {gencode_beds}
    somaticsniper:
        version: {somaticsniper_ver}
        samtools:
            version: {samtools_somaticsniper_ver}
        bam_readcount:
            version: {bamreadcount_ver}
    star_fusion:
        run: {starfusion}
        version: {star_fusion_ver}
    fusion_inspector:
        run_trinity: {run_trinity}
        version: {fusioninspector_ver}
    strelka: 
        version: {strelka_ver}
        config_file: {strelka_config}

mutation_annotation:
    snpeff:
        version: {snpeff_ver}
        index : {snpeff}
        java_Xmx : {spnff_javaxmx}

mutation_translation:
    transgene:
        version: {transgene_ver}
        gencode_peptide_fasta : {transgene_peptide_fasta}
        gencode_transcript_fasta : {transgene_transcript_fasta}
        gencode_annotation_gtf : {transgene_annotation_gtf}
        genome_fasta : {transgene_genome}

haplotyping:
    phlat:
        version: {phlat_ver}
        index : {phlat}

mhc_peptide_binding:
    mhci:
        version: {mhci_ver}
        method_file : {mhci}
        pred : {mhci_pred}

    mhcii:
        version: {mhcii_ver}
        method_file : {mhcii}
        pred : {mhcii_pred}
    netmhciipan:
        version: {netmhciipan_ver}

prediction_ranking:
    rankboost:
        version: {rankboost_ver}
        mhci_args: 
            npa: {mhci_npa}
            nph: {mhci_nph}
            nMHC: {mhci_nMHC}
            TPM: {mhci_TPM}
            overlap: {mhci_overlap}
            tndelta: {mhci_tndelta}
        mhcii_args:
            npa: {mhcii_npa}
            nph: {mhcii_nph}
            nMHC: {mhcii_nMHC}
            TPM: {mhcii_TPM}
            tndelta: {mhcii_tndelta}

reports:
    mhc_pathways_file: {mhc_pathways_file}
    itx_resistance_file: {itx_resistance_file}
    immune_resistance_pathways_file: {immune_resistance_pathways_file}
    car_t_targets_file: {car_t_targets_file}""")


def str2bool(v):
    """
    Necessary due to how argsparse works; see 
    https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    """
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == '__main__':
    import protect
    with open(os.path.join(os.path.dirname(protect.__file__), "pipeline",
                           "defaults.yaml")) as def_file:
        defaults = yaml.load(def_file)
    wrapper = PipelineWrapperBuilder('ProTECT', desc, config)
    parser = wrapper.get_args()

    parser.add_argument('--sample-name', type=str, required=True,
                        help="Name for the sample.")
    parser.add_argument('--tumor-dna',
                        type=str, required=True,
                        help='Path for the tumor fastq.')
    parser.add_argument('--normal-dna', type=str, required=True,
                        help='Path for the normal fastq.')
    parser.add_argument('--tumor-rna', type=str, required=True,
                        help='Path for the tumor RNA fastq.')
    parser.add_argument('--tumor-type', type=str, required=True,
                        help='Tumor type (such as STAD).')
    parser.add_argument('--tumor-dna2', type=str, required=True,
                        help='Path for the tumor fastq pair.')
    parser.add_argument('--normal-dna2', type=str, required=True,
                        help='Path for the normal fastq.')
    parser.add_argument('--tumor-rna2', type=str, required=True,
                        help='Path for the tumor RNA fastq.')

    parser.add_argument('--reference-build', type=str,
                        choices=['hg19', 'hg38'], default='hg19',
                        help='Reference build. Can be hg19 or hg38.')
    parser.add_argument('--ssec-encrypted', type=str2bool, default='False')
    parser.add_argument('--filter-for-OxoG', type=str2bool, default='False')

    parser.add_argument('--cutadapt-a', type=str, default="AGATCGGAAGAG")
    parser.add_argument('--cutadapt-A', type=str, default="AGATCGGAAGAG")
    parser.add_argument('--cutadapt_ver', type=str,
                        default=defaults["alignment"]["cutadapt"]["version"],
                        help='Version of cutadapt.')

    parser.add_argument('--star-type', type=str, choices=['star','starlong'], default='star',
                        help='Use starlong if your reads are > 150bp')
    parser.add_argument('--star-index', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "star_with_fusion_100bp_readlen_indexes.tar.gz",
                        help='Index for star.')
    parser.add_argument('--star-ver', type=str,
                        default=defaults["alignment"]["star"]["version"],
                        help='Version of star.')

    parser.add_argument('--bwa-index', type=str,
                        default="S3://protect-data/<reference_build>_references/bwa_index.tar.gz",
                        help='Path for bwa index.')
    parser.add_argument('--bwa-ver', type=str,
                        default=defaults["alignment"]["bwa"]["version"],
                        help='Version of bwa.')

    parser.add_argument('--samtools_alignment_ver', type=str,
                        default=defaults["alignment"]["post"]["samtools"]["version"],
                        help='Version of samtools for alignment.')

    parser.add_argument('--picard-ver', type=str,
                        default=defaults["alignment"]["post"]["picard"]["version"],
                        help='Version of picard.')

    parser.add_argument('--rsem-index', type=str,
                        default="S3://protect-data/<reference_build>_references/rsem_index.tar.gz",
                        help='Path for rsem index.')
    parser.add_argument('--rsem-ver', type=str,
                        default=defaults["expression_estimation"]["rsem"]["version"],
                        help='Version of rsem.')

    parser.add_argument('--mutect-ver', type=str,
                        default=defaults["mutation_calling"]["mutect"]["version"],
                        help='Version of mutect.')
    parser.add_argument('--mutect-javaxmx', type=str, default='2G')

    parser.add_argument('--muse-ver', type=str,
                        default=defaults["mutation_calling"]["muse"]["version"],
                        help='Version of muse.')

    parser.add_argument('--radia-ver', type=str,
                        default=defaults["mutation_calling"]["radia"]["version"],
                        help='Version of radia.')
    parser.add_argument('--cosmic-beds', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "radia_cosmic.tar.gz",
                        help='Cosmic bed file for use by Radia.')
    parser.add_argument('--dbsnp-beds', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "radia_dbsnp.tar.gz",
                        help='dbsnp bed file for use by Radia.')
    parser.add_argument('--retrogene-beds', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "radia_retrogenes.tar.gz",
                        help='Retrogene bed file for use by Radia.')
    parser.add_argument('--pseudogene-beds', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "radia_pseudogenes.tar.gz",
                        help='Psuedogene bed file for use by Radia.')
    parser.add_argument('--gencode-beds', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "radia_gencode.tar.gz",
                        help='Gencode bed file for use by Radia.')

    parser.add_argument('--somaticsniper-ver', type=str,
                        default=defaults["mutation_calling"]["somaticsniper"]["version"],
                        help='Version of somatic sniper.')

    parser.add_argument('--samtools_somaticsniper-ver', type=str,
                        default=defaults["mutation_calling"]["somaticsniper"]["samtools"]["version"],
                        help='Version of samtools for somatic sniper')

    parser.add_argument('--bamreadcount-ver', type=str,
                        default=defaults["mutation_calling"]["somaticsniper"]["bam_readcount"]["version"],
                        help='Version of bam_readcount.')

    parser.add_argument('--strelka-ver', type=str,
                        default=defaults["mutation_calling"]["strelka"]["version"],
                        help='Version of strelka.')
    parser.add_argument('--strelka-config', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "strelka_bwa_WXS_config.ini.tar.gz",
                        help='Path to config for strelka.')

    parser.add_argument('--starfusion', type=str2bool,
                        default=defaults["mutation_calling"]["star_fusion"]["run"],
                        help='Set to false to skip fusion calling.')
    parser.add_argument('--star-fusion-ver', type=str,
                        default=defaults["mutation_calling"]["star_fusion"]["version"],
                        help='Version of star fusion.')

    parser.add_argument('--run-trinity', type=str2bool,
                        default=defaults["mutation_calling"]["fusion_inspector"]["run_trinity"],
                        help='Set to false to skip de novo transcript assembly.')
    parser.add_argument('--fusioninspector-ver', type=str,
                        default=defaults["mutation_calling"]["fusion_inspector"]["version"])

    parser.add_argument('--snpeff-ver', type=str,
                        default=defaults["mutation_annotation"]["snpeff"]["version"],
                        help='Version of snpeff')
    parser.add_argument('--snpeff', type=str,
                        default="S3://protect-data/<reference_build>_references/snpeff_index.tar.gz",
                        help='Path to indexes for snpeff.')
    parser.add_argument('--spnff-javaxmx', type=str, default='20G')

    parser.add_argument('--transgene-peptide-fasta', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "gencode.<reference_gencode>.pc_translations_NOPARY.fa.tar.gz",
                        help='Path to gencode peptide fasta for transgene.')
    parser.add_argument('--transgene-transcript-fasta', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "gencode.<reference_gencode>.pc_transcripts.fa.tar.gz",
                        help='Path to gencode transcript fasta.')
    parser.add_argument('--transgene-annotation-gtf', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "gencode.<reference_gencode>.annotation_NOPARY.gtf.tar.gz",
                        help='Path to gencode annotation gtf.')
    parser.add_argument('--transgene-genome', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "<reference_build>.fa.tar.gz",
                        help='Path to genome fasta.')
    parser.add_argument('--transgene-ver', type=str,
                        default=defaults["mutation_translation"]["transgene"]["version"],
                        help='Version of transgene.')

    parser.add_argument('--phlat-ver', type=str,
                        default=defaults["haplotyping"]["phlat"]["version"],
                        help='Version of phlat.')
    parser.add_argument('--phlat', type=str,
                        default="S3://protect-data/<reference_build>_references/phlat_index.tar.gz",
                        help='Path to config for phlat.')

    parser.add_argument('--mhci-ver', type=str,
                        default=defaults["mhc_peptide_binding"]["mhci"]["version"],
                        help='Version of mhci.')
    parser.add_argument('--mhci', type=str,
                        default="S3://protect-data/<reference_build>_references"
                                "/mhci_restrictions.json.tar.gz",
                        help='Path to config for mhci.')
    parser.add_argument('--mhci-pred', type=str, default='IEDB_recommended')
    parser.add_argument('--mhcii-pred', type=str, default='IEDB_recommended')
    parser.add_argument('--mhci-npa', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhci_args']['npa'])
    parser.add_argument('--mhci-nph', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhci_args']['nph'])
    parser.add_argument('--mhci-nMHC', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhci_args']['nMHC'])
    parser.add_argument('--mhci-TPM', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhci_args']['TPM'])
    parser.add_argument('--mhci-overlap', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhci_args']['overlap'])
    parser.add_argument('--mhci-tndelta', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhci_args']['tndelta'])

    parser.add_argument('--mhcii-ver', type=str,
                        default=defaults["mhc_peptide_binding"]["mhcii"]["version"],
                        help='Version of mhcii.')
    parser.add_argument('--mhcii', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "mhcii_restrictions.json.tar.gz",
                        help='Path to config for mhcii.')
    parser.add_argument('--mhcii-npa', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhcii_args']['npa'])
    parser.add_argument('--mhcii-nph', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhcii_args']['nph'])
    parser.add_argument('--mhcii-nMHC', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhcii_args']['nMHC'])
    parser.add_argument('--mhcii-TPM', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhcii_args']['TPM'])
    parser.add_argument('--mhcii-tndelta', type=str,
                        default=defaults['prediction_ranking']['rankboost']['mhcii_args']['tndelta'])

    parser.add_argument('--netmhciipan-ver', type=str,
                        default=defaults["mhc_peptide_binding"]["netmhciipan"]["version"],
                        help='Version of netmhciipain.')

    parser.add_argument('--rankboost-ver', type=str,
                        default=defaults["prediction_ranking"]["rankboost"]["version"],
                        help='Version of rankboost.')
    parser.add_argument('--chromosomes', type=str, nargs=2, default="")
    parser.add_argument('--genome-fasta', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "<reference_build>.fa.tar.gz",
                        help='Genome fasta to be used by the mutation callers.')
    parser.add_argument('--genome-fai', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "<reference_build>.fa.fai.tar.gz",
                        help='Corresponding fai file for the genome fasta.')
    parser.add_argument('--genome-dict', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "<reference_build>.dict.tar.gz",
                        help='Corresponding dict file for the genome fasta.')

    parser.add_argument('--cosmic-vcf', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "CosmicCodingMuts.vcf.tar.gz",
                        help='vcf for cosmic coding.')
    parser.add_argument('--cosmic-idx', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "CosmicCodingMuts.vcf.idx.tar.gz",
                        help='Corresponding idx for the cosmic coding vcf.')

    parser.add_argument('--dbsnp-vcf', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "dbsnp_coding.vcf.gz",
                        help='vcf for dbsnp.')
    parser.add_argument('--dbsnp-idx', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "dbsnp_coding.vcf.idx.tar.gz",
                        help='Corresponding idx for the dbsnp vcf.')
    parser.add_argument('--dbsnp-tbi', type=str,
                        default="S3://protect-data/<reference_build>_references/"
                                "dbsnp_coding.vcf.gz.tbi",
                        help='Tabix index for dbsnp.gz.')

    parser.add_argument('--mhc-pathways-file', type=str,
                        default="S3://cgl-pipeline-inputs/protect/ci_references/"
                                "mhc_pathways.tsv.tar.gz",
                        help='JSON file containing the various genes in the MHC pathway'
                             'and their mean TPM expressions across samples in a background set.')
    parser.add_argument('--itx-resistance-file', type=str,
                        default="S3://cgl-pipeline-inputs/protect/ci_references/"
                                "itx_resistance.tsv.tar.gz",
                        help='')
    parser.add_argument('--immune-resistance-pathways-file', type=str,
                        default="S3://cgl-pipeline-inputs/protect/ci_references/"
                                "immune_resistance_pathways.json.tar.gz",
                        help='')
    parser.add_argument('--car-t-targets-file', type=str,
                        default="S3://cgl-pipeline-inputs/protect/ci_references/"
                                "car_t_targets.tsv.tar.gzz",
                        help='')

    parser.add_argument('--dockerhub', type=str, default='aarjunrao')
    parser.add_argument('--javaxmx', default='20G', type=str)
    parser.add_argument('--sse-key', type=str, default='',
                        help='Path to the desired SSE-key, if any.')
    parser.add_argument('--sse-key-is-master', type=str2bool, default='False',
                        help='Indicates if the passed sse-key is the master key.')
    parser.add_argument('--gdc-download-token', type=str,
                        help='A download token used to download files from the GDC')
    parser.add_argument('--mail-to', type=str,
                        help='Address to send an email to upon successful completion.')

    parser.add_argument('--work-mount', required=True,
                        help='Mount where intermediate files should be written. This directory '
                             'should be mirror mounted into the container.')

    args = parser.parse_args()
    args.chromosomes = ', '.join(args.chromosomes)
    for key in args.__dict__:
        try:
            args.__dict__[key] = args.__dict__[key].replace('<reference_build>',
                                                            args.reference_build)
        except AttributeError:
            pass
        try:
            args.__dict__[key] = args.__dict__[key].replace('<reference_gencode>',
                                                            'v19' if args.reference_build == 'hg19'
                                                            else 'v25')
        except AttributeError:
            pass
    command = []
    wrapper.run(args, command)
