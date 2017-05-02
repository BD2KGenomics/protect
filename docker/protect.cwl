#!/usr/bin/env cwl-runner

class: CommandLineTool
id: ProTECT-CGL
label: ProTECT CGL Pipeline
cwlVersion: v1.0

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

doc: |
    This repo contains the Python libraries for the Precision Immunology Pipeline developed at UCSC.

    src/protect/pipeline/ProTECT.py             - The python script for running the pipeline.
    src/protect/pipeline/input_parameters.yaml  - The config file for the run that contains all the
                                                  necessary parameters for the run.
    Flowchart.txt                               - A (super cool) flowchart describing the flow of
                                                  the pipeline.


    ProTECT uses sequencing information from a patient to predict the neo-epitopes produced in their
    tumors that can be used in T-cell based, or peptide vaccine based therapies.

    All docker images used in this pipeline are available at

                            https://hub.docker.com/u/aarjunrao/


    To learn how the pipeline can be run on a sample, head over to the [ProTECT Manual](
    https://github.com/BD2KGenomics/protect/blob/master/MANUAL.md)

    ProTECT is currently in its infancy and is under continuous development.  We would appreciate users sharing the level 3 data produced by ProTECT with us such that we can better train our predictive models.


dct:creator:
  '@id': http://orcid.org/0000-0002-7681-6415
  foaf:name: Brian O'Connor
  foaf:mbox: briandoconnor@gmail.com

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/ucsc_cgl/protect:2.3.0--1.12.3"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 64000
    outdirMin: 21474836480
    description: "The process requires at least 16G of RAM and we recommend 21.5GB of storage."

inputs:

  sample-name:
    type: String
    doc: "Name of sample."
    inputBinding:
      prefix: --sample-name

  tumor-dna:
    type: File
    doc: "Tumor DNA fastq"
    inputBinding:
      prefix: --tumor-dna

  tumor-rna:
    type: File
    doc: "Tumor RNA fastq"
    inputBinding:
      prefix: --tumor-rna

  normal-dna:
    type: File
    doc: "Normal DNA fastq"
    inputBinding:
      prefix: --normal-dna

  tumor-dna2:
    type: File
    doc: "Second Tumor DNA fastq"
    inputBinding:
      prefix: --tumor-dna2

  tumor-rna2:
    type: File
    doc: "Second Tumor RNA fastq"
    inputBinding:
      prefix: --tumor-rna2

  normal-dna2:
    type: File
    doc: "Second Normal DNA fastq"
    inputBinding:
      prefix: --normal-dna2

  reference_build:
    type: string
    doc: "Reference build (hg19 or hg38)"
    inputBinding:
      prefix: --reference-build

  cutadapt_ver:
    inputBinding:
      --cutadapt_ver
    type: string
  star_ver:
    inputBinding:
      --star-ver
    type: string
  bwa_ver:
    inputBinding:
      --bwa-ver
    type: string
  samtools_alignment_ver:
    inputBinding:
      --samtools_alignment_ver
    type: string
    doc: "Samtools for alignment"
  picard_ver:
    inputBinding:
      --picard-ver
    type: string
  rsem_ver:
    inputBinding:
      --rsem-ver
    type: string
  mutect_ver:
    inputBinding:
      --mutect-ver
    type: string
  muse_ver:
    inputBinding:
      --muse-ver
    type: string
  radia_ver:
    inputBinding:
      --radia-ver
    type: string
  somaticsniper_ver:
    inputBinding:
      --somaticsniper-ver
    type: string
  samtools_somaticsniper_ver:
    inputBinding:
      --samtools_somaticsniper-ver
    type: string
    doc: "Samtools for somatic sniper"
  bamreadcount_ver:
    inputBinding:
      --bamreadcount-ver
    type: string
  strelka_ver:
    inputBinding:
      --strelka-ver
    type: string
  snpeff_ver:
    inputBinding:
      --snpeff-ver
    type: string
  transgene_ver:
    inputBinding:
      --transgene-ver
    type: string
  phlat_ver:
    inputBinding:
      --phlat-ver
    type: string
  mhci_ver:
    inputBinding:
      --mhci-ver
    type: string
  mhcii_ver:
    inputBinding:
      --mhcii-ver
    type: string
  netmhciipan_ver:
    inputBinding:
      --netmhciipan-ver
    type: string
  rankboost_ver:
    inputBinding:
      --rankboost-ver
    type: string


  star_type:
    type: string
    doc: "Star type. Use starlong if reads >150bp."
    inputBinding:
      prefix: --star-type

  star:
    type: File?
    doc: "Absolute path to STAR index tarball."
    inputBinding:
      prefix: --star-index

  bwa:
    type: File?
    doc: "Absolute path to bwa index tarball."
    inputBinding:
      prefix: --bwa-index

  rsem:
    type: File?
    doc: "Absolute path to rsem reference tarball."
    inputBinding:
      prefix: --rsem-index

  genome_fasta:
    type: File?
    doc: "mut_callers.genome_fasta"
    inputBinding:
      prefix: --genome-fasta

  genome_fai:
    type: File?
    doc: "mut_callers.genome_fai"
    inputBinding:
      prefix: --genome-fai

  genome_dict:
    type: File?
    doc: "mut_callers.genome_dict"
    inputBinding:
      prefix: --genome-dict

  cosmic_vcf:
    type: File?
    doc: "mut_callers.cosmic_vcf"
    inputBinding:
      prefix: --cosmic-vcf

  cosmic_idx:
    type: File?
    doc: "mut_callers.cosmic_idx"
    inputBinding:
      prefix: --cosmic-idx

  dbsnp_vcf:
    type: File?
    doc: "mut_callers.dbsnp_vcf"
    inputBinding:
      prefix: --dbsnp-vcf

  dbsnp_idx:
    type: File?
    doc: "mut_callers.dbsnp_idx"
    inputBinding:
      prefix: --dbsnp-idx

  dbsnp_tbi:
    type: File?
    doc: "mut_callers.dbsnp_tbi"
    inputBinding:
      prefix: --dbsnp-tbi

  dbsnp_beds:
    type: File?
    doc: "beds"
    inputBinding:
      prefix: --dbsnp-beds

  cosmic_beds:
    type: File?
    doc: "Cosmic beds"
    inputBinding:
      prefix: --cosmic_beds

  retrogene_beds:
    type: File?
    doc: "Retrogene beds"
    inputBinding:
      prefix: --retrogene-beds

  psuedogene_beds:
    type: File?
    doc: "Psuedogene beds"
    inputBinding:
      prefix: --psuedogene_beds

  gencode_beds:
    type: File?
    doc: "Gencode beds
    inputBinding:
      prefix: --gencode_beds

  strelka_config:
    type: File?
    doc: "mut_callers.strelka_config"
    inputBinding:
      prefix: --strelka-config

  snpeff:
    type: File?
    doc: "snpeff"
    inputBinding:
      prefix: --snpeff

  transgene:
    type: File?
    doc: "transgene"
    inputBinding:
      prefix: --transgene

  phlat:
    type: File?
    doc: "phlat"
    inputBinding:
      prefix: --phlat

  mhci:
    type: File?
    doc: "mhci"
    inputBinding:
      prefix: --mhci

  mhcii:
    type: File?
    doc: "mhcii"
    inputBinding:
      prefix: --mhcii

  mhc_pathway_assessment:
    type: File?
    doc: "mhc_pathway_assessment"
    inputBinding:
      prefix: --mhc-pathway-assessment

  work_mount:
    type: string
    doc: "Path of the working directory to be mounted into the container"
    inputBinding:
      prefix: --work-mount

  sse_key:
    type: string?
    doc: 'Path to encryption key'
    prefix: --sse-key

  sse_key_is_master:
    type: string?
    doc: 'Required if on master'
    prefix: --sse-key-is-master

outputs:

  binding_predictions:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'binding_predictions.tar'
    doc: "Result files from ProTECT"

  expression:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'expression.tar'
    doc: "Result files from ProTECT"

  haplotyping:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'haplotyping.tar'
    doc: "Result files from ProTECT"

  merged_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'merged_perchrom.tar'
    doc: "Result files from ProTECT"

  muse_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'muse_perchrom.tar'
    doc: "Result files from ProTECT"

  mutect_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'mutect_perchrom.tar'
    doc: "Result files from ProTECT"

  peptides:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'peptides.tar'
    doc: "Result files from ProTECT"

  radia_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'radia_perchrom.tar'
    doc: "Result files from ProTECT"

  somaticsniper_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'somaticsniper_perchrom.tar'
    doc: "Result files from ProTECT"

  strelka_snv_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'strelka_snv_perchrom.tar'
    doc: "Result files from ProTECT"

  strelka_indel_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'strelka_indel_perchrom.tar'
    doc: "Result files from ProTECT"

  rankboost:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'rankboost.tar'
    doc: "Result files from ProTECT"

  reports:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'reports.tar'
    doc: "Result files from ProTECT"

  normal_alignment:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'normal_dna_fix_pg_sorted.bam'
    doc: "Result files from ProTECT"

  normal_index:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'normal_dna_fix_pg_sorted.bam.bai'
    doc: "Result files from ProTECT"

  tumor_alignment:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'tumor_dna_fix_pg_sorted.bam'
    doc: "Result files from ProTECT"

  tumor_index:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'tumor_dna_fix_pg_sorted.bam.bai'
    doc: "Result files from ProTECT"

  rna_alignment:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'rna_fix_pg_sorted.bam'
    doc: "Result files from ProTECT"

  rna_index:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'rna_fix_pg_sorted.bam.bai'
    doc: "Result files from ProTECT"

  rna_transcriptome_alignment:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'rna_transcriptome.bam'
    doc: "Result files from ProTECT"

  all_merged:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'all_merged.vcf'
    doc: "Result files from ProTECT"

  mhci_merged:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'mhci_merged_files_concise_results.tsv'
    doc: "Result files from ProTECT"

  mhcii_merged:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'mhcii_merged_files_concise_results.tsv'
    doc: "Result files from ProTECT"

  all_snpeffed:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'all_snpeffed.vcf'
    doc: "Result files from ProTECT"

  all_transgened:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'all_transgened.vcf'
    doc: "Result files from ProTECT"

baseCommand: ["--no-clean"]
