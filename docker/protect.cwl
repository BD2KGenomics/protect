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
  tumor-dna:
    type: File[]
    doc: "Tumor DNA fastq"
    inputBinding:
      prefix: --tumor-dna

  tumor-rna:
    type: File[]
    doc: "Tumor RNA fastq"
    inputBinding:
      prefix: --tumor-rna

  normal-dna:
    type: File[]
    doc: "Normal DNA fastq"
    inputBinding:
      prefix: --normal-dna

  tumor-dna2:
    type: File[]
    doc: "Second Tumor DNA fastq"
    inputBinding:
      prefix: --tumor-dna2

  tumor-rna2:
    type: File[]
    doc: "Second Tumor RNA fastq"
    inputBinding:
      prefix: --tumor-rna2

  normal-dna2:
    type: File[]
    doc: "Second Normal DNA fastq"
    inputBinding:
      prefix: --normal-dna2

  autoscale:
    type: boolean
    doc: "Indicates whether to use Toil's autoscaling capabilities"
    inputBinding:
      prefix: --autoscale

  provisioner:
    type:string
    doc: "Sets up where autoscaling is going to happen"
    inputBinding:
        prefix: --provisioner

  nodeType:
    type:string
    doc: "Sets the type of node used to help with the clustering"
    inputBinding:
        prefix: --nodeType

  storage-location:
    type:string
    doc:"Sets the destination where output should be held"
    inputBinding:
        prefix: --storage_location

  output-folder:
    type:string
    doc:"Sets the source from where output is coming."
    inputBinding:
        prefix: --output_folder

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

  strelka_perchrom:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'strelka_perchrom.tar'
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
