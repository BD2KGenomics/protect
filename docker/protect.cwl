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
  alignments:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'alignments.tar.gz'
    doc: "Result files Protect CGL pipeline"

  binding_predictions:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'binding_predictions.tar.gz'
    doc: "Result files Protect CGL pipeline"

  expression:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'expression.tar.gz'
    doc: "Result files Protect CGL pipeline"

  haplotyping:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'haplotyping.tar.gz'
    doc: "Result files Protect CGL pipeline"

  mutations:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'mutations.tar.gz'
    doc: "Result files Protect CGL pipeline"

  peptides:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'peptides.tar.gz'
    doc: "Result files Protect CGL pipeline"

  rankboost:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'rankboost.tar.gz'
    doc: "Result files Protect CGL pipeline"

  reports:
    type:
      type: array
      items: File
    outputBinding:
      glob: 'reports.tar.gz'
    doc: "Result files Protect CGL pipeline"

baseCommand: ["--no-clean"]
