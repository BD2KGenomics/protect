

# Overview

ProTECT attempts to predict the neo-antigens in the patient's tumor that are most
likely to stimulate a T-cell based response. Ideally, ProTECT runs using a "trio" of input files,
namely Illumina sequencing data from paired Tumor DNA (WGS or WXS), Normal DNA (WGS or WXS), and
Tumor RNA. In this mode, the samples are all aligned using the same genome reference and annotation,
and somatic mutants are called using a panel of callers. The samples are haplotyped and the
potential neoepitopes are assayed against the alleles in the HLA haplotype for binding. The
neoepitopes are finally subjected to a rank boosting approach to identify the most immunogenic in
the list.

Advanced run modes allow the user to run ProTECT using pre-aligned bam files, or even a
pre-computed vcf of mutations. These modes however put a certain degree of responsibility on the
user regarding the annotations used for mapping and integrity of input files. If any of the input
fastqs is not provided, a corresponding haplotype must be provided in the form of a tar.gz
archive (See below).

NOTE: Input fastq files from multiple lanes and libraries should be combined to yield just 2 files
per sample type (T_DNA, N_DNA, T_RNA) having the standard naming convention -- XYZ1.fq and XYZ2.fq.
Currently ProTECT only supports paired-end runs and only accepts one pair of fastqs per analyte
type.


ProTECT is implemented in the [Toil](https://github.com/BD2KGenomics/toil.git) framework and fully
runs the workflow described in [protect/Flowchart.txt](
https://github.com/BD2KGenomics/protect/blob/master/Flowchart.txt).

**This manual is a quick adaptation for an adaptation of ProTECT to py3** 


# Installation

ProTECT requires Toil and we recommend installing ProTECT and its requirements in a
[virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/).

~ProTECT also requires [s3am](https://github.com/BD2KGenomics/s3am.git) version 2.0.1 to download and
upload files from S3. We recommend installing s3am in its own virtualenv using the directions in
the s3am manual, then putting the s3am binary on your $PATH.  ProTECT will NOT attempt to install
s3am during installation.~ 

currently WIP. for now, **only references to local files will work**. anything that requires access to s3am (s3 buckets) will **fail**. 

~ProTECT uses pkg_resources from setuptools to verify versions of tools during install. As of setuptools
39.0.1, some modules were moved to the packaging module. If your machine has setuptools >=39.0.1, you
will need the packaging module.~

Lastly, ProTECT uses [docker](https://www.docker.com/) to run the various sub-tools in a
reproducible, platform independent manner. ProTECT will NOT attempt to install docker during
installation.

~Method 1 - Using PIP (recommended)~
### Method 2 - Installing from Source

This will install ProTECT in an editable mode.

Obtain the source from Github

    git clone https://www.github.com/Dranion/protect.git

Create and activate a virtualenv in the project folder (Important since the Makefile checks for
this and will fail if it detects that you are not in a virtual environment)

    cd protect
    virtualenv venv
    source venv/bin/activate

Install Toil and pytest

    make prepare

Install packaging (required if setuptools>=39.0.1)

    pip install packaging

Install ProTECT

    make develop

~Method 3 - Using Docker~


# Running ProTECT

Running ProTECT requires you to have activated the virtualenv where it was installed.

Activate the virtualenv with

    source ~/venvs/protect/bin/activate
    OR
    source ~/protect/venv/bin/activate

Running ProTECT is as simple as:

            ProTECT --config /path/to/config.yaml --workDir /path/to/a/workdir JobStore

Where

a) /path/to/config.yaml is the absolute path to the config file for the run.   If a relative path is
provided, ProTECT will attempt to make it an absolute path, however this **might** be an issue on
OSX because folders like `/tmp` are symlinks to `/private/tmp`.

b) /path/to/a/workdir points to a directory to use as a working directory. This directory is where
all the worker sandboxes (temp directories) will be set up during the run. If this option is omitted,
Toil will setup workdir in the system $TMPDIR

c) JobStore is the job store used for the Toil run. The job store can either be a path to a
directory on the system (e.g. /home/probox/toil/jobstore), or an aws bucket as
aws:< region >:< bucket_name > (e.g. aws:us-west-2:ProTEST). Refer to the Toil manual for more
information regarding job stores.

**NOTE:** The workdir **MUST** be on a device that has enough space to run a genomics pipeline.
STAR indexes for a human genome are usually 25Gb and running star usually takes ~60Gb per sample
(for hg19).  You need enough space on this device such that concurrent jobs do not run into
`Unable to free up enough space for caching` errors. If the job store being used is the file job
store (points to a directory on the system), and if the device hosting the job store is the same
one hosting the working directory, you need to ensure the device can handle the concurrent writes
and the persistent job store files as well.

**NOTE:** It is advisable to setup Docker such that it runs off a large volume as well. The
following link contains information on how to setup docker such that it stores images and containers
on another volume.

    http://stackoverflow.com/questions/24309526/how-to-change-the-docker-image-installation-directory

 **TL;DR:** Use large volumes for workdir, job store and Docker.

Alternatively, the Dockerized pipeline can be run from Dockstore (https://dockstore.org/).
The protect.cwl file in the /docker/ directory describes the necessary inputs and outputs and
should be used if running via Dockstore.

# Setting up a config file

A config file pre-filled with references and options for an HG19 run can be generated with

    ProTECT --generate_config

The file `ProTECT_config.yaml` will be created in the current working directory.

The config file contains all the relevant information required for a run. The file is written in the
[YAML format](http://www.yaml.org/start.html) and should be relatively simple to fill in. The file
is split into a number of element groups that describe the important flags passed to different tools
in the pipeline, and the information on the input samples. Elements before a `:` are keys in the
dictionary read into ProTECT and should **NOT** be modified (Barring the patient ID key in the
patients dictionary). Only values to the right of the `:` should be edited.

Every required reference file is provided in the AWS bucket `protect-data` under the folder
`protect/hg19_references` or `protect/hg38_references`. The `README` file in the same location
describes in detail how each file was generated. To use a file located in an s3 bucket, replace
`/path/to` in the following descriptions with `s3://<databucket>/<folder_in_bucket>`.

The following section describes the arguments that can be provided to ProTECT via the config file.

**patients**

The first group of elements is the list of patients to be processed in the run. There is no limit to
the number of patients that can be run in the same workflow and each patient group is marked with a
patient ID, an optional flag describing whether the file is encrypted with sse-c in an aws bucket
(only relevant if files are on s3), and has file combinations  that meet the following conditions:

1. The group contains 3 values which correspond to the forward fastq read file for the Tumor DNA,
Normal DNA, and Tumor RNA samples respectively.  The file extensions for fastqs can be .fq, .fq.gz,
.fastq, or .fastq.gz. If only an `<analyte_type>_fastq_1` flag is passed (without a corresponding
`<analyte_type>_fastq_2` flag), the file names should be symmetrically named `<prefix>1.extn` and
`<prefix>2.extn`.  This is the ideal way to run the pipeline.

2. The group contains a pre-aligned bam (and optionally the corresponding bai file) OR fastq(s) for
each of the analyte types. If any of the input files is not a fastq, the user must also provide a
tar.gz archive with the `hla_haplotype_files` flag that contains 2 files enclosed in a folder
(mhci_alleles.list and mhcii_alleles.list). Each alleles.list file should contain one input MHC
allele (4 digit resolution) per line. NOTE: If RNA-seq fastqs are not provided, the 2 bams, one
with RNA-seq reads mapped to the correct human reference genome build, and one with reads mapped to
a transcriptome inferred from the correct annotation must be provided.

3. The group contains a pre-computed vcf (In lieu of input sequences) and either a trio of fastqs
(for haplotyping the patient), or a pre-computed haplotype bundle and an rna-seq fastq (or pair of
bams; for expression estimation and variant filtering).

In addition to local and AWS S3 hosted files, input bams and vcfs can be pulled form the NCBI GDC
using gdc://<UUID> as the value for the coresponding key. Since GDC bams are indexed, the .bai from
the data bundle will be used.

```
patients:    -> This is the group name. Do not modify this
    PRTCT-01:   -> A string describing a unique Identifier for the patient. This is the only key
                   that can be modified by the user.
        tumor_dna_fastq_1: /path/to/Tum_1.fq
        tumor_dna_fastq_2: /path/to/Tum_2.fq
        tumor_dna_bam: /path/to/Tum.bam
        tumor_dna_bai: /path/to/Tum.bam.bai
        normal_dna_fastq_1: /path/to/Norm_1.fq.gz
        normal_dna_fastq_2: /path/to/Norm_2.fq.gz
        normal_dna_bam: /path/to/Norm.bam
        normal_dna_bai: /path/to/Norm.bam.bai
        tumor_rna_fastq_1: /path/to/Rna_1.fq.gz
        tumor_rna_fastq_2: /path/to/Rna_2.fq.gz
        tumor_rna_bam: /path/to/Rna.bam
        tumor_rna_bai: /path/to/Rna.bam.bai
        tumor_rna_transcriptome_bam: /path/to/RNA_transcriptome.bam
        hla_haplotype_files: /path/to/hla_haplotypes.tar.gz
        mutation_vcf: /path/to/mutations.vcf
        tumor_type: STAD # closest TCGA tumor type
        ssec_encrypted: False   # Are the files for this patient encrypted with ssec? (Optional)
        filter_for_OxoG: False  # Should the samples be processed for filtering OxoG artefacts (optional)
```


Encrypted data will be decrypted using per-file SSE-C keys hashed from a master sse-key
(see `Universal_Options`). The backend download and upload from S3 is done using
[s3am](https://www.github.com/BD2KGenomics/s3am.git).

**Universal_Options**

These describe options that are used universally by most tools/jobs in the workflow.


    Universal_Options:
        dockerhub: aarjunrao                    -> All tools used in the pipeline are dockerized to
                                                   allow for easy reproducibility across operating
                                                   systems. ProTECT was tested using the hub owned
                                                   by `aarjunrao` but will work with any hub given
                                                   it contains the required tools described in
                                                   `required_docker_tools.txt`.
        java_Xmx: 20G                           -> The default Java heap space to be provided to
                                                   tools.  Per-tool heap space can be specified for
                                                   some tools to override this value.
        reference_build: hg19                   -> The reference build used in this run. Can be
                                                   hg19, hg38, GRCh37 or GRCh38.
        sse_key: /path/to/master.key            -> Used to create per-file SSE-C keys for decrypting
                                                   S3 hosted input files.  It is highly recommended
                                                   that the files be uploaded to S3 using s3am using
                                                   the --sse-key-is-master flag. This file should be
                                                   present at the same location on all workers in
                                                   the workflow.
        sse_key_is_master: True                 -> Were the sample files all encrypted with sse_key
                                                   (`False`), or were they encrypted with individual
                                                   per-file keys hashed from the master sse_key
                                                   (`True`)
        gdc_download_token: /path/to/token.txt  -> If any of the input files are being pulled from
                                                   the NCBI GDC, this token is required to access
                                                   the data. This file should be present at the same
                                                   location on all workers in the workflow.
        storage_location: aws:protect-out       -> Should the intermediate files be stored locally
                                                   (`Local`) or in an Amazon S3 bucket using SSE-C
                                                   per-file encryption from the provided master sse
                                                   key (`aws:<bucket_name>`)?
        output_folder: /path/to/out             -> A path to a folder where intermediate, and output
                                                   files (mutation calls, MHC haplotype, etc). This
                                                   folder will be created if it doesn't exist.
                                                   However, files in the folder from a previous run
                                                   will not be overwritten with updated values.
        mail_to: test.email@host.com            -> Optionally send a success email to the specified email
                                                   address. Ensure that results@protect.cgl.genomics.ucsc.edu
                                                   is in your address book.



In reality, for the inexperienced user, these are the only values that need to be modified for a
run.  Experienced users can provide their own references below for better control over the results.
Look at the help for each tool for details. All ProTECT provided indexes were generated using
[Gencode release 19 for hg19](http://www.gencodegenes.org/releases/19.html). Any and all paths can
be substituted with S3 links. Descriptions for creating all files can be found in the file
`S3://cgl-protect-data/hg19_references/README`.

    alignment:
        cutadapt:                                             -> RNA-Seq QC
            a: AGATCGGAAGAG                                       -> Adapter for trimming forward
                                                                     read
            A: AGATCGGAAGAG                                       -> Adapter for trimming reverse
                                                                     read
            version: 1.2
        star:                                                 -> RNA-Seq Alignment
            type: star                                            -> Can be `star` or `starlong`
                                                                     depending on the read size.
                                                                     Refer to the STAR manual for
                                                                     details.
            index_tar: /path/to/star_index.tar.gz                 -> The indexes to use for STAR.
                                                                     Protect provides indexes
                                                                     created using edge sizes 50,
                                                                     75, 100, 150 and 250 in the S3
                                                                     bucket `cgl-protect-data`
                                                                     under the folder
                                                                     `hg19_references`.
            version: 2.4.2a
        bwa:                                                  -> DNA-Seq alignment
            index_tar: /path/to/bwa_index.tar.gz                  -> The indexes to use for bwa.
                                                                     Protect provides indexes in the
                                                                     S3 bucket `cgl-protect-data`
                                                                     under the folder
                                                                     `hg19_references`.
            version: 0.7.9a
        post:                                                 -> Post-alignment processing tools
            samtools:                                           -> Indexing, sam to bam conversion,
                                                                   etc
                version: 1.2
            picard:                                             -> Fixing bam headers
                version: 1.135

    expression_estimation:
        rsem:
            index_tar: /path/to/rsem_index.tar.gz                 -> The indexes to use for rsem.
                                                                     Protect provides indexes in the
                                                                     S3 bucket `cgl-protect-data`
                                                                     under the folder
                                                                     `hg19_references`.
            version: 1.2.0

    mutation_calling:
        indexes:
            chromsosomes: canonical, canonical_chr, chr1, Y       -> A list of chromosomes to process.
                                                                     This options overrides the
                                                                     chromosomes in the fai and dict
                                                                     files. canonical implies
                                                                     1,2,..21,22,X,Y and canonical_chr
                                                                     is the same but with the chr
                                                                     prefix. Cannot have a mix of
                                                                     chr and non-chr. This value is
                                                                     optional.
            genome_fasta: /path/to/hg19.fa.tar.gz                 -> The genome fasta to use for the
                                                                     mutation callers.
            genome_fai: /path/to/hg19.fa.fai.tar.gz               -> The corresponding .fai file for
                                                                     the genome fasta
            genome_dict: /path/to/hg19.dict.tar.gz                -> The corresponding .dict file
                                                                     for the genome fasta
            cosmic_vcf: /path/to/CosmicCodingMuts.vcf.tar.gz      -> The Cosmic Coding vcf
            cosmic_idx: /path/to/CosmicCodingMuts.vcf.idx.tar.gz  -> The corresponding .idx file for
                                                                     the Cosmic vcf
            dbsnp_vcf: /path/to/dbsnp_coding.vcf.gz               -> The dbSNP vcf
            dbsnp_idx: /path/to/dbsnp_coding.vcf.idx.tar.gz       -> The corresponding .idx file for
                                                                     the dbSNP vcf
            dbsnp_tbi : /path/to/dbsnp_coding.vcf.gz.tbi          -> The tabix index for dbsnp.gz
        mutect:
            java_Xmx: 5G                                          -> The heap size to use for MuTect
                                                                     per job (i.e. per chromosome)
            version: 1.1.7
        muse:
            version: 1.0rc_submission_b391201
        radia:                                                -> Radia uses perchrom bed files in
                                                                 folders as references.
            cosmic_beds: /path/to/radia_cosmic.tar.gz
            dbsnp_beds: /path/to/radia_dbsnp.tar.gz
            retrogene_beds: /path/to/radia_retrogenes.tar.gz
            pseudogene_beds: /path/to/radia_pseudogenes.tar.gz
            gencode_beds: /path/to/radia_gencode.tar.gz
            version: bcda721fc1f9c28d8b9224c2f95c440759cd3a03
        somaticsniper:
            version: 1.0.4
            samtools:                                           -> pileup reads
                version: 0.1.8
            bam_readcount:                                      -> obtain readcounts
                version: 0.7.4
        strelka:
            version: 1.0.15
            config_file: /path/to/strelka_config.ini.tar.gz       -> The Strelka config file for a
                                                                     bwa run (modified for a WXS run
                                                                     if necessary)
        star_fusion:
            run: True                                           -> Switch to skip fusion calling
            version: 1.0.0

        fusion_inspector:
            run_trinity: True                                   -> Switch to skip de novo transcript assembly
            version: 1.0.1

    mutation_annotation:
        snpeff:
            index_tar: /path/to/snpeff_index.tar.gz               -> The indexes to use for snpeff.
                                                                     Protect provides indexes in the
                                                                     S3 bucket `cgl-protect-data`
                                                                     under the folder
                                                                     `hg19_references`.
            java_Xmx: 5G                                          -> The heap size to use for SnpEFF
            version: 3.6

    mutation_translation:                                     -> Translate events from genomic to
                                                                 proteomic space
        transgene:
            gencode_peptide_fasta: /path/to/gencode.faa           -> The peptide file for the
                                                                     gencode gtf used in the
                                                                     analysis. (If not gencode, the
                                                                     file must be made to follow the
                                                                     gencode format for fasta record
                                                                     names

            gencode_transcript_fasta : /path/to/gencode_transcripts.faa     -> The transcript file for the gencode gtf.
            gencode_annotation_gtf : /path/to/gencode_annotation.gtf        -> The gencode genome annotation file.
            genome_fasta : /path/to/hg19.faa                                -> The gencode genome fasta file
            version: 2.2.2

    haplotyping:
            phlat:
                index_tar: /path/to/snpeff_index.tar.gz           -> The indexes to use for PHLAT.
                                                                     Protect provides indexes in the
                                                                     S3 bucket `cgl-protect-data`
                                                                     under the folder
                                                                     `hg19_references`.

    mhc_peptide_binding:
        mhci:
            method_file: /path/to/mhci_restrictions.json.tar.gz   -> A json list of allowable MHCs
                                                                     each predictors can handle.
            pred: IEDB_recommended                                -> The IEDB method to use.
            version: 2.13
        mhcii:
            method_file: /path/to/mhcii_restrictions.json.tar.gz  -> A json list of allowable MHCs
                                                                     the predictors can handle.
            pred: IEDB_recommended                                -> The IEDB method to use.
            version: 2.13
        netmhciipan:
            version: 3.1

    prediction_ranking:
        rank_boost:
            mhci_args:                                        -> Weights for the mhci ranking
                npa: 0.0                                          -> Number of peptides in an IAR
                nph: 0.0                                          -> Number of "good" (top 1%)
                                                                     peptides in an IAR
                nMHC: 0.32                                        -> Number of MHCs stimulated by
                                                                     the IAR
                TPM: 0.0                                          -> Expression of the gene
                                                                     containing the IAR
                overlap: 0.68                                     -> Number of 9/10-mer overlap
                                                                     events in the IAR
                tndelta: 0.0                                      -> The number of "good" (>1.5%)
                                                                     altered-self events seen in the
                                                                     IAR
            mhcii_args:                                       -> Weights for the mhci ranking
                npa: 0.2
                nph: 0.2
                nMHC: 0.2
                TPM: 0.2
                tndelta: 0.2
            version: 2.0.1

    reports:
        mhc_pathways_file: /path/to/mhc_pathways.tsv.tar.gz       -> A tsv file containing information
                                                                     about the various genes in the MHC
                                                                     pathway and their mean TCGA and GTEX
                                                                     expressions.
        itx_resistance_file: /path/to/itx_resistance.tsv.tar.gz   -> A tsv file containing information
                                                                     about the mean TCGA and GTEX values
                                                                     for genes involved in immunotherapy
                                                                     resistance pathways.
        immune_resistance_pathways_file: /path/to/immune_resistance_pathways.json.tar.gz
                                                                  -> A json file containing information
                                                                     about immunotherapy resistance
                                                                     pathways.
        car_t_targets_file: /path/to/car_t_targets.tsv.tar.gz     -> A tsv file containiing information
                                                                     about various genes investigated as
                                                                     possible targets for CAR T cell gene
                                                                     therapy, including DOI for scientific
                                                                     literature and clinical trials IDs.




# Default values for config entries (Only in source installations)

If you installed from source, we provide a way to specify default values for all entries in the
config file (except for the `patients` tab). The file `src/protect/pipeline/defaults.yaml` will be
used to fill in default values for config entries if the user does not specify them in the input
config file.


# A note on dockerized tools used in pipeline

ProTECT uses dockerized versions of every tool used during the run to ensure reproducibility of
results. Every docker image required for the run is described in required_docker_tools.txt. Every
required tool is also hosted freely on the dockerhub `aarjunrao` (Hence the default value in the
config). If you wish to use a personal repo instead, ensure that every required version of every
tool is made available on the repo.

ProTECT is registered on Dockstore.org and can be run using Dockstore with the following command:
`
dockstore tool launch --debug --entry quay.io/ucsc_cgl/protect:<version tag>  --json protect.json
`

The example protect.json file provided in this repo describes a run where the input reference files
are on the local file system at `/home/ubuntu/protect-work/protect-reference-files/hg38-reference-files`
and the input fastqs are in `/home/ubuntu/protect-work/samples/`.

# Troubleshooting
## Installation
#### 1. packaging related issues
If you see

    Traceback (most recent call last):
      File "<string>", line 1, in <module>
      File "/mnt/dstew/protect/setup.py", line 6, in <module>
        from packaging.version import LegacyVersion as _LegacyVersion
    ImportError: No module named packaging.version


Then it is mostly due to a version of setuptools >= 39.0.1. The fix for this is to run the
following line before reattempting to install protect.

    pip install packaging

## Errors during Runtime
### How to read and understand an error arising during/at the end of the run
Since ProTECT is written in TOIL, the error printed on a failed job is the logfile for the worker
that failed to execute the python function that was the failed job.

An example error is shown here with line numbers manually added at the beginning for annotation
purposes:

    01: The job seems to have left a log file, indicating failure: 'merge_mhc_peptide_calls' g/f/jobO4yiE4
    02: g/f/jobO4yiE4    ---TOIL WORKER OUTPUT LOG---
    03: g/f/jobO4yiE4    INFO:toil:Running Toil version 3.8.0-4c83830e4f42594d995e01ccc07b47396b88c9e7.
    04: g/f/jobO4yiE4    INFO:toil.fileStore:Starting job ('merge_mhc_peptide_calls' g/f/jobO4yiE4) with ID (972acf421c864831d756dec528bb9cc2a4d3c281).
    05: g/f/jobO4yiE4    INFO:toil.fileStore:LOG-TO-MASTER: Merging MHC calls
    06: g/f/jobO4yiE4    Traceback (most recent call last):
    07: g/f/jobO4yiE4      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/worker.py", line 340, in main
    08: g/f/jobO4yiE4        job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore)
    09: g/f/jobO4yiE4      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/job.py", line 1289, in _runner
    10: g/f/jobO4yiE4        returnValues = self._run(jobGraph, fileStore)
    11: g/f/jobO4yiE4      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/job.py", line 1234, in _run
    12: g/f/jobO4yiE4        return self.run(fileStore)
    13: g/f/jobO4yiE4      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/job.py", line 1406, in run
    14: g/f/jobO4yiE4        rValue = userFunction(*((self,) + tuple(self._args)), **self._kwargs)
    15: g/f/jobO4yiE4      File "/home/ucsc/arjun/tools/protect_toil_clean/local/lib/python3/site-packages/protect/binding_prediction/common.py", line 566, in merge_mhc_peptide_calls
    16: g/f/jobO4yiE4        raise RuntimeError('No peptides available for ranking')
    17: g/f/jobO4yiE4    RuntimeError: No peptides available for ranking
    18: g/f/jobO4yiE4    ERROR:toil.worker:Exiting the worker because of a failed job on host sjcb10st7
    19: g/f/jobO4yiE4    WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'merge_mhc_peptide_calls' g/f/jobO4yiE4 with ID g/f/jobO4yiE4 to 0
    20: g/f/jobO4yiE4    WARNING:toil.jobGraph:We have increased the default memory of the failed job 'merge_mhc_peptide_calls' g/f/jobO4yiE4 to 2147483648 bytes

* Line 1 tells you which function failed.
* Line 6 is where the actual error starts (the line starting with `Traceback`)
* Line 16 tells you the tool failed with a RuntimeError and contains the actual error message. In
this case, the tool found no peptides that were available for ranking.
* The line numbers will vary with different tools and dependin gon how jobs were scheduled, however
a rule of thumb is to look for the string `RuntimeError` (the most common error thrown by ProTECT)
or just `Error` and then read the error message.

### A note on non-deterministic errors with RADIA
It has been noticed that RADIA sometimes fails non-deterministically on one (or rarely more)
chromosome in a sample during a ProTECT run. The error is highly non-informative since we currently
do not store logs from tools (see BD2KGenomics/protect#275). The error looks similar to this:

    Z/O/job1uH92D    ---TOIL WORKER OUTPUT LOG---
    Z/O/job1uH92D    INFO:toil:Running Toil version 3.8.0-4c83830e4f42594d995e01ccc07b47396b88c9e7.
    Z/O/job1uH92D    INFO:toil.fileStore:Starting job ('run_filter_radia' Z/O/job1uH92D) with ID (0e8bb5ea5d785324b5abce7b384bf5b99140aa6c).
    Z/O/job1uH92D    WARNING:toil.fileStore:LOG-TO-MASTER: Job used more disk than requested. Please reconsider modifying the user script to avoid the chance  of failure due to incorrec
    tly requested resources. Job 'run_filter_radia' Z/O/job1uH92D used 104.10% (34.3 GB [36816449536B] used, 32.9 GB [35367908263B] requested) at the end of its run.
    Z/O/job1uH92D    Traceback (most recent call last):
    Z/O/job1uH92D      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/worker.py", line 340, in main
    Z/O/job1uH92D        job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore)
    Z/O/job1uH92D      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/job.py", line 1289, in _runner
    Z/O/job1uH92D        returnValues = self._run(jobGraph, fileStore)
    Z/O/job1uH92D      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/job.py", line 1234, in _run
    Z/O/job1uH92D        return self.run(fileStore)
    Z/O/job1uH92D      File "/home/ucsc/arjun/tools/dev/toil_clean/src/toil/job.py", line 1406, in run
    Z/O/job1uH92D        rValue = userFunction(*((self,) + tuple(self._args)), **self._kwargs)
    Z/O/job1uH92D      File "/home/ucsc/arjun/tools/protect_toil_clean/local/lib/python3/site-packages/protect/mutation_calling/radia.py", line 238, in run_filter_radia
    Z/O/job1uH92D        tool_version=radia_options['version'])
    Z/O/job1uH92D      File "/home/ucsc/arjun/tools/protect_toil_clean/local/lib/python3/site-packages/protect/common.py", line 138, in docker_call
    Z/O/job1uH92D        'for command \"%s\"' % ' '.join(call),)
    Z/O/job1uH92D    RuntimeError: docker command returned a non-zero exit status (1)for command "docker run --rm=true -v /scratch/bio/ucsc/toil-681c097c-61da-4687-b734-c5051f0aa19f/tmped2fnu/f041f939-5c0d-40be-a884-68635e929d09:/data --log-driver=none aarjunrao/filterradia:bcda721fc1f9c28d8b9224c2f95c440759cd3a03 TCGA-CH-5788 17 /data/radia.vcf /data /home/radia/scripts -d /data/radia_dbsnp -r /data/radia_retrogenes -p /data/radia_pseudogenes -c /data/radia_cosmic -t /data/radia_gencode --noSnpEff --noBlacklist --noTargets --noRnaBlacklist -f /data/hg38.fa --log=INFO -g /data/radia_filtered_chr17_radia.log"
    Z/O/job1uH92D    ERROR:toil.worker:Exiting the worker because of a failed job on host sjcb10st1
    Z/O/job1uH92D    WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'run_filter_radia' Z/O/job1uH92D with ID Z/O/job1uH92D to 0

The fix for this is to just rerun protect with the same jobstore (if you haven't specified
`--clean always`) with the `--restart` flag. I.e. run exactly the same command as the original run
with `--restart` appended at the end of the command.

For example, if your command was:

    ProTECT --config_file test.yaml --workDir /tmp/working/protect_test /tmp/jobstores/protect_test

then your restart command becomes:

    ProTECT --config_file test.yaml --workDir /tmp/working/protect_test /tmp/jobstores/protect_test \
            --restart


