

#Overview

ProTECT is implemented in the [Toil](https://github.com/BD2KGenomics/toil.git) framework and fully
runs the workflow described in [protect/Flowchart.txt](
https://github.com/BD2KGenomics/protect/blob/master/Flowchart.txt).

ProTECT accepts Illumina sequencing data from paired Tumor DNA (WGS or WXS), Normal DNA (WGS or
WXS), and Tumor RNA, and attempts to predict the neo-antigens in the patient's tumor that are most
likely to stimulate a T-cell based response. Input fastq files from multiple lanes and libraries
should be combined to yield just 2 files per sample type (T_DNA, N_DNA, T_RNA) having the standard
naming convention -- XYZ1.fq and XYZ2.fq . Currently ProTECT only supports paired-end runs.

#Installation

ProTECT requires Toil and we recommend installing ProTECT and its requirements in a
[virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/).

ProTECT also requires [s3am](https://github.com/BD2KGenomics/s3am.git) version 2.0 to download and
upload files from S3. We recommend installing s3am in its own virtualenv using the directions in
the s3am manual, then putting the s3am binary on your $PATH.

###Method 1 - Using PIP (recommended)

First create a virtualenv at your desired location (Here we create it in the folder ~/venvs)

    virtualenv ~/venvs/protect

Activate the virtualenv

    source ~/venvs/protect/bin/activate

NOTE: Installation was tested using pip 7.1.2 and 8.1.1. We have seen issues with the installation
of pyYAML with lower versions of pip and recommend upgrading pip before installing ProTECT.

    pip install --upgrade pip

Install Toil

    pip install toil[aws]==3.3.0

Install ProTECT and all dependencies in the virtualenv

    pip install protect

###Method 2 - Installing from Source

This will install ProTECT in an editable mode.

Obtain the source from Github

    git clone https://www.github.com/BD2KGenomics/protect.git

Create a virtualenv in the project folder

    cd protect
    virtualenv venv

Activate the virtualenv

    source venv/bin/activate

Install Toil

    pip install toil[aws]==3.3.0

Install ProTECT

    make develop


#Running ProTECT

Running ProTECT requires you to be in the virtualenv where it was installed.

Activate the virtualenv with

    source ~/venvs/protect/bin/activate
    OR
    source ~/protect/venv/bin/activate

Running ProTECT is as simple as:

            ProTECT --config /path/to/config.txt --workDir /path/to/a/workdir JobStore

Where

a) /path/to/config.txt is the absolute path to the config file for the run.   If a relative path is
provided, ProTECT will attempt to make it an absolute path, however this **might** be an issue on
OSX because folders like `/tmp` are symlinks to `/private/tmp`.

b) /path/to/a/workdir points to a directory to use as a working directory. This directory is where
all the worker sandboxes (temp directories) will be set up during the run. If this option is omitted,
Toil will setup workdir in the sytem $TMPDIR

c) JobStore is the job store used for the Toil run. The job store can either be a path to a
directory on the system (e.g. /home/probox/toil/jobstore), or an aws bucket as
aws:< region >:< bucket_name > (e.g. aws:us-west-2:ProTEST). Refer to the Toil manual for more
information regarding job stores.

**NOTE:** The workdir **MUST** be on a device that has enough space to run a genomics pipeline.
STAR indexes are usually 25Gb and running star usually takes ~60Gb per sample. You need enough space
on this device such that concurrent jobs do not run into `Unable to free up enough space for caching`
errors. If the job store being used is the file job store (points to a directory on the system), and
if the device hosting the job store is the same one hosting the working directory, you need to
ensure the device can handle the concurrent writes and the persistent job store files as well.

**NOTE:** It is advisable to setup Docker such that it runs off a large volume as well. The
following link contains information on how to setup docker such that it stores images and containers
on another volume.

    http://stackoverflow.com/questions/24309526/how-to-change-the-docker-image-installation-directory

 **TL;DR:** Use large volumes for workdir, job store and Docker.

#Setting up a config file

A config file pre-filled with references and options for an HG19 run can be generated with

    ProTECT --generate_config

The file `ProTECT_config.yaml` will be created in the current working directory.

The config file contains all the relevant information required for a run. The file is written in the
[YAML format](http://www.yaml.org/start.html) and should be relatively simple to fill in. The file
is split into a number of element groups that describe the important flags passed to different tools
in the pipeline, and the information on the input samples. Elements before a `:` are keys in the
dictionary read into ProTECT and should **NOT** be modified (Barring the patient ID key in the
patients dictionary). Only values to the right of the `:` should be edited.

Every required reference file is provided in the AWS bucket `cgl-protect-data` under the folder
`hg19_references`. The `README` file in the same location describes in detail how each file was
generated.

The following section describes the arguments that can be provided to ProTECT via the config file.

**patients**

The first group of elements is the list of patients to be processed in the run. Each patient group
is marked with a patient ID, and has 3 values which correspond to the forward fastq read file for
the Tumor DNA, Normal DNA, and Tumor RNA samples respectively. There is no limit to the number of
patients that can be run in the same workflow. The paths to the forward read containing file for the
sample types can either be absolute links to files on the system, or S3 links. The file extensions
can be .fq, .fq.gz. .fastq, or .fastq.gz and the file names should be symmetrically named
< prefix >1.extn and < prefix >2.extn where extn is one of the 4 available extensions mentioned
above.

    patients:    -> This is the group name. Do not modify this
        PRTCT-01:   -> A string describing a unique Identifier for the patient. This is the only key
                       that can be modified by the user.
            tumor_dna_fastq_1: /path/to/Tum_1.fq
            normal_dna_fastq_1: /path/to/Norm_1.fq.gz
            tumor_rna_fastq_1: S3://databucket/datafolder/Rna_1.fq.gz

Encrypted data will be decrypted using per-file SSE-C keys hashed from a master sse-key
(see `Universal_Options`). The backend download and upload from S3 is done using [s3am](https://www.github.com/BD2KGenomics/s3am.git).

**Universal_Options**

These describe options that are used universally by most tools/jobs in the workflow.

    Universal_Options:
        dockerhub: aarjunrao              -> All tools used in the pipeline are dockerized to allow
                                             for easy reproducibility across operating systems.
                                             ProTECT was tested using the hub owned by `aarjunrao`
                                             but will work with any hub given it contains the
                                             required tools described in
                                             `required_docker_tools.txt`.
        java_Xmx: 20G                     -> The default Java heap space to be provided to tools.
                                             Per-tool heap space can be specified for some tools to
                                             override this value.
        sse_key: /path/to/master.key      -> Used to create per-file SSE-C keys for decrypting S3
                                             hosted input files.  It is highly recommended that the
                                             files be uploaded to S3 using s3am using the
                                             --sse-key-is-master flag.
        sse_key_is_master: True           -> Were the sample files all encrypted with sse_key
                                             (`False`), or were they encrypted with individual
                                             per-file keys hashed from the master sse_key (`True`)
        cghub_key: /path/to/cghub.key     -> The CGHub credentials to use if the files must be
                                             downloaded from CGHub (This will be phased out soon)
        storage_location: aws:protect-out -> Should the intermediate files be stored locally
                                             (`Local`) or in an Amazon S3 bucket using SSE-C
                                             per-file encryption from the provided master sse key
                                             (`aws:<bucket_name>`)?
        output_folder: /path/to/out       -> A path to a folder where intermediate, and output files
                                             (mutation calls, MHC haplotype, etc). This folder will
                                             be created if it doesn't exist. However, files in the
                                             folder from a previous run will not be overwritten with
                                             updated values.


In reality, for the inexperienced user, these are the only values that need to be modified for a
run.  Experienced users can provide their own references below for better control over the results.
Look at the help for each tool for details. All ProTECT provided indexes were generated using
[Gencode release 19 for hg19](http://www.gencodegenes.org/releases/19.html). Any and all paths can
be substituted with S3 links. Descriptions for creating all files can be found in the file
`S3://cgl-protect-data/hg19_references/README`.

**cutadapt**

These flags describe the sequences to use for adapter trimming in the RNA Seq data

    cutadapt:
        a: AGATCGGAAGAG  -> Adapter for trimming forward read
        A: AGATCGGAAGAG  -> Adapter for trimming reverse read


**star**

These flags describe the arguments for aligning the RNA Seq data

    star:
        type: star                             -> Can be `star` or `starlong` depending on the read
                                                  size. Refer to the STAR manual for details.
        index_tar: /path/to/star_index.tar.gz  -> The indexes to use for STAR. Protect provides
                                                  indexes created using edge sizes 50, 75, 100, 150
                                                  and 250 in the S3 bucket `cgl-protect-data` under
                                                  the folder `hg19_references`.

**bwa**

These flags describe the arguments for aligning the DNA Seq data

    bwa:
        index_tar: /path/to/bwa_index.tar.gz  -> The indexes to use for bwa. Protect provides
                                                 indexes in the S3 bucket `cgl-protect-data` under
                                                 the folder `hg19_references`.


**rsem**

These flags describe the arguments for conducting expression calling.

    rsem:
        index_tar: /path/to/rsem_index.tar.gz  -> The indexes to use for rsem. Protect provides
                                                  indexes in the S3 bucket `cgl-protect-data` under
                                                  the folder `hg19_references`.

**mut_callers**

These flags describe the arguments for conducting mutation calling.

    mut_callers:
        genome_fasta: /path/to/hg19.fa.tar.gz                 -> The genome fasta to use for the
                                                                 mutation callers.
        genome_fai: /path/to/hg19.fa.fai.tar.gz               -> The corresponding .fai file for the
                                                                 genome fasta
        genome_dict: /path/to/hg19.dict.tar.gz                -> The corresponding .dict file for
                                                                 the genome fasta
        cosmic_vcf: /path/to/CosmicCodingMuts.vcf.tar.gz      -> The Cosmic Coding vcf
        cosmic_idx: /path/to/CosmicCodingMuts.vcf.idx.tar.gz  -> The corresponding .idx file for the
                                                                 Cosmic vcf
        dbsnp_vcf: /path/to/dbsnp_coding.vcf.gz               -> The dbSNP vcf
        dbsnp_idx: /path/to/dbsnp_coding.vcf.idx.tar.gz       -> The corresponding .idx file for the
                                                                 dbSNP vcf
        dbsnp_tbi : /path/to/dbsnp_coding.vcf.gz.tbi          -> The tabix index for dbsnp.gz
        java_Xmx: 5G                                          -> The heap size to use for MuTect
                                                                 per job (i.e. per chromosome)
        strelka_config: /path/to/strelka_config.ini.tar.gz    -> The Strelka config file for a bwa
                                                                 run (modified for a WXS run if
                                                                 necessary)


**snpeff**

These flags describe the arguments for conducting mutation calling.

    snpeff:
        index_tar: /path/to/snpeff_index.tar.gz  -> The indexes to use for snpeff. Protect provides
                                                    indexes in the S3 bucket `cgl-protect-data`
                                                    under the folder `hg19_references`.
        java_Xmx: 5G                             -> The heap size to use for SnpEFF

**transgene**

These flags describe the arguments for conducting mutation translation. Transgene is our in-house
tool to go from mutations to peptides.

    gencode_peptide_fasta: /path/to/gencode.faa   -> The corresponding peptide file for the gencode
                                                     gtf used for the analysis. (If not gencode, the
                                                     file must be made to follow the gencode format
                                                     for fasta record names)

**phlat**

These flags describe the arguments for conducting MHC Haplotyping.

    phlat:
        index_tar: /path/to/snpeff_index.tar.gz  -> The indexes to use for PHLAT. Protect provides
                                                    indexes in the S3 bucket `cgl-protect-data`
                                                    under the folder `hg19_references`.

**mhci**

These flags describe the arguments for conducting MHCI:peptide binding prediction.

    mhci:
        method_file: /path/to/mhci_restrictions.json.tar.gz -> A json list of allowable MHCs
                                                               different predictors can handle.
        pred: IEDB_recommended                              -> The IEDB method to use.

**mhcii**

These flags describe the arguments for conducting MHCII:peptide binding prediction.

    mhcii:
        method_file: /path/to/mhcii_restrictions.json.tar.gz -> A json list of allowable MHCs
                                                                different predictors can handle.
        pred: IEDB_recommended                               -> The IEDB method to use.

**rank_boost**

These flags describe the arguments for ranking the peptide binding calls from the various
predictors.

    rank_boost:
        mhci_combo: V,W,X,Y,Z  -> Weights used for ranking the predicted MHCI epitopes
                                  (V+W+X+Y+Z = 1)
        mhcii_combo W,X,Y,Z    -> Weights used for ranking the predicted MHCII epitopes
                                  (W+X+Y+Z = 1)

**mhc_pathway_assessment**

These flags describe the arguments for assessing the state of the MHC pathway in the patient.

    mhc_pathway_assessment:
        genes_file: /path/to/mhc_pathway_genes.json.tar.gz -> A json file containing the various
                                                              genes in the MHC pathway, and their
                                                              TPM expressions in a background set.
