

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

ProTECT also requires [s3am](https://github.com/BD2KGenomics/s3am.git) version 2.0.1 to download and
upload files from S3. We recommend installing s3am in its own virtualenv using the directions in
the s3am manual, then putting the s3am binary on your $PATH.  ProTECT will NOT attempt to install
s3am during installation.

Lastly, ProTECT uses [docker](https://www.docker.com/) to run the various sub-tools in a
reproducible, platform independent manner. ProTECT will NOT attempt to install docker during
installation.

###Method 1 - Using PIP (recommended)

First create a virtualenv at your desired location (Here we create it in the folder ~/venvs)

    virtualenv ~/venvs/protect

Activate the virtualenv

    source ~/venvs/protect/bin/activate

NOTE: Installation was tested using pip 7.1.2 and 8.1.1. We have seen issues with the installation
of pyYAML with lower versions of pip and recommend upgrading pip before installing ProTECT.

    pip install --upgrade pip

Install Toil

    pip install toil[aws]==3.5.2

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

    pip install toil[aws]==3.5.2

Install ProTECT

    make develop

##Method 3 - Using Docker

Dockerized versions of ProTECT releases can be found at https://quay.io/organization/ucsc_cgl. These
Docker containers run the ProTECT pipeline in single machine mode. The only difference between the
Docker and Python versions of the pipeline is that the Docker container takes the config options,
described below, as command line arguments as opposed to a config file. Running the container
without any arguments will list all the available options.

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
The protect.cwl file in the /docker/ directory describes the necessary inputs and outputs and should be used if running via Dockstore.

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
        reference_build: hg19             -> The reference build used in this run. Can be hg19,
                                             hg38, GRCh37 or GRCh38.
        sse_key: /path/to/master.key      -> Used to create per-file SSE-C keys for decrypting S3
                                             hosted input files.  It is highly recommended that the
                                             files be uploaded to S3 using s3am using the
                                             --sse-key-is-master flag.
        sse_key_is_master: True           -> Were the sample files all encrypted with sse_key
                                             (`False`), or were they encrypted with individual
                                             per-file keys hashed from the master sse_key (`True`)
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
            version: 2.1.0

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
            verison: 3.1

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
    mhc_pathway_assessment:
        genes_file: /path/to/mhc_pathway_genes.json.tar.gz        -> A json file containing the
                                                                     various genes in the MHC
                                                                     pathway, and their mean TPM
                                                                     expressions across samples in
                                                                     a background set.

#Default values for config entries (Only in source installations)

If you installed from source, we provide a way to specify default values for all entries in the
config file (except for the `patients` tab). The file `src/protect/pipeline/defaults.yaml` will be
used to fill in default values for config entries if the user does not specify them in the input
config file.


#A note on dockerised tools used in pipeline

ProTECT uses dockerised versions of every tool used during the run to ensure reproduciblity of
results. Every docker image required for the run is described in required_docker_tools.txt. Every
required tool is also hosted freely on the dockerhub `aarjunrao` (Hence the default value in the
config). If you wish to use a personal repo instead, ensure that every required version of every
tool is made available on the repo.