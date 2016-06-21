ProTECT is implemented in the [Toil](https://github.com/BD2KGenomics/toil.git) framework and fully
runs the workflow described in [protect/Flowchart.txt](
https://github.com/BD2KGenomics/protect/blob/master/Flowchart.txt).
*Create an Overview section that lists the inputs, outputs, and the programs run (see e.g. https://github.com/BD2KGenomics/toil-scripts/tree/master/src/toil_scripts/rnaseq_cgl/README.md). Then point to the flow chart for details*

#Installation

ProTECT requires Toil and we recommend installing ProTECT and its requirements in a [virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/).

###Method 1 - Using PIP (recommended)

First create a virtualenv at your desired location (Here we create it in the folder ~/venvs)

    virtualenv ~/venvs/protect

Activate the virtualenv

    source ~/venvs/protect/bin/activate

Install ProTECT and all dependencies in the virtualenv

    pip install protect

###Method 2 - Installing from Source

This will install ProTECT in an editable mode.

Obtain the source form Github

    git clone https://www.github.com/BD2KGenomics/protect.git

Create a virtualenv in the project folder

    cd protect
    virtualenv venvs

Activate the virtualenv

    source venvs/bin/activate

Install ProTECT

    make develop


#Running ProTECT

Running ProTECT requires you to be in the virtualenv where it was installed.

Activate the virtualenv with

    source ~/venvs/protect/bin/activate

Running ProTECT is as simple as:

            python ProTECT.py --config config.txt

#Setting up a config file

The config file contains all the relevant information required for a run. The file is written in the
[YAML format](http://www.yaml.org/start.html) and should be relatively simple to fill in. The file is split
into a number of element groups that describe the important flags passed to different tools in the
pipeline, and the information on the input samples. All elements before a `:` are keys in the
dictionary read into ProTECT and should **NOT** be modified. Only values to the right of the `:`
should be edited.


**patients**

The first group of elements is the list of patients to be processed in the run. Each patient group
is marked with a patient ID, and has 3 values which correspond to the forward fastq read file for
the Tumor DNA, Normal DNA, and Tumor RNA samples respectively. There is no limit to the number of
patients that can be run in the same workflow. The paths to the forward read containing file for the
sample types can either be absolute links to files on the system, or S3 links. The file extensions
can be .fq, .fq.gz. .fastq, or .fastq.gz and the filenames should be symmetrically named
< prefix >1.extn and < prefix >2.extn where extn is one of the 4 available extensions mentioned above.

    patients:    -> This is the group name. Do not modify this
        PRTCT-01:   -> A string describing a unique identifier for the patient
        *FIXME: In input_parameters.yaml, this is PATIENT_ID1. It's not clear if this value IS the patient ID (this violates the 'all elements before a `:` should not be modified' rule), or if it should be PATIENT_ID1: my_actual_patient_id*
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
                                             Per-tool heap space *FIXME, unclear:* as well overrides this value.
        sse_key: /path/to/master.key      -> Used to create per-file SSE-C keys for decrypting S3
                                             hosted input files.  It is highly recommeded that the
                                             files be uploaded to S3 using s3am using the
                                             --sse-key-is-master flag.
        sse_key_is_master: True           -> Were the sample files all encrypted with sse_key
                                             (`False`), or were they encrypted with individual
                                             per-file keys hashed from the master sse_key (`True`)
        cghub_key: /path/to/cghub.key     -> The CGHub credentials to use if the files must be
                                             downloaded from CGHub ( will be phased out soon)
        storage_location: aws:protect-out -> Should the intermediate files be stored locally
                                             (`Local`) or in an Amazon S3 bucket using SSE-C
                                             per-file encryption from the provided master sse key
                                             (`aws:<bucket_name>`)?
        output_folder: /path/to/out       -> A path to a folder where intermediate, and output files
                                             (mutation calls, MHC haplotype, etc). *FIXME: will folder be created if it doesn't exist? Will the program die if the folder already exists?*


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
        A: AGATCGGAAGAG  -> Adapter for trimming reverse read *JUST CHECKING: these really are the same sequence?*


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
        genome_fai: /path/to/hg19.fa.fai.tar.gz               -> The fai index for the genome fasta
        genome_dict: /path/to/hg19.dict.tar.gz                -> The dict for the genome fasta
        cosmic_vcf: /path/to/CosmicCodingMuts.vcf.tar.gz      -> The Cosmic Coding vcf
        cosmic_idx: /path/to/CosmicCodingMuts.vcf.idx.tar.gz  -> The index for the Cosmic vcf
        dbsnp_vcf: /path/to/dbsnp_coding.vcf.tar.gz           -> The dbSNP vcf
        dbsnp_idx: /path/to/dbsnp_coding.vcf.idx.tar.gz       -> The idx for the dbSNP vcf
        java_Xmx: 5G                                          -> The heap size to use for MuTect
                                                                 per job (i.e. per chromosome)

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
