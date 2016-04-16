ProTECT is implemented in the [Toil](https://github.com/BD2KGenomics/toil.git) framework and fully
runs the workflow described in [protect/Flowchart.txt](
https://github.com/BD2KGenomics/protect/blob/master/Flowchart.txt).

Running ProTECT is as simple as:

            python ProTECT.py --config config.txt

The config file contains all the relevant information required for a run.

The file is read starting at the `BEGIN` keyword. Following the `BEGIN` keyword are a number of
element groups that describe the important flags passed to different tools in the pipeline. Each
group is prefixed by `#<groupname>` and it is important to not modify and/or delete any of them.

**# patient**

The first group of elements is the list of patients to be processed in the run. Each patient group
is preceded by `#patient` and contains 4 values. There is no limit to the number of patients that
can be run in the same workflow.

    patient_id              - A string describing a unique Identifier for the patient
    tumor_dna_fastq_prefix  - The path to the forward read containing file for the tumor dna
                              (WGS or WXS)
    normal_dna_fastq_prefix - The path to the forward read containing file for the normal dna
                              (WGS or WXS)
    tumor_rna_fastq_prefix  -  The path to the forward read containing file for the tumor rna

The allowed values for the 3 input files are `XYZ_1.fastq`, `XYZ_1.fastq.gz`, `XYZ_1.fq`,
`XYZ_1.fz.gz`, or `XYZ.xml` (Only for samples to be downloaded from CGHub).  All files can also
either be local files on the system, or links to files in Amazon S3 buckets (should be of the form
 https://s3-us-west-2.amazonaws.com/\<BUCKET\>/KEY).  Encrypted data will be decrypted using
per-file SSE-C keys hashed from a master sse-key (see `Universal Options`).

**# Universal_Options**

These describe options that are used universally by almost all tools/jobs in the workflow.

    dockerhub               - All tools used in the pipeline are dockerized to alow for easy
                              reproducibility across operating systems. ProTECT was rested using the
                              hub owned by `aarjunrao` but will work with any hub given it contains
                              the required tools.
    java_Xmx                - The default Java heap space to be provided to tools. You can specify
                              per-tool heap space as well.
    sse_key                 - Used to create per-file SSE-C keys for decrypting S3 hosted input
                              files.  It is highly recommeded that the files be uploaded to S3 using
                              [s3am](https://github.com/BD2KGenomics/s3am.git) using the
                              --sse-key-is-master flag.
    cghub_key               - The CGHub credentials to use if the files must be downloaded from
                              CGHub
    storage_location        - Should the intermediate files be stored locally (use `Local`) or in an
                              Amazon S3 bucket using SSE-C per-file encryption from the provided
                              master sse key (use `aws:<bucket_name>`)?
    output_folder           - A path to a folder where intermediate, and output files (mutation
                              calls, MHC haplotype, etc).


In reality, for the inexperienced user, those are the only values that need to be modified for a
run.  Experienced users can provide their own references for better control over the results. Look
at the help for each tool for details. All ProTECT provided indexes were generated using [Gencode
release 19 for hg19](http://www.gencodegenes.org/releases/19.html).

**# cutadapt**

The flags describe the sequences to use for adapter trimming in the RNA Seq data


**# star**

    type                    - Can be `star` or `starlong` depending on the read size. Refer the STAR
                              manual for details.
    index_tar               - The Indexes to use for STAR.  Protect provided indexes created using
                              edge sizes 100, 150 and 250 in the S3 bucket.

**# bwa**

    index_tar - The indexes to use for BWA.

**# rsem**

    index_tar - The indexes to use for RSEM.

**# mut_callers**

    genome_fasta            - The genome fasta to use for the mutation callers.
    genome_fai              - The fai for the genome fasta
    genome_dict             - The dict for the genome fasta (Refere MuTect)
    cosmic_vcf              - The Cosmic Coding vcf
    cosmic_idx              - The idx for the Cosmic vcf
    dbsnp_vcf               - The dbSNP vcf
    dbsnp_idx               - The idx for the dbSNP vcf
    java_Xmx                - The Heap size to use for MuTect per job (i.e. per chromosome)

**# snpeff**

    index_tar               - The indexes to use for SnpEFF.
    java_Xmx                - The Heap size to use for SnpEFF

**# transgene**

Transgene is our in-house tools to go from mutations to peptides.

    gencode_peptide_fasta   - The corresponding peptide file for the gencode gtf used for the
                                analysis. (If not gencode, the file must be made to follow the
                                gencode format for fasta record names)

**# phlat**

    index_tar               - The indexes to use for PHLAT.

**# mhci**

    method_file             - A json list of allowable methods different MHCs can handle.
    pred                    - The IEDB method to use (`IEDB_recommended`)

**# mhcii**

    method_file             - A json list of allowable methods different MHCs can handle.
    pred                    - The IEDB method to use (`IEDB_recommended`)

**# rank_boost**

    mhci_combo              - Weights used for ranking the predicted MHCI epitopes
    mhcii_combo             - Weights used for ranking the predicted MHCII epitopes