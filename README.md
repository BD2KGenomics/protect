# ProTECT
### **Pr**ediction **o**f **T**-Cell **E**pitopes for **C**ancer **T**herapy

Adapation of ProTECT to use python 3.8 instead of 2.7. Currently have tested a complete run using fastq files from [HCC1395 WGS Exome RNA Seq Data](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data), but have not checked results against the [original ProTECT](https://github.com/BD2KGenomics/protect) with TCGA PRAD yet. 

Adaptation done using 2to3 and manual bug testing. Manual changes recorded [at changes.md](https://github.com/Dranion/protect/blob/master/changes.md). Continuing to the original README: 

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
