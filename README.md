[![Stories in Ready](https://badge.waffle.io/BD2KGenomics/protect.png?label=ready&title=Ready)](https://waffle.io/BD2KGenomics/protect)
# ProTECT
### **Pr**ediction **o**f **T**-Cell **E**pitopes for **C**ancer **T**herapy

This repo contains the Python script for the Precision Immunology Pipeline developed at UCSC.

    ProTECT.py             - The python script for running the pipeline.
    ProTECT_large.py       - The python script for running the pipeline on larger input files
                             (>15Gb/fastq).
    input_parameters.list  - The config file for the run that contains all the necessary parameters
                             for the run.
    Flowchart.txt          - A (super cool) flowchart describing the flow of the pipeline.


All docker images used in this pipeline are available at

                        https://hub.docker.com/u/aarjunrao/


To learn how the pipeline can be run on a sample, head over to the [ProTECT Manual](
https://github.com/BD2KGenomics/protect/blob/master/MANUAL.md)

ProTECT is currently in its infancy and is under continuous development.  We would appreciate users sharing the level 3 data produced by ProTECT with us such that we can better train our predictive models.
