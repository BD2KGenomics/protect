[![Stories in Ready](https://badge.waffle.io/BD2KGenomics/protect.png?label=ready&title=Ready)](https://waffle.io/BD2KGenomics/protect)
# ProTECT
### **Pr**ediction **o**f **T**-Cell **E**pitopes for **C**ancer **T**herapy

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
