Nanopore paper analysis
=======================

This repository contains the pipeline we used for our E. coli nanopore paper. The full-pipeline.make Makefile will download the input data from the ENA, download the versions of all required software and run the complete assembly pipeline. This pipeline requires GCC >= 4.8 and should be run within a Python virtualenv. It takes a single parameter, ```CORES```, which determines how many CPU cores will be used for the assembly. This command will generate the assembly:

    make -f full-pipeline.make CORES=16 polished_assembly.fasta

The IPython notebook ```full-walkthrough``` in this repository demonstrates the use of the pipeline, starting from a stock Ubuntu 14.04 system. The notebook additionally contains the commands we ran to generate the results and figures for our paper.
