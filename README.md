Nanopore paper analysis
=======================

This repository contains the pipeline we used for our [2015 Nature Methods article](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3444.html). The version tagged [paper](https://github.com/jts/nanopore-paper-analysis/tree/paper) is the code we ran to generate the results for our article. Since publication we have merged patches to clean up the code but no changes to functionality have been made.

The ```full-pipeline.make``` Makefile will download the input data from the ENA, download the versions of all required software and run the complete assembly pipeline. This pipeline requires GCC >= 4.8 and should be run within a Python virtualenv. It takes a single parameter, ```CORES```, which determines how many CPU cores will be used for the assembly. This command will generate the assembly:

    make -f full-pipeline.make CORES=64 polished_genome.fasta

The IPython notebook ```full-walkthrough``` in this repository demonstrates the use of the pipeline, starting from a stock Ubuntu 14.04 system. The notebook additionally contains the commands we ran to generate the results and figures for our paper.
