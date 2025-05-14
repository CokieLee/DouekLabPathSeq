# Pathseq
Pathseq is a pipeline for identifying microbial content (bacteria or viruses) from metagenomic sequencing reads.\
\
Specifically, pathseq takes paired-end shotgun sequencing reads ->\
filters reads that are likely to be contamination (human, primate, mammalian) ->\
assembles remaining reads into contigs ->\
filters contigs for likely contamination ->\
taxonomically classifies and quantifies remaining contigs ->\
searches remaining unclassified contigs for RNA viral motifs (to find possible novel viruses)

# How to run pathseq
## System requirements
Pathseq is made to run on a Linux system. We have tested specifically on Linux kernel 5.14.0-284.30.1.el9_2.x86_64.\
We also provide dependencies with a conda environment. We use conda version 24.1.2.\
\
Pathseq can be run either on a single computer, or on a cluster.\
To run Pathseq on a cluster, such that each rule launches as its own job, you must define a "profile" for your cluster. A profile maps the concept of launching a rule as a cluster job to concrete commands run by your system.
The profile settings can be given to the snakemake command directly as flags, but it is recommended to persist the settings in a profile configuration file (which is passed to the snakemake command as the argument of the --profile flag).\
More information on snakemake 7 cluster execution and profiles: https://snakemake.readthedocs.io/en/v7.22.0/executing/cluster.html\
We provide a working profile configuration file in "skyline_profile/config.yaml". This configuration file is tested on our Linux cluster, running Slurm versino 23.02.6. It may need to be modified for other systems. More detail on explaining and modifying this file is in the "Software dependency details" section below.\
\
To run Pathseq on a single computer, no profile setting files are necessary.


TODO:
- explain how the sample fastqs must be named, how their names are used, and how output folders are named.
## Input requirements and output format
Pathseq must also be given a configuration file to run. This file specifies all other arguments, including the input and output directory full paths, the sequencing origin (RNA or DNA), the mininum contig length (ordinarily set to the read length), path to the rule scripts (in SlurmBaseCode), other external programs used, reference databases, and a scratch directory (directory used to store temporary data while computing).\
More details about arguments required are below:\
Pathseq takes its metagnomics sequencing data as paired-end reads in fastq.gz format (gzipped fastq).\
\
#### Input data
The config file's "inputpath" parameter specifies the path to a directory containing any such fastq.gz files. This directory must also contain a text file "file_list", which lists the file names of each fastq.gz file you wish to input to Pathseq. Pathseq treats all fastq.gz files listed on "file_list" as separate samples, and will run the complete pipeline for each sample separately, and output information related to each sample into separate output subdirectories, within a main output subdirectory. \
\
#### Output
The config file's "outputpath" specifies the full path to the directory where the user wishes for output data to be written. For each separate sample, a separate subdirectory within the output directory will be created. Within each subdirectory, a directory for each rule is created. More detail on the outputs of each rule can be foudn below in the Software details section.\
\
#### Other arguments needed
TODO: reference databases, external software.

## Testing your installation
1. Install 
2. There is a folder, "Tests/CondaEnvTest", which contains a test config file

## Running Pathseq on your own data
To summarize, there are two files which may be customized in order to run on your own data.
1. A configuration file specifying all arguments to the run (test example at "Tests/CondaEnvTest/config.yaml"). You must make sure all paths are valid full paths in your system. Full paths you to house (potentially very large) input data, output data, and reference databases outside your code area. ALL runs require a configuration file.
2. A profile file (test example at "skyline_profile/config.yaml") defining what commands on your system are used to launch a batch job. This file is only required if you want to run the pipeline in cluster mode. More details are above in the System requirements section.

# Software details
## Programs contained in global conda environments:
snakemake 7.22.0
python 3.11
## Details on each rule, dependencies and output:
## What the cluster profile means, and how to modify it
