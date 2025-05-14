# Pathseq
Pathseq is a pipeline for identifying microbial content (bacteria or viruses) from metagenomic sequencing reads.\
\
Specifically, pathseq takes paired-end shotgun sequencing reads ->\
filters reads that are likely to be contamination (human, primate, mammalian) ->\
assembles remaining reads into contigs ->\
filters contigs for likely contamination ->\
taxonomically classifies and quantifies remaining contigs ->\
searches remaining unclassified contigs for RNA viral motifs (to find possible novel viruses)

## How to run pathseq
#### System requirements
Pathseq is made to run on a Linux system. We have tested specifically on Linux kernel 5.14.0-284.30.1.el9_2.x86_64.\
We also provide dependencies with a conda environment. We use conda version 24.1.2.\
Pathseq can be run either on a single computer, or on a cluster. To run Pathseq on a cluster, such that each rule launches as its own job,
you must define a "profile" for your cluster. A profile maps the concept of launching a rule as a cluster job to concrete commands run by your system.
The profile settings can be given to the snakemake command directly as flags, but it is recommended to persist the settings in a profile configuration file (which is passed to the snakemake command as the argument of the --profile flag).
We provide a working profile configuration file in "skyline_profile/config.yaml". This configuration file is tested on our Linux cluster, running Slurm versino 23.02.6. It may need to be modified for other systems. More detail on modifying this file is in the "Software dependency details" section below.

#### Input requirements and output format
Pathseq takes paired-end reads in fastq.gz format (gzipped fastq).\
PathSeq can be given separate samples at once, and will return separate output directories for each such sample. To give Pathseq separate samples, the reads for each sample must be in separate fastq.gz files, and must be housed in the same directory as a "file_list" file, naming the fastq.gz file name of each separate sample file.\
Pathseq must also be given a 
#### Testing your installation
1. Install 
2. There is a folder, "Tests/CondaEnvTest", which contains a test config file

## Software dependency details
#### Programs contained in global conda environments:
snakemake 7.22.0
python 3.11
#### Programs contained in 
#### What the cluster profile means, and how to modify it
