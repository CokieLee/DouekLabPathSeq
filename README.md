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
Currently, PathSeq is made to run on a linux cluster (Rocky Linux v9.5) which uses the workload manager Slurm (version 23.02.6), and a module system. \
\
Plans are in development to containerize rules such that the pipeline is HPC agnostic.
So, to run on a non-Skyline HPC, you may need to change the "module load" statements, or else otherwise ensure the correct dependencies are available to each rule.
\
Pathseq can be run either on a single computer, or on a cluster. \
To run Pathseq on a cluster, such that each rule launches as its own job, you must define a "profile" for your cluster. A profile maps the concept of launching a rule as a cluster job to concrete commands run by your system.
The profile settings can be given to the snakemake command directly as flags, but it is recommended to persist the settings in a profile configuration file (which is passed to the snakemake command as the argument of the --profile flag).\
More information on snakemake 7 cluster execution and profiles: https://snakemake.readthedocs.io/en/v7.22.0/executing/cluster.html \
We provide a working profile configuration file in "skyline_profile/config.yaml" (specific to our HPC). \
\
To run Pathseq on a single computer, no profile setting files are necessary. You must also remove the --profile flag from the "runSnakemake.sh" script when running.

## Input and output formatting
Pathseq must also be given a configuration file to run. This file specifies all other arguments, including the input and output directory full paths, the sequencing origin (RNA or DNA), the mininum contig length (ordinarily set to the read length), path to the rule scripts (in SlurmBaseCode), other external programs used, reference databases, and a scratch directory (directory used to store temporary data while computing).\
More details about arguments required are below:\
Pathseq takes its metagnomics sequencing data as paired-end reads in fastq.gz format (gzipped fastq).\

#### Input data requirements
The config file's "inputpath" parameter specifies the path to a directory containing any such fastq.gz files. This directory must also contain a text file "file_list", which lists the file names of each fastq.gz file you wish to input to Pathseq. Pathseq treats all fastq.gz files listed on "file_list" as separate samples, and will run the complete pipeline for each sample separately, and output information related to each sample into separate output subdirectories, within a main output subdirectory. \

#### Output format
The config file's "outputpath" specifies the full path to the directory where the user wishes for output data to be written. For each separate sample, a separate subdirectory within the output directory will be created. Within each subdirectory, a directory for each rule is created. More detail on the outputs of each rule can be foudn below in the Software details section.\

#### Other arguments needed
TODO: reference databases, external software.

## Testing your installation on Skyline
1. After installing Pathseq
2. In the folder "Tests/TestInput/", there is a set of left and right paired-end fastq files, containing about 125 reads.
3. There is a folder, "Tests/Skyline_Full_Test_00/". This folder contains:
     a. A "config.yaml" file. This file already contains all input arguments required to run pathseq on the "Tests/TestInput/" input data.
     b. A "runSnakemake.sh" bash script. This bash script uses the "config.yaml" file in the "Skyline_Full_Test_00" directory, and runs Pathseq on it, by running snakemake with "SlurmBaseCode/Snakefile" (the central snakefile).
4. Run the following command from within the SkylineTest0 directory to start the pipeline on your current compute node:
   ```./runSnakemake.sh```
   Or, run it as a batch job on your cluster however you normally submit bash scripts as batch jobs (you may need to change the resource allocation comments at the top to match your system, and the profile contents).
   ```sbatch runSnakemake.sh```
5. Once the code is finished running, check that the outputs in "Tests/Skyline_Full_Test_00/Output/" matches the contents of "Tests/Skyline_Full_Test_00/ExpectedOutput/"

## Running Pathseq on your own data
There are three files which may be customized in order to run on your own data.
1. A configuration file specifying all arguments to the run. There is a template at "DouekLabPathSeq/configTemplate/config.yaml". \
   You must make sure all paths are valid full paths in your system. Full paths allow you to house (potentially very large) input data, output data, and reference databases outside your code area. ALL runs require a configuration file. \
   In the configuration file, you can also set resource usage for each rule's batch job. "runtime" is the allowed job time specified in minutes, "mem_mb" is the allowed job memory usage in megabytes.
2. A run script. There is a template at "DouekLabPathSeq/configTemplate/runSnakemake.sh". \
   This script takes a path to the Snakefile ("DouekLabPathSeq/SlurmBaseCode/Snakefile"), and a path to the profile ("DouekLabPathSeq/skyline_profile"). Depending on where you have stored the DouekLabPathSeq repo, you may need to customize these paths. They can be relative or full paths.
3. A profile file(s). An example profile specific to the Skyline HPC is located at "DouekLabPathSeq/skyline_profile/". \
   This defines what commands on your system are used to launch a batch job (i.e. how your system should translate the contents of the "resources" section of each rule in the Snakefile into a command that can be understand by your cluster's workload management software). This file is only required if you want to run the pipeline in cluster mode. \
   More details are above in the System requirements section.
4. A snakemake unlock script. A template is in "DouekLabPathSeq/configTemplate/unlock.sh". \
   This file should run the same snakefile as the snakemake command in "runSnakemake.sh" does. This script is only necessary if your snakemake run is interrupted in the middle (i.e. time out). In this case, it "unlocks" the snakefile (the snakemake process places a lock on the working directory at the start of a run, and this lock may become stale if the process is interrupted).

# Software details
## Program versions used in development:
snakemake 7.22.0
python 3.11
## Details on each rule, dependencies and output:
## What the cluster profile means, and how to modify it