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
Pathseq takes its metagnomics sequencing data as paired-end reads in fastq/fasta or fastq.gz format. For each sample, the paired end reads must be given as one fasta file of forward reads and one fasta file of backward reads. \
\
Pathseq must also be given a configuration file to run. This file specifies all arguments (besides the cluster profile, which must be given in the run script. See below for "Running PathSeq on your own data" and "what the cluster profile means"), including the metagenomics sequencing reads you wish to analyze. \
More details about arguments required are below:

#### Input data requirements
The config file requires the following parameters to be filled in: \
1. "inputlist" \
     The full path to a text file containing a list of the file names of all forward-read input fasta files. Each file name should be on its own line. PathSeq treats each fasta file given as its own separate sample, and will run the full pipeline analysis separately for each sample. In other words, the output would be the same as if you ran the pipeline separately for each input file, but giving them to PathSeq together allows PathSeq to run them in parallel if enough cluster nodes are given (number of cluster nodes to use is specified in the profile. More detail in "What the cluster profile means").
3. "inputdir"
     The full path to a directory containing all the forward-read fasta files listed in "inputlist", as well as their backward-read counterparts. The directory can contain other files, including the "inputlist" file, but "inputlist" does not have to be kept inside "inputdir".
5. "origin"
     "RNA" or "DNA" depending on whether your input sequencing reads are from RNA or DNA.
7. "minContigLen"
     The minimum contig length for contig assembly (can stay the same as your read lengths, or be longer).
9. "codepath"
     The full path to the rule scripts, which in this repository are found inside "DouekLabPathSeq/SlurmBaseCode/". These are scripts which will be run as a part of each snakemake rule, and correspond to each step of the pipeline.
11. "outputpath"
     The full path to the directory where you would like house the outputs of the scripts. Within this directory, PathSeq will create separate sub-directories for the output of each sample you give it. More detail below in "Output format".
13. Reference databases
    The full paths to certain reference database locations that are required by some of the bioinformatics softwares used. More detail in "Other arguments needed".
15. External programs
    The full paths to certain executable files that need to be given to certain steps of the pipeline. More detail in "Other arguments needed".
17. Resource allocation
     If running on cluster mode, then for each rule, you must specify the maximum time allowed to the rule ("runtime", which is specified in minutes), and the maximum memory allowed to the rule ("mem_mb", in megabytes). This is so the cluster can properly launch batch jobs for each rule. In "DouekLabPathSeq/configTemplate/config.yaml", resource allocations have been given which we found generally worked well for a microbiome dataset with input fasta files generally under 2-3 GB. You may need to adjust these resource allocations based on your dataset (larger input files, higher coverage, and higher microbial diversity may all increase the amount of resources needed).

#### Output format
The config file's "outputpath" specifies the full path to the directory where the user wishes for output data to be written. \
1. For each separate sample, a separate subdirectory within the output directory will be created. \
2. Within each subdirectory, a directory for each rule is created. \

#### Other arguments needed
Some software must be given to rules in the snakefile as paths to executables, rather than as "module load" statements or as scripts in the "SlurmBaseCode" folder. \
There is no technical reason for doing so, it is mostly for legacy reasons. These executables are given in "DouekLabPathSeq/PathseqExternalPrograms/". Paths to each one must be supplied to the config.yaml file. \
In-house executables: \
1. PathSeqRemoveHostForKaiju
2. PathSeqKaijuConcensusSplitter2
3. PathseqmergeQIIME2TaxAndSalmon
4. PathseqSplitOutputTableByTaxonomy
Other executables: \
1. Picard (2.18.14)
2. Prodigal (2.6.3)
3. Kaiju (1.9.0)
\
Additionally, the following reference databases are required: \
1. ERCC92 spike-in controls for bowtie
2. human genome reference dataset index for bowtie
3. human genome reference dataset for star
4. primate genome reference dataset for bowtie
5. mammal genome database for blast
6. nodes and fmi for kaiju
7. reference taxonomy for kaiju

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
1. A configuration file specifying all arguments to the run. \
   There is a template at "DouekLabPathSeq/configTemplate/config.yaml". \
   You must make sure all paths are valid full paths in your system. Full paths allow you to house (potentially very large) input data, output data, and reference databases outside your code area. ALL runs require a configuration file. \
   In the configuration file, you can also set resource usage for each rule's batch job. "runtime" is the allowed job time specified in minutes, "mem_mb" is the allowed job memory usage in megabytes.
2. A run script. \
   There is a template at "DouekLabPathSeq/configTemplate/runSnakemake.sh". \
   This script takes a path to the Snakefile ("DouekLabPathSeq/SlurmBaseCode/Snakefile"), and a path to the profile ("DouekLabPathSeq/skyline_profile"). Depending on where you have stored the DouekLabPathSeq repo, you may need to customize these paths. They can be relative or full paths.
3. A profile file(s). An example profile specific to the Skyline HPC is located at "DouekLabPathSeq/skyline_profile/". \
   This defines what commands on your system are used to launch a batch job (i.e. how your system should translate the contents of the "resources" section of each rule in the Snakefile into a command that can be understand by your cluster's workload management software). This file is only required if you want to run the pipeline in cluster mode. \
   More details are above in the System requirements section.
4. A snakemake unlock script. A template is in "DouekLabPathSeq/configTemplate/unlock.sh". \
   This file should run the same snakefile as the snakemake command in "runSnakemake.sh" does. This script is only necessary if your snakemake run is interrupted in the middle (i.e. time out). In this case, it "unlocks" the snakefile (the snakemake process places a lock on the working directory at the start of a run, and this lock may become stale if the process is interrupted).

# Software details
## Program versions used in development:
snakemake 7.22.0 \
python 3.11

## Details on each rule and its dependencies:
1. BowtieUnmasked \
      Dependencies:
      1. openjdk (java) v 17.0.11
      2. 
2. StarAfterBowtie
3. BowtiePrimate
4. Trinity
5. FilterHostBlast
6. kaiju
7. BuildSalmon
8. SalmonQuant
9. MergeTaxAndQuant

## What the cluster profile means
