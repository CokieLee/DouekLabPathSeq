#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH -t 4:00:00

module load snakemake/7.22.0

snakemake -s ../../SlurmBaseCode/Snakefile --unlock

echo "RAN SNAKEMAKE"
