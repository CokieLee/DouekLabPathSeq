#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH -t 24:00:00

module load snakemake
which snakemake

snakemake -p -c $SLURM_CPUS_PER_TASK --latency-wait 240

echo "RAN SNAKEMAKE"