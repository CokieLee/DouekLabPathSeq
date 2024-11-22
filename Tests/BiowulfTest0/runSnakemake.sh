#!/bin/bash
#SBATCH --partition=quick
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --gres=lscratch:500
#SBATCH --cpus-per-task=16
#SBATCH --mem=175G
#SBATCH -t 1:00:00

module load snakemake
snakemake --version

# snakemake --cleanup-metadata
locale

snakemake -p -c $SLURM_CPUS_PER_TASK --latency-wait 120

echo "RAN SNAKEMAKE"