#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=175G
#SBATCH -t 24:00:00

module load snakemake
snakemake --version

snakemake --unlock

snakemake -p -c1 --latency-wait 120

echo "RAN SNAKEMAKE"