#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G

module load snakemake

snakemake -p -c1 --latency-wait 120

echo "RAN SNAKEMAKE"