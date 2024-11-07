#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH -t 24:00:00

module load snakemake

snakemake -p -c1 --latency-wait 240

echo "RAN SNAKEMAKE"