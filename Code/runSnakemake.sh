#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

module load snakemake

snakemake -p -c1 --latency-wait 30

echo "RAN SNAKEMAKE"