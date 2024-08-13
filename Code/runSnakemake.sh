#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g

module load snakemake

snakemake -p -c1

echo "RAN SNAKEMAKE"