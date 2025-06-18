#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -t 15-00:00:00

module load snakemake/7.22.0
which snakemake

snakemake --snakefile /data/vrc_his/douek_lab/projects/PathSeq/cokie_pathseq_rewrite/DouekLabPathSeq/SlurmBaseCode/Snakefile \
	  --profile /data/vrc_his/douek_lab/projects/PathSeq/cokie_pathseq_rewrite/DouekLabPathSeq/skyline_profile/ \
	  -p -c $SLURM_CPUS_PER_TASK --latency-wait 120

echo "RAN SNAKEMAKE"
