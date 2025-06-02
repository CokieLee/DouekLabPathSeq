#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH -t 4:00:00

module load snakemake/7.22.0

snakemake testRuleAll -s ../../SlurmBaseCode/Snakefile \
	  --profile /data/vrc_his/douek_lab/projects/PathSeq/Krystelle/skyline_profile/ \
	  -p -c $SLURM_CPUS_PER_TASK --latency-wait 240

echo "RAN SNAKEMAKE"
