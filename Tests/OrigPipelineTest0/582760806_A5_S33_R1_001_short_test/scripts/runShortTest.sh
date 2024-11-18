#!/bin/bash
#SBATCH --partition=quick
#SBATCH -J pathSeq
#SBATCH --output=pathseq.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 1:00:00

./PathSeqSubmitter.sh
