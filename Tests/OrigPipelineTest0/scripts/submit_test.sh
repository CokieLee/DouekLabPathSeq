#!/bin/bash
#SBATCH --partition=quick
#SBATCH -J pathSeq
#SBATCH --output=pathseq.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 1:00:00

PathSeqSubmitter.sh 582760806_A5_S33_R1_001_short_test ../Input_Data/582760806_A5_S33_R1_001_short_test.fastq.gz  ../Input_Data/582760806_A5_S33_R2_001_short_test.fastq.gz RNA 300
