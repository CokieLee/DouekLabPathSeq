#!/bin/bash
#SBATCH --cpus-per-task=16

module load star

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir /data/vrc_his/douek_lab/reference_sets/hg38/Sequence/STAR/ \
     --genomeFastaFiles /data/vrc_his/douek_lab/reference_sets/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --sjdbGTFfile /data/vrc_his/douek_lab/reference_sets/hg38/Sequence/STAR/ \
     --sjdbOverhang 299
