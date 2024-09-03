#!/bin/sh
#SBTACH -J starAfterBowtie
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

module load star
module load samtools

COUNTER=$SGE_TASK_ID

leftUnalignedFile=$1
rightUnalignedFile=$2
left_read_file_base_name=$3
codePath=$4
outputPath=$5
hg38_starDB=$6

## print input args
echo "STAR AFTER BOWTIE INPUTS:"

echo "leftUnalignedFile:"
echo $leftUnalignedFile
echo "rightUnalignedFile:"
echo $rightUnalignedFile

echo "output_base_name:"
echo $left_read_file_base_name

echo "codePath:"
echo $codePath

echo "hg38_starDB:"
echo $hg38_starDB

echo "outputPath:"
echo $outputPath

## Source script for directory checking function
dos2unix $codePath"/dir_check.sh"
source $codePath"/dir_check.sh"

##Change directories to output folder
cd $outputPath

mkdir "alignment_stats"
mkdir Generated_Data_Star_Alignment

cd Generated_Data_Star_Alignment/

## Confirm that we are in the folder Generated_Data_Star_Alignment
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="Generated_Data_Star_Alignment"
dir_check  $correct_cur_Dir

echo "CHECKPOINT 1: star command:"
starCmd="STAR --genomeDir $hg38_starDB --readFilesIn $leftUnalignedFile $rightUnalignedFile  --outFileNamePrefix ./$left_read_file_base_name --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt"
echo starCmd

echo "$starCmd" | bash

echo "STAR run finished"

bam_out_file_name=${left_read_file_base_name}"Aligned.sortedByCoord.out.bam"
log_out_file_name=${left_read_file_base_name}"Log.final.out"
echo "BAM out file name"
echo $bam_out_file_name
echo "Log.out file name"
echo $log_out_file_name
samtools index $bam_out_file_name
samtools idxstats $bam_out_file_name > ${left_read_file_base_name}"_star_alignment_stats.txt"
cp ${left_read_file_base_name}"_star_alignment_stats.txt" $alignments
mv $log_out_file_name ${left_read_file_base_name}"_star_align_summary.txt"
cp ${left_read_file_base_name}"_star_align_summary.txt" $alignments

mv $left_read_file_base_name"Unmapped.out.mate1" $left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq"
mv $left_read_file_base_name"Unmapped.out.mate2" $right_read_file_base_name"_unalignedRead2AgainstTranscriptome.fq"

## remove unaligned files from bowtie after aligning them
rm $leftUnalignedFile
rm $rightUnalignedFile

lineCount="$(wc -l $left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq" | cut -d' ' -f1)"
fastqCount=$(expr $lineCount / 4)
echo $left_read_file_base_name","$fastqCount >> "../finished_bowtie_star.csv"

echo "END OF STAR FILE"