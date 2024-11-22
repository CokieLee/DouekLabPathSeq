#!/bin/sh
#SBTACH -J starAfterBowtie
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# star
# samtools

## indices ##
# star human index
##########################

codePath=$1
leftUnalignedFile=$2
rightUnalignedFile=$3
left_read_file_base_name=$4
right_read_file_base_name=$5
outputPath=$6
hg38_starDB=$7

## print input args
## check that all input args exist and are non-zero
#####################################################
echo "STAR AFTER BOWTIE INPUTS:"

echo "codePath:"
echo $codePath
if [[ -d "$codePath" ]]; then
    echo "Directory exists: $codePath"
  else
    echo "Directory does not exist: $codePath. FAILED. QUITTING."
    exit 1
fi

## Source script for directory checking function
dos2unix $codePath"dir_check.sh"
source $codePath"dir_check.sh"

echo "leftUnalignedFile:"
echo $leftUnalignedFile
file_exist_check $leftUnalignedFile
echo "rightUnalignedFile:"
echo $rightUnalignedFile
file_exist_check $rightUnalignedFile

echo "output_base_name:"
echo $left_read_file_base_name
variable_is_empty $left_read_file_base_name
echo $right_read_file_base_name
variable_is_empty $right_read_file_base_name

echo "hg38_starDB:"
echo $hg38_starDB
file_exist_check $hg38_starDB

echo "outputPath:"
echo $outputPath
directory_exists $outputPath
#####################################################

startDir=pwd

alignStatsDir=$outputPath"alignment_stats"
starGeneratedDir=$outputPath"Generated_Data_Star_Alignment"

mkdir $alignStatsDir
mkdir $starGeneratedDir

cd $starGeneratedDir

## Confirm that we are in the folder Generated_Data_Star_Alignment
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="Generated_Data_Star_Alignment"
dir_check  $correct_cur_Dir

## TODO: remove the relative directory pathing in star cmd

echo "CHECKPOINT 1: star command:"
starCmd="STAR --genomeDir $hg38_starDB --readFilesIn $leftUnalignedFile $rightUnalignedFile --outFileNamePrefix ./$left_read_file_base_name --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt"
echo starCmd

echo "$starCmd" | bash

exitCode=$?

bam_out_file_name=${left_read_file_base_name}"Aligned.sortedByCoord.out.bam"
log_out_file_name=${left_read_file_base_name}"Log.final.out"
echo "BAM out file name:"
echo $bam_out_file_name
echo "Log.out file name:"
echo $log_out_file_name

if [ $exitCode -eq 0 ]
	then
		echo "FINISHED star command"
	else
        echo "Star command had exit code 1. If the error log star_stderr.txt is empty,\
        consider checking if memory is sufficient (human genome typically requires at least 16G for STAR). FAILED. QUITTING."
        exit 1
fi

## TODO: add checking that the output of star exists, and if possible, that is it well formed

samtools index $bam_out_file_name
samtools idxstats $bam_out_file_name > ${left_read_file_base_name}"_star_alignment_stats.txt"

## TODO: add checking that samtools index succeeded

cp ${left_read_file_base_name}"_star_alignment_stats.txt" $alignStatsDir
mv $log_out_file_name ${left_read_file_base_name}"_star_align_summary.txt"
cp ${left_read_file_base_name}"_star_align_summary.txt" $alignStatsDir

mv $left_read_file_base_name"Unmapped.out.mate1" $left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq"
mv $left_read_file_base_name"Unmapped.out.mate2" $right_read_file_base_name"_unalignedRead2AgainstTranscriptome.fq"

## remove unaligned files from bowtie after aligning them
# rm $leftUnalignedFile
# rm $rightUnalignedFile

lineCount="$(wc -l $left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq" | cut -d' ' -f1)"
fastqCount=$(expr $lineCount / 4)
echo $left_read_file_base_name","$fastqCount >> "../finished_bowtie_star.csv"


cd $starDir
echo "END OF STAR FILE"