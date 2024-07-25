#!/bin/sh
#$ -N bowtiePrimate
#$ -S /bin/bash
#$ -M cokie.parker@nih.gov
#$ -m be
#$ -pe threaded 12
#$ -l quick
#$ -cwd
#$ -j y

module load bowtie2
module load samtools

COUNTER=$SGE_TASK_ID

projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
bowtiePrimateIndex=$4

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

echo $left_read_file_base_name
echo $right_read_file_base_name
echo "line 25"



##########################################
#####
left="unalignedRead1AgainstTranscriptome.fq"
right="unalignedRead2AgainstTranscriptome.fq"
##########################################
outputDir_Primate="../primate_alignment_rates/"

outSam_Primate="primate_alignment_bowtie_"$left_read_file_base_name".sam"

##Confirm that we are  in scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

##Change directories to project folder
cd "../"

##Confirm that we are in project folder
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir=$projectID
dir_check  $correct_cur_Dir


sample_folder_name="Sample_"$left_read_file_base_name

## Confirm that sample folder exists
file_exist_check $sample_folder_name

## Change into sample folder
cd $sample_folder_name

## Confirm that we are in sample folder
correct_cur_Dir=$sample_folder_name
dir_check  $correct_cur_Dir

## Make folder to store primate alignment rates
mkdir primate_alignment_rates

##Change directories into that folder
cd primate_alignment_rates

## Confirm that we are in primate_alignment_rates folder
correct_cur_Dir="primate_alignment_rates"
dir_check  $correct_cur_Dir

pwd

starDir="../Generated_Data_Star_Alignment/"
pwd
rel_path_to_star_left_out_from_bowtiePrimate_folder=$starDir$left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq"
rel_path_to_star_right_out_from_bowtiePrimate_folder=$starDir$right_read_file_base_name"_unalignedRead2AgainstTranscriptome.fq"

##Confirm that the relative path to read files works (that the files exist)
file_exist_check $rel_path_to_star_left_out_from_bowtiePrimate_folder
file_exist_check $rel_path_to_star_right_out_from_bowtiePrimate_folder

## Do bowtie alignment to primate reference genome
bowtieAlignRate_Primate="primate_alignment_rate_"$left_read_file_base_name".txt"
bowtie2 -p 12 --no-mixed --no-discordant -x $bowtiePrimateIndex -1 $rel_path_to_star_left_out_from_bowtiePrimate_folder -2 $rel_path_to_star_right_out_from_bowtiePrimate_folder 1>$outSam_Primate 2>$bowtieAlignRate_Primate
#cp $bowtieAlignRate_Primate $outputDir_Primate

#From SAM file, view with 12 cores, select only unmapped reads (-f 4 option), select lines that are not
## part of the header (do not start with @), select only even or odd lines for read 1 and read 2 
## respectively. 

samtools view -@ 12 -f 4 $outSam_Primate | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstPrimate_"$left_read_file_base_name".fq"
samtools view -@ 12 -f 4 $outSam_Primate | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstPrimate_"$right_read_file_base_name".fq"

## Return to scripts folder at end of script
cd ../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir
