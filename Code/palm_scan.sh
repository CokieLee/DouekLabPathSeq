#!/bin/sh
#$ -N palm_scan
#$ -S /bin/bash
#$ -M rahul.subramanian@nih.gov
#$ -m n
#$ -l h_vmem=20G
#$ -l quick
#$ -cwd

export TMPDIR=/hpcdata/scratch/
export _JAVA_OPTIONS="-Djava.io.tmpdir=/hpcdata/scratch"
export PATH="/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/programs/palmscan-main/bin:$PATH"

projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4 ## RNA, DNA or all

outDir="../../palmscan/"
origin_sample_unique_id_tag=$origin"_Sample_"$left_read_file_base_name

outFile=$outDir$origin_sample_unique_id_tag"_palmscan.fa"
#cd "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

## Confirm that we are in the scripts directory
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

##Change directories to project folder
cd "../"

##Confirm that we are in project folder
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


mkdir palmscan
file_exist_check palmscan

cd $origin"_trinity_output/salmon/"

#cd "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/"$origin"_trinity/salmon/"

## Arguments:  (source: excerpted from Palmscan  documentation palmscan -help: ) 
## ppout = Amino acid palmprint sequences
## search_pp = Input file (can be nt or aa sequence)
## report = Human readable alignment format
## rt = Report reverse transcriptases
## rdrp = Report RdRPs (RNA dependent RNA polymerase)
## fevout = Field-equals-value format
palmscan -threads 1 -search_pp unclassified_sequences.fa -rt -rdrp -ppout $origin_sample_unique_id_tag"_pp.fa" -report $origin_sample_unique_id_tag"_pp.txt" -fevout $origin_sample_unique_id_tag"_pp.fev"

grep ">" $origin_sample_unique_id_tag"_pp.fa" > $origin_sample_unique_id_tag"_palmscan_headers.txt"

for line in $(cat $origin_sample_unique_id_tag"_palmscan_headers.txt" ); 
do
	query=$line
	grep -A 1 $query unclassified_sequences.fa >> $outFile
done

## Return to scripts folder at end of script
cd ../../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir