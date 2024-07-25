#!/bin/sh
#SBATCH -J mergeTaxonomyAndQuantificationSingleLevel
#SBATCH --mem=25G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## load Locus modules
module load javafx


COUNTER=$SGE_TASK_ID
projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4 ## RNA, DNA or all
mytaxLevel=$5
PathSeqMergeQIIME2TaxAndSalmon_program=$6

salmonQuantBase="/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/"$origin"_salmon_quant/"

echo $projectID
echo $left_read_file_base_name
echo $right_read_file_base_name
echo $origin
echo $taxLevel

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

salmon_quant_folder_name=$origin"_salmon_quant"
## Confirm that folder exists a
file_exist_check $salmon_quant_folder_name
## Change into salmon quant folder 
cd $salmon_quant_folder_name

origin_sample_unique_id_tag=$origin"_Sample_"$left_read_file_base_name

## part1 is repetitive so let's define and use a function 
salmon_merge() {
	
	taxLevel=$1
	origin_sample_unique_id_tag=$2

	##Confirm folder exists
	file_exist_check $taxLevel"_quant_"$origin_sample_unique_id_tag

	##Change into it
	cd $taxLevel"_quant_"$origin_sample_unique_id_tag

	##Verify current directory
	correct_cur_Dir=$taxLevel"_quant_"$origin_sample_unique_id_tag
	dir_check  $correct_cur_Dir
	
	taxFile="../../"$origin"_trinity_output/salmon/"$taxLevel"_table.csv"
	## Confirm that the file exists
	file_exist_check $taxFile

	finalFile=$origin_sample_unique_id_tag"_quant.sf"
	
	## Confirm that the file exists
	file_exist_check $finalFile

	## Merge salmon output (which uses Trinity contig labels and Salmon counts) with Kaiju output (which has taxonomic classifications
	## for each Trinity contig)
	java -Xmx20G -jar $PathSeqMergeQIIME2TaxAndSalmon_program $finalFile $origin_sample_unique_id_tag $taxFile $origin"_"$taxLevel"_"$origin_sample_unique_id_tag 2> $taxLevel"_merge_stderr.txt"
	
	##Confirm that output files were created
	file_exist_check $origin"_"$taxLevel"_"$origin_sample_unique_id_tag"_pseudocounts.csv"
	file_exist_check $origin"_"$taxLevel"_"$origin_sample_unique_id_tag"_tpm.csv"

	mv $origin"_"$taxLevel"_"$origin_sample_unique_id_tag"_pseudocounts.csv" ../
	mv $origin"_"$taxLevel"_"$origin_sample_unique_id_tag"_tpm.csv" ../
}

salmon_merge $mytaxLevel $origin_sample_unique_id_tag

## Return to scripts folder at end of script
cd ../../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir