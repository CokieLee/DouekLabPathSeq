#!/bin/sh
#SBATCH -J salmonQuantSingleLevel
#SBATCH --mem=100G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov


## load Locus modules
module load salmon

COUNTER=$SGE_TASK_ID
projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4 ## RNA, DNA or all
taxLevel=$5
left_read_file=$6
right_read_file=$7

echo $projectID
echo $left_read_file_base_name
echo $right_read_file_base_name
echo $origin
echo $taxLevel
echo $left_read_file
echo $right_read_file

##########################################

## Get relative paths to original base reads
rel_path_to_proj_dir="../"

left_read_file_rel_path_from_salmon=$left_read_file
right_read_file_rel_path_from_salmon=$right_read_file

echo "Read path relative to salmon folder"

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

## Make salmon quant folder if it doesn't already exist
salmon_quant_folder_name=$origin"_salmon_quant"
mkdir $salmon_quant_folder_name

## Confirm that the directory exists (will cd into it later on)
file_exist_check $salmon_quant_folder_name

trinity_output_folder=$origin"_trinity_output"
## Confirm that the trinity output folder exists
file_exist_check $trinity_output_folder

## Change directories into that folder
cd $trinity_output_folder

## Confirm that we are in that folder
correct_cur_Dir=$trinity_output_folder
dir_check  $correct_cur_Dir

## Confirm that salmon folder exists
file_exist_check salmon

## Change into salmon folder
cd salmon

## Confirm that we are in that folder
correct_cur_Dir=salmon
dir_check  $correct_cur_Dir

origin_sample_unique_id_tag=$origin"_Sample_"$left_read_file_base_name

## part1 is repetitive so let's define and use a function 
salmon_quantification() {
  my_taxonomy_level=$1
  path_to_left_read_file=$2
  path_to_right_read_file=$3
  origin_sample_unique_id_tag=$4
  salmon_quant_folder_name_in_func=$5
  salmon_index=$my_taxonomy_level"_salmon/"
  ####################################
  ## Just in case you need to resubmit
  rm -r $my_taxonomy_level"_quant"

  echo "Path to left read file inside function"
  echo $path_to_left_read_file

  echo "Path to right read file inside function"
  echo $path_to_right_read_file

  #Confirm that we can get to Salmon quant folder from where we are
  ## Confirm that the directory exists (will cd into it later on)
  path_to_salmon_quant_folder="../../"$salmon_quant_folder_name_in_func"/"
  file_exist_check $path_to_salmon_quant_folder 

  salmon_output_file_name_and_path=$path_to_salmon_quant_folder$my_taxonomy_level"_quant_"$origin_sample_unique_id_tag

  ####################################
  ### The quant command quantifies transcripts, with the index provided by the -i option, -1 and -2 denote the left
  ### and right original paired reads
  ## The -l A option asks Salmon to automatically infer the library type to determine if library should be treated 
  ## as single-end or paired end.  See here for more info: https://salmon.readthedocs.io/en/latest/salmon.html
  ## The -p option specifies the number of threads to be used, in this case 1. 
  ## The validateMappings option enables selective alignment
  ## The -o option specifies the output file
  salmon quant -i $salmon_index -l A -1 $path_to_left_read_file -2 $path_to_right_read_file -p 1 --validateMappings -o $salmon_output_file_name_and_path

  #Confirm that output file was created
  file_exist_check $salmon_output_file_name_and_path 
  file_exist_check $salmon_output_file_name_and_path"/"quant.sf

  mv $salmon_output_file_name_and_path"/"quant.sf $salmon_output_file_name_and_path"/"$origin_sample_unique_id_tag"_quant.sf"

  file_exist_check $salmon_output_file_name_and_path"/"$origin_sample_unique_id_tag"_quant.sf"

  ## Return to scripts folder at end of script
  cd ../../../scripts/

  ## Confirm that we are in fact in the scripts folder
  cur_Dir=$(basename $(pwd))
  #echo $cur_Dir
  correct_cur_Dir="scripts"
  dir_check  $correct_cur_Dir
  
}

salmon_quantification $taxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag $salmon_quant_folder_name
