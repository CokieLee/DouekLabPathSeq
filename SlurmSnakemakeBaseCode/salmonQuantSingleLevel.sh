#!/bin/sh
#SBATCH -J buildSalmon
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# java
# salmon

## indices ##
# kaiju indices (what reference?)
##########################

codePath=$1
left_read_file=$2
right_read_file=$3
baseName=$4
origin=$5 ## RNA, DNA or all
taxLevel=$6
salmon_index=$7
outPath=$8

## print input args
## check that all input args exist and are non-zero
#####################################################
echo "CHECKPOINT 1: SALMON INPUTS:"

echo "1. codePath:"
if [[ -d "$codePath" ]]; then
    echo "Directory exists: $codePath"
  else
    echo "Directory does not exist: $codePath. FAILED. QUITTING."
    exit 1
fi
echo $codePath

## Source script for directory checking function
dos2unix $codePath"dir_check.sh"
source $codePath"dir_check.sh"

echo "2. left read file:"
variable_is_empty $left_read_file
file_exist_check $left_read_file
echo $left_read_file
echo "3. right read file:"
variable_is_empty $right_read_file
file_exist_check $right_read_file
echo $right_read_file
echo "4. baseName:"
variable_is_empty $baseName
echo $baseName

echo "5. origin:"
variable_is_empty $origin
echo $origin
echo "6. taxLevel:"
variable_is_empty $taxLevel
echo $taxLevel

echo "7. salmon index:"
variable_is_empty $salmon_index
file_exist_check $salmon_index
echo $salmon_index
echo "8. outPath:"
variable_is_empty $outPath
directory_exists $outPath
echo $outPath
##########################################

## part1 is repetitive so let's define and use a function 
salmon_quantification() {
  my_taxonomy_level=$1
  path_to_left_read_file=$2
  path_to_right_read_file=$3
  baseName=$4
  salmon_output_dir=$5
  salmon_index=$6
  ####################################
  echo "Path to left read file inside function:"
  echo $path_to_left_read_file

  echo "Path to right read file inside function:"
  echo $path_to_right_read_file

  #Confirm that we can get to Salmon quant folder
  file_exist_check $salmon_output_dir

  ####################################
  ### The quant command quantifies transcripts, with the index provided by the -i option, -1 and -2 denote the left
  ### and right original paired reads
  ## The -l A option asks Salmon to automatically infer the library type to determine if library should be treated 
  ## as single-end or paired end.  See here for more info: https://salmon.readthedocs.io/en/latest/salmon.html
  ## The -p option specifies the number of threads to be used, in this case 1. 
  ## The validateMappings option enables selective alignment
  ## The -o option specifies the output file
  salmon quant -i $salmon_index -l A -1 $path_to_left_read_file -2 $path_to_right_read_file -p 1 --validateMappings -o $salmon_output_dir

  #Confirm that output file was created
  file_exist_check $salmon_output_dir 
  file_exist_check $salmon_output_dir"/"quant.sf

  mv $salmon_output_dir"/"quant.sf $salmon_output_dir"/"$baseName"_quant.sf"

  file_exist_check $salmon_output_dir"/"$baseName"_quant.sf"
}

echo "CHECKPOINT 2: CALL SALMON QUANTIFICATION"

salmon_output_dir=$outPath"/"$taxLevel"_quant_"$baseName
echo "salmon quant output dir:"
echo $salmon_output_dir
mkdir $salmon_output_dir
directory_exists $salmon_output_dir

salmonQuantCmd="salmon_quantification $taxLevel \
                $left_read_file \
                $right_read_file \
                $baseName \
                $salmon_output_dir \
                $salmon_index"

echo "Salmon quant call:"
echo $salmonQuantCmd
eval "$salmonQuantCmd"
process_fail_check "Salmon quantfication call failed. QUITTING."

echo "END OF SALMON"