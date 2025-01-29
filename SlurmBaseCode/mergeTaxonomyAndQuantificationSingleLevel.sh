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
##########################

codePath=$1
taxInputFile=$2		##("$origin"_trinity_output/salmon/"$taxLevel"_table.csv")
salmonQuantInputFile=$3		##($baseName_quant.sf)
baseName=$4
origin=$5 ## RNA, DNA or all
taxLevel=$6
PathSeqMergeQIIME2TaxAndSalmon_program=$7
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

echo "2. taxonomy input file:"
variable_is_empty $taxInputFile
file_exist_check $taxInputFile
echo $taxInputFile
echo "3. salmon quantification input file:"
variable_is_empty $salmonQuantInputFile
file_exist_check $salmonQuantInputFile
echo $salmonQuantInputFile

echo "4. baseName:"
variable_is_empty $baseName
echo $baseName
echo "5. origin:"
variable_is_empty $origin
echo $origin
echo "6. taxLevel:"
variable_is_empty $taxLevel
echo $taxLevel

echo "7. pathseq merge program:"
variable_is_empty $PathSeqMergeQIIME2TaxAndSalmon_program
file_exist_check $PathSeqMergeQIIME2TaxAndSalmon_program
echo $PathSeqMergeQIIME2TaxAndSalmon_program
echo "8. outPath:"
variable_is_empty $outPath
directory_exists $outPath
echo $outPath
##########################################


## part1 is repetitive so let's define and use a function 
salmon_merge() {
	taxFile=$1
	salmonQuantFile=$2
	taxLevel=$3
	baseName=$4
	outPath=$5

	outputBaseName=$origin"_"$taxLevel"_"$baseName

	## Merge salmon output (which uses Trinity contig labels and Salmon counts) with Kaiju output (which has taxonomic classifications
	## for each Trinity contig)
	java -Xmx20G -jar $PathSeqMergeQIIME2TaxAndSalmon_program $salmonQuantFile $baseName $taxFile $outputBaseName 2> $outPath"/"$taxLevel"_merge_stderr.txt"
	
	##Confirm that output files were created
	file_exist_check $origin"_"$taxLevel"_"$baseName"_pseudocounts.csv"
	file_exist_check $origin"_"$taxLevel"_"$baseName"_tpm.csv"

	mv $origin"_"$taxLevel"_"$baseName"_pseudocounts.csv" $outPath
	mv $origin"_"$taxLevel"_"$baseName"_tpm.csv" $outPath
}

echo "CHECKPOINT 2: CALL MERGE FUNCTION"

mergeOutputDir=$outPath"/"$origin"_merge_TaxAndQuant/"

mergeFuncCmd="salmon_merge $taxInputFile $salmonQuantInputFile $taxLevel $baseName $mergeOutputDir"
echo $mergeFuncCmd
eval $mergeFuncCmd
process_fail_check

echo "merge completed successfully" > $mergeOutputDir"finished_mergeTaxAndQuant.txt"

echo "FINISHED MERGE SCRIPT"