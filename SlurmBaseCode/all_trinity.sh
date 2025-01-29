#!/bin/sh
#SBATCH -J all_trinity
#SBATCH --cpus-per-task=12
#SBATCH --mem=17G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# trinity
# fastx-toolkit (fasta_formatter)
# fastp

## indices ##
# bowtie primate index
##########################

codePath=$1
unalignedInputLeft=$2
unalignedInputRight=$3
left_read_file_base_name=$4
right_read_file_base_name=$5
origin=$6
MIN_CONTIG_LENGTH=$7
numTrinityThreads=$8
maxTrinityMemory=$9
outPath=${10}
scratchDir=${11}

## print input args
## check that all input args exist and are non-zero
#####################################################
echo "CHECKPOINT 1: ALL_TRINITY INPUTS:"

echo "1. codePath:"
echo $codePath
if [[ -d "$codePath" ]]; then
    echo "Directory exists: $codePath"
  else
    echo "Directory does not exist: $codePath. FAILED. QUITTING."
    exit 1
fi
## Source script for directory checking function
dos2unix $codePath"dir_check.sh"
source $codePath"/dir_check.sh"

echo "2. unalignedInputLeft:"
echo $unalignedInputLeft
file_exist_check $unalignedInputLeft
echo "3. unalignedInputRight:"
echo $unalignedInputRight
file_exist_check $unalignedInputRight

echo "4. left_read_file_base_name:"
echo $left_read_file_base_name
variable_is_empty $left_read_file_base_name
echo "5. right_read_file_base_name:"
echo $right_read_file_base_name
variable_is_empty $right_read_file_base_name

echo "6. origin:"
echo $origin
variable_is_empty $origin
echo "7. min_contig_length:"
echo $MIN_CONTIG_LENGTH
variable_is_empty $MIN_CONTIG_LENGTH
echo "8. numTrinityThreads:"
echo $numTrinityThreads
variable_is_empty $numTrinityThreads
echo "9. maxTrinityMemory:"
echo $maxTrinityMemory
variable_is_empty $maxTrinityMemory

echo "10. outPath:"
echo $outPath
directory_exists $outPath
echo "11. scratchDir:"
echo $scratchDir
directory_exists $scratchDir

export TMPDIR=$scratchDir
export _JAVA_OPTIONS="-Djava.io.tmpdir=$scratchDir"
#####################################################

# Confirm input exists
file_exist_check $unalignedInputLeft
file_exist_check $unalignedInputRight

## Make directory to store trinity output
trinity_folder_name=$origin"_trinity_output"
trinity_folder=$outPath/$trinity_folder_name
mkdir $trinity_folder

## Confirm that the directory exists
file_exist_check $trinity_folder

## Change directories into that folder
cd $trinity_folder

## Confirm that we are in that folder
correct_cur_Dir=$trinity_folder_name
dir_check  $correct_cur_Dir

fastp_output_file_name_for_trinity_left=$trinity_folder"/fastp_unalignedRead1AgainstPrimate_"$left_read_file_base_name".fq"
fastp_output_file_name_for_trinity_right=$trinity_folder"/fastp_unalignedRead1AgainstPrimate_"$right_read_file_base_name".fq"

echo "CHECKPOINT 2: FASTP COMMAND:"
fastpCmd="fastp -i $unalignedInputLeft -o $fastp_output_file_name_for_trinity_left -I $unalignedInputRight -O $fastp_output_file_name_for_trinity_right --dedup --thread 12"
echo $fastpCmd

echo "$fastpCmd" | bash

echo "fastp command complete"

trinityOutDirectory="myTrinity_Origin_"$origin"_Sample_"$left_read_file_base_name
if [ -e $trinityOutDirectory ]
then
	rm -r $trinityOutDirectory
fi

Trinity_fa_out_file_name="$trinity_folder/$trinityOutDirectory.Trinity.fasta"
formatted_Trinity_fa_out_file_name="$trinity_folder/formatted_"$trinityOutDirectory".Trinity.fasta"
trinityRunOutput="$trinity_folder/trinity_out_"$MIN_CONTIG_LENGTH"_Sample_"$left_read_file_base_name".txt"
trinityRunError="$trinity_folder/trinity_err_"$MIN_CONTIG_LENGTH"_Sample_"$left_read_file_base_name".txt"

## TODO: using $trinityoutDirectory in trinityCmd below is confusing, fix
echo "CHECKPOINT 3: TRINITY COMMAND: "

trinityCmd="Trinity --CPU $numTrinityThreads --output $trinityOutDirectory --min_contig_length $MIN_CONTIG_LENGTH \
--seqType fq --max_memory $maxTrinityMemory  --min_kmer_cov 1 --left $unalignedInputLeft --right $unalignedInputRight \
--no_version_check 1>$trinityRunOutput 2>$trinityRunError"

echo $trinityCmd

echo "$trinityCmd" | bash

exitCode=$?

readsPerFile=20000

## TODO: substitute the following error checking functions for something in dir_checks.sh.
## May need to rewrite some of those functions to take arguments such as which file we are checking.

if [ $exitCode -eq 0 ]
	then
		echo "Trinity sample ($left_read_file_base_name) successful ($exitCode)" >> "../finished_Trinity.txt"

		if [ ! -s "$Trinity_fa_out_file_name" ];
			then
				echo "Error: Input file for fasta_formatter not found: $Trinity_fa_out_file_name"
				exit 1
			else
				echo "Trinity output file present."
		fi
	else
		echo "Trinity sample ("$left_read_file_base_name") failed ("$exitCode") retrying and STOPPING! You will need to resume the pipeline manually" > "../finished_Trinity.txt"
fi

#changes the width of sequences line in a FASTA file (all nucleotide sequences appear on a single line)
		format_cmd="fasta_formatter -i $Trinity_fa_out_file_name -o $formatted_Trinity_fa_out_file_name -e"
		echo "Format command: $format_cmd"
		echo "$format_cmd" | bash > fasta_formatter_error.log 2>&1

		format_exit_code=$?
		if [ $format_exit_code -ne 0 ];
			then
				echo "Fasta_formatter failed with exit code: $format_exit_code"
				exit 1
			else
				echo "Fasta_formatter finished with successful exit code: $format_exit_code"
		fi

echo "END OF TRINITY SCRIPT"