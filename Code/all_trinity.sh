#!/bin/sh
#SBATCH -J all_trinity
#SBATCH --cpus-per-task=12
#SBATCH --mem=17G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

export TMPDIR=/hpcdata/scratch/
export _JAVA_OPTIONS="-Djava.io.tmpdir=/hpcdata/scratch"


module load py-setuptools
module load trinity
export PATH="/data/vrc_his/douek_lab/projects/PathSeq/programs/fastp:$PATH"

unalignedInputLeft=$1
unalignedInputRight=$2
left_read_file_base_name=$3
right_read_file_base_name=$4
origin=$5
MIN_CONTIG_LENGTH=$6
codePath=$7
outPath=$8

## Source script for directory checking function
dos2unix $codePath"dir_check.sh"
source $codePath"/dir_check.sh"

echo "CHECKPOINT 1: ALL_TRINITY INPUTS:"
echo "1. unalignedInputLeft:"
echo $unalignedInputLeft
echo "2. unalignedInputRight:"
echo $unalignedInputRight

echo "3. left_read_file_base_name:"
echo $left_read_file_base_name
echo "4. right_read_file_base_name:"
echo $right_read_file_base_name

echo "5. origin:"
echo $origin
echo "6. min_contig_length:"
echo $MIN_CONTIG_LENGTH

echo "7. codePath:"
echo $codePath
echo "8. outPath:"
echo $outPath

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

echo "CHECKPOINT 2: FASTP COMMAND:"
fastp -i $unalignedInputLeft -o $fastp_output_file_name_for_trinity_left -I $unalignedInputRight -O $fastp_output_file_name_for_trinity_right --dedup --thread 12

trinityOutDirectory="myTrinity_Origin_"$origin"_Sample_"$left_read_file_base_name
if [ -e $trinityOutDirectory ]
then
	rm -r $trinityOutDirectory
fi

Trinity_fa_out_file_name=$trinityOutDirectory".Trinity.fasta"
formatted_Trinity_fa_out_file_name="formatted_"$trinityOutDirectory".Trinity.fasta"
trinityRunOutput="trinity_out_"$MIN_CONTIG_LENGTH"_"$left_read_file_base_name".txt"
trinityRunError="trinity_err_"$left_read_file_base_name".txt"

echo "CHECKPOINT 3: TRINITY COMMAND: "

Trinity --CPU 12 --output $trinityOutDirectory --min_contig_length $MIN_CONTIG_LENGTH --seqType fq --max_memory 175G  --min_kmer_cov 1 --left $unalignedInputLeft --right $unalignedInputRight --no_version_check 1>$trinityRunOutput 2>$trinityRunError

exitCode=$?


readsPerFile=20000

if [ $exitCode -eq 0 ]
	then
		echo "Trinity sample (" $left_read_file_base_name ") successful ("$exitCode")" >> "../finished_Trinity.txt"
		module purge
		module load fastx-toolkit
		##changes the width of sequences line in a FASTA file (all nucleotide sequences appear on a single line)
		fasta_formatter -i $Trinity_fa_out_file_name > $formatted_Trinity_fa_out_file_name
	else
		echo "Trinity sample ("$left_read_file_base_name") failed ("$exitCode") retrying and STOPPING! You will need to resume the pipeline manually" > "../finished_Trinity.txt"
		##Trinity --CPU 12 --FORCE --output $trinityOutDirectory --min_contig_length $MIN_CONTIG_LENGTH --seqType fq --max_memory 175G  --min_kmer_cov 1 --left fastp_unmapped_left.fq --right fastp_unmapped_right.fq --no_version_check 1>"trinity_retry_out_"$MIN_CONTIG_LENGTH".txt" 2>trinity_err.txt

fi

## command for running next script
# qsub $full_path_to_start_filter_script $projectID $left_read_file_base_name $right_read_file_base_name $MIN_CONTIG_LENGTH $origin $readsPerFile $program_PathSeqRemoveHostForKaiju $blastDB_Mammalia $filterScript $kaiju_nodes $kaiju_fmi $kaijuScript $parseKaijuScript $PathSeqKaijuConcensusSplitter2_program $NCBI_nt_kaiju_ref_taxonomy $mergeScript $prepDiversityScript $salmonQuantScript $left_read_file  $right_read_file $PathSeqMergeQIIME2TaxAndSalmon_program $PathSeqSplitOutputTableByTaxonomy_program $palmScanScript $rScriptDiv
