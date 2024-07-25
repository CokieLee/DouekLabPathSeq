#!/bin/sh
#SBATCH -J all_trinity
#SBATCH --cpus-per-task=12
#SBATCH --mem=16.5G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

export TMPDIR=/hpcdata/scratch/
export _JAVA_OPTIONS="-Djava.io.tmpdir=/hpcdata/scratch"

module load trinity
export PATH="/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/programs/fastp:$PATH"

projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4
MIN_CONTIG_LENGTH=$5
split_filter_submit_script=$6
program_PathSeqRemoveHostForKaiju=$7
blastDB_Mammalia=$8
filterScript=$9
kaiju_nodes=${10}
kaiju_fmi=${11}
kaijuScript=${12}
parseKaijuScript=${13}
PathSeqKaijuConcensusSplitter2_program=${14}
NCBI_nt_kaiju_ref_taxonomy=${15}
mergeScript=${16}
prepDiversityScript=${17}
salmonQuantScript=${18}
left_read_file=${19}
right_read_file=${20}
PathSeqMergeQIIME2TaxAndSalmon_program=${21}
PathSeqSplitOutputTableByTaxonomy_program=${22}
palmScanScript=${23}
rScriptDiv=${24}

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

echo $left_read_file_base_name
echo $right_read_file_base_name

##Confirm that we are  in scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

## Confirm correct script file inputs
file_exist_check $split_filter_submit_script
file_exist_check $filterScript
file_exist_check $kaijuScript
file_exist_check $parseKaijuScript
file_exist_check $mergeScript
file_exist_check $prepDiversityScript
file_exist_check $salmonQuantScript
file_exist_check $palmScanScript
file_exist_check $rScriptDiv

##Change directories to project folder
cd "../"

##Confirm that we are in project folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
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

## Make directory to store trinity output


trinity_folder_name=$origin"_trinity_output"
mkdir $trinity_folder_name

## Confirm that the directory exists
file_exist_check $trinity_folder_name

## Change directories into that folder
cd $trinity_folder_name

## Confirm that we are in that folder
correct_cur_Dir=$trinity_folder_name
dir_check  $correct_cur_Dir

# Get relative path to read files
bowtiePrimateDir="../primate_alignment_rates/"
rel_path_to_bowtie_Primate_left_out_from_trinity_folder=$bowtiePrimateDir"unalignedRead1AgainstPrimate_"$left_read_file_base_name".fq"
rel_path_to_bowtie_Primate_right_out_from_trinity_folder=$bowtiePrimateDir"unalignedRead1AgainstPrimate_"$right_read_file_base_name".fq"

##Confirm that the relative path to read files works (that the files exist)
file_exist_check $rel_path_to_bowtie_Primate_left_out_from_trinity_folder
file_exist_check $rel_path_to_bowtie_Primate_right_out_from_trinity_folder

fastp_output_file_name_for_trinity_left="fastp_unalignedRead1AgainstPrimate_"$left_read_file_base_name".fq"
fastp_output_file_name_for_trinity_right="fastp_unalignedRead1AgainstPrimate_"$right_read_file_base_name".fq"

fastp -i $rel_path_to_bowtie_Primate_left_out_from_trinity_folder -o $fastp_output_file_name_for_trinity_left -I $rel_path_to_bowtie_Primate_right_out_from_trinity_folder -O $fastp_output_file_name_for_trinity_right --dedup --thread 12

trinityOutDirectory="myTrinity_Origin_"$origin"_Sample_"$left_read_file_base_name
if [ -e $trinityOutDirectory ]
then
	rm -r $trinityOutDirectory
fi

Trinity_fa_out_file_name=$trinityOutDirectory".Trinity.fasta"

formatted_Trinity_fa_out_file_name="formatted_"$trinityOutDirectory".Trinity.fasta"

Trinity --CPU 12 --output $trinityOutDirectory --min_contig_length $MIN_CONTIG_LENGTH --seqType fq --max_memory 175G  --min_kmer_cov 1 --left $fastp_output_file_name_for_trinity_left --right $fastp_output_file_name_for_trinity_right --no_version_check 1>"trinity_out_"$MIN_CONTIG_LENGTH"_"$left_read_file_base_name".txt" 2>"trinity_err_"$left_read_file_base_name".txt"

exitCode=$?


readsPerFile=20000

if [ $exitCode -eq 0 ]
	then
		echo "Trinity sample (" $left_read_file_base_name ") successful ("$exitCode")" >> "../finished_Trinity.txt"
		module purge
		module load fastx-toolkit
		##changes the width of sequences line in a FASTA file (all nucleotide sequences appear on a single line)
		fasta_formatter -i $Trinity_fa_out_file_name > $formatted_Trinity_fa_out_file_name

		

		## Return to scripts folder at end of script
		cd ../../scripts/

		## Confirm that we are in fact in the scripts folder
		cur_Dir=$(basename $(pwd))
		#echo $cur_Dir
		correct_cur_Dir="scripts"
		dir_check  $correct_cur_Dir

		echo "Current directory is"
		echo $cur_Dir

		#start_filter_script="start_filter_host_kaiju.sh"


		# module purge
		# module load uge

		full_dir=$(pwd)

		full_path_to_start_filter_script=$full_dir"/"$split_filter_submit_script

		sbatch $full_path_to_start_filter_script $projectID $left_read_file_base_name $right_read_file_base_name $MIN_CONTIG_LENGTH $origin $readsPerFile $program_PathSeqRemoveHostForKaiju $blastDB_Mammalia $filterScript $kaiju_nodes $kaiju_fmi $kaijuScript $parseKaijuScript $PathSeqKaijuConcensusSplitter2_program $NCBI_nt_kaiju_ref_taxonomy $mergeScript $prepDiversityScript $salmonQuantScript $left_read_file  $right_read_file $PathSeqMergeQIIME2TaxAndSalmon_program $PathSeqSplitOutputTableByTaxonomy_program $palmScanScript $rScriptDiv

		#sbatch "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/scripts/start_filter_host_kaiju.sh" $projectID $left_read_file_base_name $right_read_file_base_name $MIN_CONTIG_LENGTH $origin $readsPerFile $program_PathSeqRemoveHostForKaiju $blastDB_Mammalia


		#holdID5=$(sbatch $start_filter_script $projectID $left_read_file_base_name $right_read_file_base_name $MIN_CONTIG_LENGTH  $readsPerFile $origin | cut -d' ' -f3)
		
		#echo $holdID5
	else
		echo "Trinity sample ("$left_read_file_base_name") failed ("$exitCode") retrying and STOPPING! You will need to resume the pipeline manually" > "../finished_Trinity.txt"
		##Trinity --CPU 12 --FORCE --output $trinityOutDirectory --min_contig_length $MIN_CONTIG_LENGTH --seqType fq --max_memory 175G  --min_kmer_cov 1 --left fastp_unmapped_left.fq --right fastp_unmapped_right.fq --no_version_check 1>"trinity_retry_out_"$MIN_CONTIG_LENGTH".txt" 2>trinity_err.txt

fi

