#!/bin/sh
#SBATCH -J PathSeqSubmitter
#SBATCH --mem=2G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

readsPerFile=20000 ## Trinity reads per file (used when splitting output)
## Define absolute paths to reference sets

inputpath=/home/parkercol/PathSeq/Tests/TestInput/
origin=RNA
minContigLen=300

## configurable parameters
codepath=/home/parkercol/PathSeq/SlurmSnakemakeBaseCode/
outputpath=/data/parkercol/PathseqOutput/
scratchDir=/lscratch/

# reference sets and databases
bowtieERCCIndex=/data/parkercol/PathseqDependencyData/bowtieERCCIndex/ERCC92   ## local 9M
bowtieUnmaskedGenomeIndex=/fdb/htgts/genomes/bowtie2_indexes/hg19/hg19   ## Biowulf hosted

bowtiePrimatePaths=/data/parkercol/PathseqDependencyData/bowtiePrimateIndex/primates  ## local 30G
blastDB_Mammalia=/data/parkercol/blastDBMammalia/Mammalia.fa   ## local 65G
kaiju_nodes=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20220512/kaiju/nodes.dmp    ## TODO=update local. 159M
kaiju_fmi=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20220512/kaiju/nr_euk/kaiju_db_nr_euk.fmi   ##TODO=update local. 149G

hg38_starDB=/fdb/STAR_indices/2.7.11b/UCSC/hg38/ref    ## Biowulf hosted

## external programs
# all individual executable files have been copied to biowulf home directory
program_picard=/home/parkercol/PathseqDependencyPrograms/picard.jar
program_RemoveHostForKaiju=/home/parkercol/PathseqDependencyPrograms/PathSeqRemoveHostForKaiju.jar
program_PathSeqKaijuConcensusSplitter2=/home/parkercol/PathseqDependencyPrograms/PathSeqKaijuConcensusSplitter2.jar
program_PathSeqMergeQIIME2TaxAndSalmon=/home/parkercol/PathseqDependencyPrograms/PathSeqMergeQIIME2TaxAndSalmon.jar
program_PathSeqSplitOutputTableByTaxonomy=/home/parkercol/PathseqDependencyPrograms/PathSeqSplitOutputTableByTaxonomy.jar
program_prodigal=/home/parkercol/PathseqDependencyPrograms/prodigal/bin
program_kaiju=/home/parkercol/PathseqDependencyPrograms/kaiju-v1.9.0-linux-x86_64-static
program_fastp=/home/parkercol/PathseqDependencyPrograms/fastp
program_palmscan=/home/parkercol/PathseqDependencyPrograms/palmscan-main/bin

## Names of scripts to be called (should be stored in scripts folder)
bowtieScript="bowtieUnmaskedGenome.sh"
starScript="starAfterBowtie.sh"
primateScript="bowtiePrimate.sh"
trinityScript="all_trinity.sh"
split_filter_submit_script="start_filter_host_kaiju.sh"
filterScript="array_filter_host.sh"
kaijuScript="protein_kaiju.sh"
parseKaijuScript="parse_protein_kaiju_build_salmon.sh"
mergeScript="mergeTaxonomyAndQuantificationSingleLevel.sh"
prepDiversityScript="prepForDiversity.sh"
salmonQuantScript="salmonQuantSingleLevel.sh"
palmScanScript="palm_scan.sh"

## temporary
sampleNameLeft="582760806_A5_S33_R1_001_short_test"
sampleNameRight="582760806_A5_S33_R2_001_short_test"
right_read_file="$inputpath/582760806_A5_S33_R1_001_short_test.fastq.gz"
left_read_file="$inputpath/582760806_A5_S33_R2_001_short_test.fastq.gz"

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

## Confirm that we are in the scripts folder
#correct_cur_Dir="test" #Unit testing of dir check function (comment out for regular use)
correct_cur_Dir="scripts" ## Regular operation
dir_check  $correct_cur_Dir 

## Check if script files with those names actually exist and are not empty (size >0) in the scripts folder
#file_exist_check "unit_test_file" ## Uncomment if unit testing
file_exist_check $bowtieScript
file_exist_check $starScript
file_exist_check $primateScript
file_exist_check $trinityScript
file_exist_check $split_filter_submit_script
file_exist_check $filterScript
file_exist_check $kaijuScript
file_exist_check $parseKaijuScript
file_exist_check $mergeScript
file_exist_check $prepDiversityScript
file_exist_check $salmonQuantScript
file_exist_check $palmScanScript
file_exist_check $rScriptDiv

## Sources the functions for extracting the unique IDs for each read from the original read file name
## The "./" makes the file "get_unique_sample_ID_names.sh" executable
dos2unix get_unique_sample_ID_names.sh
source ./get_unique_sample_ID_names.sh

###########################################

###########################################
## submit first alignment pass and capture the array jobID
bowtie_cmd="sbatch --parsable $bowtieScript $codepath $sampleNameLeft $sampleNameRight $left_read_file $right_read_file $outputpath $bowtieERCCIndex $bowtieUnmaskedGenomeIndex $program_picard"
echo "Bowtie cmd: bowtie_cmd"
bowtieJobID=$(sbatch --parsable $bowtieScript $codepath $sampleNameLeft $sampleNameRight $left_read_file $right_read_file $outputpath $bowtieERCCIndex $bowtieUnmaskedGenomeIndex $program_picard)
###########################################
## get the short jobID from array jobID
echo "holdID is "
echo $bowtieJobID

echo "Ready for star Alignment"
###########################################
## submit second alignment pass, holding for first pass to finish and capture the array jobID
starJobID=$(sbatch --parsable --dependency=afterok:$bowtieJobID $starScript $projectID $left_read_file_base_name  $right_read_file_base_name $hg38_starDB| cut -d' ' -f3)
###########################################
## get the short jobID from array jobID
echo "holdID is "
echo $starJobID

echo "Ready for bowtie Primate"
###########################################
## submit third alignment pass, holding for second pass to finish and capture the array jobID
bowtie2JobID=$(sbatch --parsable --dependency=afterok:$starJobID $primateScript $projectID $left_read_file_base_name  $right_read_file_base_name $bowtiePrimateIndex | cut -d' ' -f3)
# ###########################################
## get the short jobID from array jobID
echo "holdID is "
echo $bowtie2JobID

echo "Ready for Trinity"
###########################################
## submit the trinity script, holding for the all alignment passes to finish 
trinityJobID=$(sbatch --parsable --dependency=afterok:$bowtie2JobID $trinityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $MIN_CONTIG_LENGTH $split_filter_submit_script $program_PathSeqRemoveHostForKaiju $blastDB_Mammalia $filterScript $kaiju_nodes $kaiju_fmi $kaijuScript $parseKaijuScript $PathSeqKaijuConcensusSplitter2_program $NCBI_nt_kaiju_ref_taxonomy $mergeScript $prepDiversityScript $salmonQuantScript $left_read_file  $right_read_file $PathSeqMergeQIIME2TaxAndSalmon_program $PathSeqSplitOutputTableByTaxonomy_program $palmScanScript $rScriptDiv | cut -d' ' -f3)

echo "holdID is "
echo $trinityJobID


# ## Sources the functions for splitting trinity output files
# ## The "./" makes the file "split_trinity_output.sh" executable
# dos2unix split_trinity_output.sh
# source ./split_trinity_output.sh 


# fileCount=0


# split_trinity_out_files  $projectID $left_read_file_base_name $right_read_file_base_name $$MIN_CONTIG_LENGTH $origin $readsPerFile $fileCount
    

# echo "Function output"
# echo $fileCount




#holdID5=$(sbatch -hold_jid $holdID $split_filter_submit_script $projectID $left_read_file_base_name $right_read_file_base_name $MIN_CONTIG_LENGTH  $readsPerFile $origin | cut -d' ' -f3)


# sbatch -hold_jid $holdID $parseKaijuScript $projectID $file $origin
