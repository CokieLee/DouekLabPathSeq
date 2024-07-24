#!/bin/sh
#SBATCH -J PathSeqSubmitter
#SBATCH --mem=2G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

projectID=$1
left_read_file=$2
right_read_file=$3
origin=$4
MIN_CONTIG_LENGTH=$5
readsPerFile=20000 ## Trinity reads per file (used when splitting output)
## Define absolute paths to reference sets

### Bowtie ERCC reference index
bowtieERCCIndex=/data/vrc_his/douek_lab/reference_sets/ERCC/ERCC92
### Bowtie unmasked human genome reference index
bowtieUnmaskedGenomeIndex="/data/vrc_his//douek_lab/reference_sets/hg19/unmasked/genome"
## starDB reference databases 
hg38_starDB="/data/vrc_his/douek_lab/reference_sets/hg38/Sequence/STAR/"
#hg_38_gtf="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg38/Annotation/Genes/genes.gtf"
#hg38_referenceFasta="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg38/Sequence/WholeGenomeFasta/genome.fa"

### Bowtie primate reference index
bowtiePrimateIndex=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20210302/primates

## Paths to programs
picard=/data/vrc_his/douek_lab/projects/PathSeq/programs/picard/2.18.14/picard.jar
program_PathSeqRemoveHostForKaiju=/data/vrc_his/douek_lab/projects/PathSeq/programs/PathSeqRemoveHostForKaiju/dist/PathSeqRemoveHostForKaiju.jar

PathSeqKaijuConcensusSplitter2_program=/data/vrc_his/douek_lab/projects/PathSeq/programs/PathSeqKaijuConcensusSplitter2/dist/PathSeqKaijuConcensusSplitter2.jar

PathSeqMergeQIIME2TaxAndSalmon_program=/data/vrc_his/douek_lab/projects/PathSeq/programs/PathSeqMergeQIIME2TaxAndSalmon/dist/PathSeqMergeQIIME2TaxAndSalmon.jar

PathSeqSplitOutputTableByTaxonomy_program=/data/vrc_his/douek_lab/projects/PathSeq/programs/PathSeqSplitOutputTableByTaxonomy/dist/PathSeqSplitOutputTableByTaxonomy.jar

## Blast reference databases
blastDB_Mammalia=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20210716/Mammalia.fa

## Files for kaiju reference database  (.dmp is a memory dump file format, fmi= functional 
## mockup interface, format used for model simulations)

kaiju_nodes=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20220512/kaiju/nodes.dmp
kaiju_fmi=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20220512/kaiju/nr_euk/kaiju_db_nr_euk.fmi

NCBI_nt_kaiju_ref_taxonomy=/data/vrc_his/douek_lab/reference_sets/NCBI_nt/20220512/kaiju/nr_euk/qiime2_formatted_taxonomy.tab

## Rscript for calculating diversity metrics
rScriptDiv="pathDiv.R"

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


left_read_file_base_name=""
right_read_file_base_name=""
get_left_read_unique_ID  $projectID $left_read_file $left_read_file_base_name
get_right_read_unique_ID  $projectID $right_read_file $right_read_file_base_name

echo "Function output"
echo $left_read_file_base_name
echo $right_read_file_base_name

###########################################

###########################################
## submit first alignment pass and capture the array jobID
holdID=$(qsub $bowtieScript $projectID $left_read_file $right_read_file $left_read_file_base_name  $right_read_file_base_name $bowtieERCCIndex $bowtieUnmaskedGenomeIndex $picard | cut -d' ' -f3)
###########################################
## get the short jobID from array jobID
echo $holdID

OLD_IFS=$IFS
IFS="."
newArray=($holdID)
holdID=${newArray[0]}
IFS=$OLD_IFS

echo "holdID is "
echo $holdID
echo "Ready for star Alignment"
###########################################
## submit second alignment pass, holding for first pass to finish and capture the array jobID
holdID2=$(qsub -hold_jid $holdID $starScript $projectID $left_read_file_base_name  $right_read_file_base_name $hg38_starDB| cut -d' ' -f3)
###########################################
## get the short jobID from array jobID
echo $holdID2

OLD_IFS=$IFS
IFS="."
newArray=($holdID2)
holdID=${newArray[0]}
IFS=$OLD_IFS

echo "holdID is "
echo $holdID

echo "Ready for bowtie Primate"
###########################################
## submit third alignment pass, holding for second pass to finish and capture the array jobID
holdID3=$(qsub -hold_jid $holdID $primateScript $projectID $left_read_file_base_name  $right_read_file_base_name $bowtiePrimateIndex | cut -d' ' -f3)
# ###########################################
## get the short jobID from array jobID
echo $holdID3

OLD_IFS=$IFS
IFS="."
newArray=($holdID3)
holdID=${newArray[0]}
IFS=$OLD_IFS

echo $holdID
###########################################
## submit the trinity script, holding for the all alignment passes to finish 
holdID4=$(qsub -hold_jid $holdID $trinityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $MIN_CONTIG_LENGTH $split_filter_submit_script $program_PathSeqRemoveHostForKaiju $blastDB_Mammalia $filterScript $kaiju_nodes $kaiju_fmi $kaijuScript $parseKaijuScript $PathSeqKaijuConcensusSplitter2_program $NCBI_nt_kaiju_ref_taxonomy $mergeScript $prepDiversityScript $salmonQuantScript $left_read_file  $right_read_file $PathSeqMergeQIIME2TaxAndSalmon_program $PathSeqSplitOutputTableByTaxonomy_program $palmScanScript $rScriptDiv | cut -d' ' -f3)

echo $holdID4

OLD_IFS=$IFS
IFS="."
newArray=($holdID4)
holdID=${newArray[0]}
IFS=$OLD_IFS

echo $holdID


# ## Sources the functions for splitting trinity output files
# ## The "./" makes the file "split_trinity_output.sh" executable
# dos2unix split_trinity_output.sh
# source ./split_trinity_output.sh 


# fileCount=0


# split_trinity_out_files  $projectID $left_read_file_base_name $right_read_file_base_name $$MIN_CONTIG_LENGTH $origin $readsPerFile $fileCount
    

# echo "Function output"
# echo $fileCount




#holdID5=$(qsub -hold_jid $holdID $split_filter_submit_script $projectID $left_read_file_base_name $right_read_file_base_name $MIN_CONTIG_LENGTH  $readsPerFile $origin | cut -d' ' -f3)


# qsub -hold_jid $holdID $parseKaijuScript $projectID $file $origin
