#!/bin/sh
#$ -N concat_kaiju_build_salmon
#$ -S /bin/bash
#$ -M rahul.subramanian@nih.gov
#$ -m be
#$ -l h_vmem=100G
#$ -cwd

module load java

projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4 ## RNA, DNA or all
kaiju_nodes=$5
PathSeqKaijuConcensusSplitter2_program=$6
NCBI_nt_kaiju_ref_taxonomy=$7
mergeScript=$8
prepDiversityScript=$9
salmonQuantScript=${10}
left_read_file=${11}
right_read_file=${12}
PathSeqMergeQIIME2TaxAndSalmon_program=${13}
PathSeqSplitOutputTableByTaxonomy_program=${14}
palmScanScript=${15}
rScriptDiv=${16}

#baseDir="/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/"$origin"_trinity/"



echo "left_read_file_base_name"
echo $left_read_file_base_name

echo "origin"
echo $origin

echo "merge script"
echo $mergeScript

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

## Confirm that we are in the scripts directory
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

echo "salmonQuantScript"
echo $salmonQuantScript

file_exist_check $mergeScript
file_exist_check $prepDiversityScript
file_exist_check $salmonQuantScript
file_exist_check $palmScanScript
file_exist_check $rScriptDiv

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


rm -r salmon
mkdir salmon
cd salmon

origin_sample_unique_id_tag=$origin"_Sample_"$left_read_file_base_name

sorted_protein_kaiju_output_file_tab="../sorted_protein_kaiju_"$origin_sample_unique_id_tag".tab"

formatted_non_host_proteins_nucleotide_sequences_file_fa="../formatted_non_host_proteins_nucleotide_"$origin_sample_unique_id_tag".fa"

formatted_non_host_proteins_translations_output_file_faa="../formatted_non_host_proteins_"$origin_sample_unique_id_tag".faa"



java -Xmx90G -jar $PathSeqKaijuConcensusSplitter2_program $sorted_protein_kaiju_output_file_tab $formatted_non_host_proteins_nucleotide_sequences_file_fa $formatted_non_host_proteins_translations_output_file_faa $NCBI_nt_kaiju_ref_taxonomy $kaiju_nodes 2>stderr.txt

module purge
module load salmon

## Get indexes for salmon using the salmon indexer (no pre-defined reference, denoted by the -t flag)
salmon index -t kingdom_sequences.fa -i kingdom_salmon
salmon index -t phylum_sequences.fa -i phylum_salmon
salmon index -t class_sequences.fa -i class_salmon
salmon index -t order_sequences.fa -i order_salmon
salmon index -t family_sequences.fa -i family_salmon
salmon index -t genus_sequences.fa -i genus_salmon
salmon index -t species_sequences.fa -i species_salmon

module purge
module load uge


## Get relative paths to original base reads
rel_path_to_proj_dir_from_salmon="../../../"

left_read_file_rel_path_from_salmon=${rel_path_to_proj_dir_from_salmon}$left_read_file
right_read_file_rel_path_from_salmon=${rel_path_to_proj_dir_from_salmon}$right_read_file

## Confirm that we can find the files
file_exist_check $left_read_file_rel_path_from_salmon
file_exist_check $right_read_file_rel_path_from_salmon

# salmon_quantification() {
#   my_taxonomy_level=$1
#   path_to_left_read_file=$2
#   path_to_right_read_file=$3
#   origin_sample_unique_id_tag=$4
#   salmon_index=$my_taxonomy_level"_salmon/"
#   ####################################
#   ## Just in case you need to resubmit
#   rm -r $my_taxonomy_level"_quant"

#   #Confirm that we can get to Salmon quant folder from where we are
#   ## Confirm that the directory exists (will cd into it later on)
#   path_to_salmon_quant_folder="../../"$salmon_quant_folder_name"/"
#   file_exist_check $path_to_salmon_quant_folder 

#   salmon_output_file_name_and_path=$path_to_salmon_quant_folder$my_taxonomy_level"_quant_"$origin_sample_unique_id_tag

#   ####################################
#   ### The quant command quantifies transcripts, with the index provided by the -i option, -1 and -2 denote the left
#   ### and right original paired reads
#   ## The -l A option asks Salmon to automatically infer the library type to determine if library should be treated 
#   ## as single-end or paired end.  See here for more info: https://salmon.readthedocs.io/en/latest/salmon.html
#   ## The -p option specifies the number of threads to be used, in this case 1. 
#   ## The validateMappings option enables selective alignment
#   ## The -o option specifies the output file
#   salmon quant -i $salmon_index -l A -1 $left -2 $right -p 1 --validateMappings -o $salmon_output_file_name_and_path

#   #Confirm that output file was created
#   file_exist_check $salmon_output_file_name_and_path 
  
# }

# mytaxLevel=kingdom
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

# mytaxLevel=phylum
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

# mytaxLevel=class
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

# mytaxLevel=order
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

# mytaxLevel=family
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

# mytaxLevel=genus
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

# mytaxLevel=species
# salmon_quantification $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon $origin_sample_unique_id_tag

### now kick off the salmon quantification scripts as an array job


## Return to scripts folder at end of script
cd ../../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

file_exist_check $mergeScript
file_exist_check $prepDiversityScript
file_exist_check $salmonQuantScript

mytaxLevel=kingdom
holdID_1=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
echo "holdID post 207"
echo $holdID_1
#tempJobsHoldID=$(getShortHoldID $holdID)
#echo "temp Jobs Hold ID"
#echo $tempJobsHoldID

hold_jid_trunc=$holdID_1
hold_jid=$holdID_1
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_1_trunc=$hold_jid_trunc

echo "holdID_1_trunc"
echo $holdID_1_trunc


echo "merge script before call"
echo $mergeScript
holdID_2=$(qsub -hold_jid $holdID_1_trunc $mergeScript $projectID $left_read_file_base_name $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)

hold_jid_trunc=$holdID_2
hold_jid=$holdID_2
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_2_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_2_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv

mytaxLevel=phylum
holdID_4=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name  $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
hold_jid_trunc=$holdID_4
hold_jid=$holdID_4
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_4_trunc=$hold_jid_trunc

holdID_5=$(qsub -hold_jid $holdID_4_trunc $mergeScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)

hold_jid_trunc=$holdID_5
hold_jid=$holdID_5
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_5_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_5_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv


mytaxLevel=class
holdID_7=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name  $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
hold_jid_trunc=$holdID_7
hold_jid=$holdID_7
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_7_trunc=$hold_jid_trunc

holdID_8=$(qsub -hold_jid $holdID_7_trunc $mergeScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)
hold_jid_trunc=$holdID_8
hold_jid=$holdID_8
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_8_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_8_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv


mytaxLevel=order
holdID_10=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name  $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
hold_jid_trunc=$holdID_10
hold_jid=$holdID_10
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_10_trunc=$hold_jid_trunc

holdID_11=$(qsub -hold_jid $holdID_10_trunc $mergeScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)
hold_jid_trunc=$holdID_11
hold_jid=$holdID_11
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_11_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_11_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv


mytaxLevel=family
holdID_13=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name  $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
hold_jid_trunc=$holdID_13
hold_jid=$holdID_13
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_13_trunc=$hold_jid_trunc

holdID_14=$(qsub -hold_jid $holdID_13_trunc $mergeScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)
hold_jid_trunc=$holdID_14
hold_jid=$holdID_14
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_14_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_14_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv


mytaxLevel=genus
holdID_16=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name  $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
hold_jid_trunc=$holdID_16
hold_jid=$holdID_16
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_16_trunc=$hold_jid_trunc


holdID_17=$(qsub -hold_jid $holdID_16_trunc $mergeScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)
hold_jid_trunc=$holdID_17
hold_jid=$holdID_17
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_17_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_17_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv


mytaxLevel=species
holdID_19=$(qsub -t 1-1 -tc 100 $salmonQuantScript $projectID $left_read_file_base_name  $right_read_file_base_name  $origin $mytaxLevel $left_read_file_rel_path_from_salmon $right_read_file_rel_path_from_salmon | cut -d' ' -f3)
hold_jid_trunc=$holdID_19
hold_jid=$holdID_19
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_19_trunc=$hold_jid_trunc


holdID_20=$(qsub -hold_jid $holdID_19_trunc $mergeScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqMergeQIIME2TaxAndSalmon_program | cut -d' ' -f3)
hold_jid_trunc=$holdID_20
hold_jid=$holdID_20
get_trunc_hold_jid  $hold_jid $hold_jid_trunc
holdID_20_trunc=$hold_jid_trunc

qsub -hold_jid $holdID_20_trunc $prepDiversityScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin $mytaxLevel $PathSeqSplitOutputTableByTaxonomy_program $rScriptDiv

if [[ $origin == *"RNA"* ]]; then
	qsub $palmScanScript $projectID $left_read_file_base_name  $right_read_file_base_name $origin
	echo "Submitting palmscan"
fi



# ###########################################
# ## count the number of lines in the sub_sheet, which corresponds to the number of jobs to submit
# lineCount="$(wc -l $file | cut -d' ' -f1)"
# echo $lineCount ## just to also have the information in the stdout on Locus
# ## just have to fix it because it counts return lines
# lineCount=`expr $lineCount + 1`
# echo $lineCount ## just to also have the information in the stdout on Locus
# ###########################################
# function getShortHoldID()
# {
# 	local myLongHoldID=$1
# 	OLD_IFS=$IFS
# 	IFS="."
# 	newArray=($myLongHoldID)
# 	myShortHoldID=${newArray[0]}
# 	IFS=$OLD_IFS

# 	echo $myShortHoldID
# }

# mytaxLevel=kingdom
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel $left_read_file  $right_read_file | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# mytaxLevel=phylum
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# mytaxLevel=class
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# mytaxLevel=order
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# mytaxLevel=family
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# mytaxLevel=genus
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# mytaxLevel=species
# holdID=$(qsub -t 1-$lineCount -tc 100 $salmonQuantScript $file $origin $mytaxLevel | cut -d' ' -f3)
# tempJobsHoldID=$(getShortHoldID $holdID)
# holdID=$(qsub -hold_jid $tempJobsHoldID $mergeScript $projectID $file $origin $mytaxLevel | cut -d' ' -f3)
# qsub -hold_jid $holdID $prepDiversityScript $projectID $origin $mytaxLevel

# if [[ $origin == *"RNA"* ]]; then
# 	#qsub "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/scripts/palm_scan.sh" $projectID $origin
# 	#echo "Submitting palmscan"
# fi

