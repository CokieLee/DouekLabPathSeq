#!/bin/sh
#$ -N start_filter_host_kaiju
#$ -S /bin/bash
#$ -M rahul.subramanian@nih.gov
#$ -m be
#$ -l quick
#$ -cwd
projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
MIN_CONTIG_LENGTH=$4
origin=$5  ## RNA, DNA or all
readsPerFile=$6
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

cur_Dir=$(basename $(pwd))
full_dir=$(pwd)

echo "Current basename work dir is "
echo $cur_Dir


echo "Full working dir  is "
echo $full_dir

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

## Confirm that we are in the scripts directory
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

## Confirm correct script file inputs
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

trinity_output_folder=$origin"_trinity_output"
## Confirm that the trinity output folder exists
file_exist_check $trinity_output_folder

## Change directories into that folder
cd $trinity_output_folder

## Confirm that we are in that folder
correct_cur_Dir=$trinity_output_folder
dir_check  $correct_cur_Dir


## Define sample specific trinity split output folder
trinitySplitDirectory="splitTrinity_"$left_read_file_base_name

## If folder already exists, delete it
if [ -e $trinitySplitDirectory ]
then
rm -r $trinitySplitDirectory
fi

## Make the folder and change directories into it
mkdir  $trinitySplitDirectory
cd $trinitySplitDirectory

    
trinityOutDirectory="myTrinity_Origin_"$origin"_Sample_"$left_read_file_base_name

formatted_Trinity_fa_out_file_name="formatted_"$trinityOutDirectory".Trinity.fasta"


trinityFile=$formatted_Trinity_fa_out_file_name

## The number of lines per file that we would like will be twice the number of reads we want to 
## include in the file (one line in the formatted fasta file
## is the header, and the other line is the sequence)
linesPerFile=$(expr $readsPerFile \* 2)
## Count the number of lines in the formatted trinity FASTA output file that we currently have
lineCount="$(cat ../$trinityFile | wc -l | cut -d' ' -f1)"

## Add the number of lines that we have in the file currently to the number of lines per file we would like, and subtract 1
newNumerator=$(expr $lineCount + $linesPerFile - 1)
fileCount=$(expr $newNumerator / $linesPerFile)

## Reasoning: Division in unix shell scripts with the "/" operator rounds down. 
## For example, if lineCount=11 and $linesPerFile=2, we would like to make 6 files. 
## However, if we just did fileCount=$(expr $lineCount / $linesPerFile), we would obtain only 5 files.
## Adding $linesPerFile to $lineCount in the numerator partially corrects this isssue 
## newNumerator=$(expr $lineCount + $linesPerFile ), and then fileCount=$(expr $newNumerator / $linesPerFile)
## In this case, the new numerator is 12, so we will obtain 6 as desired. Incrementing the numerator by $linesPerFile
## ensures that most edge cases are dealt with. However, the only problem now is when $lineCount is a direct multiple 
## of $linesPerFile, for example, if lineCount=10 and $linesPerFile=2. We would like to make 5 files in this scenario, 
## but because we are incrementing by $linesPerFile, we obtain 6 as the number of files. Subtracting 1 from the numerator
## before dividing resolves this issue. 

## Estimate the number of files by summing the total number of lines in the file  dividing the number of lines in the filoe

trinity_prefix="trinity_"$origin"_Sample_"$left_read_file_base_name"_split_"

## Print out (via cat) the formatted Trinity output file, then pipe contents and split them into
## files of line count $linesPerFile via the split command, giving all files the prefix trinity
## and a numeric suffix.
## The default format is split [options] filename prefix,
## or -l linenumber instead of options
## The -l $linesPerFile dentoes that we want to split the Trinity file into files where each 
## has only $linesPerFile number of lines. 
## The -a 4 means generate suffixes of length 4
## The argument "--numeric-suffixes=1" means use numeric suffixes starting at 1
## The - trinity means that the prefix "trinity" is used for all files.

cat ../$trinityFile | split -a 4 --numeric-suffixes=1 -l $linesPerFile - $trinity_prefix

echo $fileCount

## Return to scripts folder at end of script
cd ../../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir
echo "Line 141"
echo $filterScript
holdID=$(qsub -t 1-$fileCount -tc 100 $filterScript $projectID $left_read_file_base_name $right_read_file_base_name $origin $program_PathSeqRemoveHostForKaiju $blastDB_Mammalia  | cut -d' ' -f3)

OLD_IFS=$IFS
IFS="."
newArray=($holdID)
holdID=${newArray[0]}
IFS=$OLD_IFS

echo "line 151"
echo $kaijuScript
holdID=$(qsub -hold_jid $holdID $kaijuScript $projectID  $left_read_file_base_name $right_read_file_base_name $origin $fileCount $kaiju_nodes $kaiju_fmi | cut -d' ' -f3)
qsub -hold_jid $holdID $parseKaijuScript $projectID $left_read_file_base_name $right_read_file_base_name $origin $kaiju_nodes $PathSeqKaijuConcensusSplitter2_program $NCBI_nt_kaiju_ref_taxonomy $mergeScript $prepDiversityScript $salmonQuantScript $left_read_file  $right_read_file $PathSeqMergeQIIME2TaxAndSalmon_program $PathSeqSplitOutputTableByTaxonomy_program $palmScanScript $rScriptDiv

