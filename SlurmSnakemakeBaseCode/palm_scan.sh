#!/bin/sh
#SBATCH -J buildSalmon
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# palmscan
##########################

codePath=$1
unclassified_sequences=$2
program_palmscan=$3
scratchDir=$4
baseName=$5
origin=$6 ## RNA, DNA or all
outPath=$7

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
dos2unix $codePath/dir_check.sh
source $codePath/dir_check.sh

echo "2. Consensus splitter (build salmon) unclassified sequences:"
variable_is_empty $unclassified_sequences
file_exist_check $unclassified_sequences
echo $unclassified_sequences

echo "3. program palmscan:"
variable_is_empty $program_palmscan
directory_exists $program_palmscan
echo $program_palmscan
echo "4. scratch directory:"
variable_is_empty $scratchDir
directory_exists $scratchDir
echo $scratchDir

echo "5. baseName:"
variable_is_empty $baseName
echo $baseName
echo "6. origin:"
variable_is_empty $origin
echo $origin

echo "7. outPath:"
variable_is_empty $outPath
directory_exists $outPath
echo $outPath

export TMPDIR=$scratchDir
export _JAVA_OPTIONS="-Djava.io.tmpdir=$scratchDir"
export PATH="$program_palmscan:$PATH"
#####################################################

outDir=$outPath"/palmscan/"
mkdir $outDir
directory_exists $outDir

echo "CHECKPOINT 2: RUN PALMSCAN"

outFile=$outDir$baseName"_palmscan.fa"

## Arguments:  (source: excerpted from Palmscan  documentation palmscan -help: ) 
## ppout = Amino acid palmprint sequences
## search_pp = Input file (can be nt or aa sequence)
## report = Human readable alignment format
## rt = Report reverse transcriptases
## rdrp = Report RdRPs (RNA dependent RNA polymerase)
## fevout = Field-equals-value format
palmscan -threads 1 -search_pp $unclassified_sequences -rt -rdrp -ppout $outDir$baseName"_pp.fa" -report $outDir$baseName"_pp.txt" -fevout $outDir$baseName"_pp.fev"

grep ">" $outDir$baseName"_pp.fa" > $outDir$baseName"_palmscan_headers.txt"

for line in $(cat $outDir$baseName"_palmscan_headers.txt" ); 
do
	query=$line
	grep -A 1 $query $unclassified_sequences >> $outFile
done

echo "FINISHED PALMSCAN"