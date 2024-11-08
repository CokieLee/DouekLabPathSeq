#!/bin/sh
#SBATCH -J palm_scan
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

program_palmscan=$1
scratchDir=$2
unclassified_sequences=$3
salmonDataDir=$4
origin=$5 ## RNA, DNA or all
outPath=$6
codePath=$7

export TMPDIR=/hpcdata/scratch/
export _JAVA_OPTIONS="-Djava.io.tmpdir=/hpcdata/scratch"
export PATH="$program_palmscan:$PATH"

## Source script for directory checking function
dos2unix $codePath/dir_check.sh
source $codePath/dir_check.sh

outDir=$outPath"palmscan/"
mkdir $outDir
file_exist_check $outDir

origin_sample_unique_id_tag=$salmonDataDir/$origin"_Sample_"$left_read_file_base_name
outFile=$outDir$origin_sample_unique_id_tag"_palmscan.fa"

## Arguments:  (source: excerpted from Palmscan  documentation palmscan -help: ) 
## ppout = Amino acid palmprint sequences
## search_pp = Input file (can be nt or aa sequence)
## report = Human readable alignment format
## rt = Report reverse transcriptases
## rdrp = Report RdRPs (RNA dependent RNA polymerase)
## fevout = Field-equals-value format
palmscan -threads 1 -search_pp $unclassified_sequences -rt -rdrp -ppout $origin_sample_unique_id_tag"_pp.fa" -report $origin_sample_unique_id_tag"_pp.txt" -fevout $origin_sample_unique_id_tag"_pp.fev"

grep ">" $origin_sample_unique_id_tag"_pp.fa" > $origin_sample_unique_id_tag"_palmscan_headers.txt"

for line in $(cat $origin_sample_unique_id_tag"_palmscan_headers.txt" ); 
do
	query=$line
	grep -A 1 $query unclassified_sequences.fa >> $outFile
done