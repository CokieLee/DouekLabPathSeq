#!/bin/sh
#SBATCH -J buildSalmon
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# java
# R
##########################

codePath=$1
mergedPseudocountsFile=$2
mergedTPMFile=$3
baseName=$4
origin=$5 ## RNA, DNA or all
taxLevel=$6
PathSeqSplitOutputTableByTaxonomy_program=$7
rScriptDiv=$8
outPath=$9

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

echo "Merged pseudocounts file:"
variable_is_empty $mergedPseudocountsFile
file_exist_check $mergedPseudocountsFile
echo $mergedPseudocountsFile
echo "Merged TPM file:"
variable_is_empty $mergedTPMFile
file_exist_check $mergedTPMFile
echo $mergedTPMFile

echo "baseName:"
variable_is_empty $baseName
echo $baseName
echo "origin"
variable_is_empty $origin
echo $origin
echo "taxLevel"
variable_is_empty $taxLevel
echo $taxLevel

echo "Split output table by taxonomy program:"
variable_is_empty $PathSeqSplitOutputTableByTaxonomy_program
file_exist_check $PathSeqSplitOutputTableByTaxonomy_program
echo $PathSeqSplitOutputTableByTaxonomy_program
echo "rScriptDiv"
variable_is_empty $rScriptDiv
file_exist_check $rScriptDiv
echo $rScriptDiv

echo "outPath:"
variable_is_empty $outPath
directory_exists $outPath
echo $outPath
#####################################################

##Make diversity metric output folder
diversityOutDir=$outPath"/"$origin"_diversity/"$taxLevel
echo "diversity output directory:"
echo $diversityOutDir

mkdir $diversityOutDir
directory_exists $diversityOutDir

## TODO: fix this directory movement
startDir=$(pwd)
cd $diversityOutDir
echo "current dir: $diversityOutDir"
cp $mergedPseudocountsFile .
cp $mergedTPMFile .

mergedPseudocountsLocal=$origin"_"$taxLevel"_"$baseName"_pseudocounts.csv"
mergedTPMLocal=$origin"_"$taxLevel"_"$baseName"_tpm.csv"

echo "CHECKPOINT 2: CALL SPLIT OUTPUT BY TAXONOMY PROGRAM"

java -Xmx4G -jar $PathSeqSplitOutputTableByTaxonomy_program $mergedPseudocountsLocal k
java -Xmx4G -jar $PathSeqSplitOutputTableByTaxonomy_program $mergedTPMLocal k

function generateDiversity()
{
	local kingdom=$1

	## input file, output of merge
	file=$kingdom"_"$mergedPseudocountsLocal

	## output file names
	diversityOut=$kingdom"_"$origin"_"$taxLevel"_"$baseName"_diversity.csv"
	echo "diversityOut: $diversityOut"
	distMatrixOut=$kingdom"_"$origin"_"$taxLevel"_"$baseName"_distance.csv"
	echo "distmatrixOut: $distMatrixOut"
	distTreeOut=$kingdom"_"$origin"_"$taxLevel"_"$baseName"_distance_tree.pdf"
	echo "distTreeOut: $distTreeOut"
	distHeatMapOut=$kingdom"_"$origin"_"$taxLevel"_"$baseName"_distance_heatmap.pdf"
	echo "distHeatMapOut: $distHeatMapOut"

	if [ -f "$file" ]; then
		echo "$file exists."
		Rscript $rScriptDiv $file $diversityOut $distMatrixOut $distTreeOut $distHeatMapOut
	else 
		echo "$file does not exist."
	fi
}


echo "CHECKPOINT 3: CALL GENERATEDIVERSITY ON EUKARYOTA, BACTERIA, VIRUSES, ARCHAEA"

# kingdom=Eukaryota
# generateDiversity $kingdom $baseName

kingdom=Bacteria
generateDiversity $kingdom $baseName

# kingdom=Viruses
# generateDiversity $kingdom $baseName

# kingdom=Archaea
# generateDiversity $kingdom $baseName

cd $startDir

echo "DIVERSITY FINISHED"