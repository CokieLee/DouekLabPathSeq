#!/bin/sh
#SBATCH -J bowtiePrimate
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

module load bowtie2
module load samtools

unalignedLeftFile=$1
unalignedRightFile=$2
outPath=$3
scriptsPath=$4
left_read_file_base_name=$5
right_read_file_base_name=$6
bowtiePrimateIndex=$7

##############################
# print input paths
echo "CHECKPOINT 1: BOWTIE PRIMATE INPUTS"
echo "1: unalignedLeftFile: "
echo $unalignedLeftFile
echo "2: unalignedRightFile:"
echo $unalignedRightFile

echo "3: outPath: "
echo $outPath
echo "4: scriptsPath: "
echo $scriptsPath

echo "5: left_read_file_base_name: "
echo $left_read_file_base_name
echo "6: right_read_file_base_name: "
echo $right_read_file_base_name

echo "7: bowtiePrimateIndex: "
echo $bowtiePrimateIndex
##############################

dos2unix $scriptsPath"dir_check.sh"
source $scriptsPath"dir_check.sh"

$startDir=pwd

## Confirm that output folder exists
file_exist_check $outPath
##Confirm that the relative path to read files works (that the files exist)
file_exist_check $unalignedLeftFile
file_exist_check $unalignedRightFile

outputDir_Primate=$outPath"/primate_alignment_rates/"
mkdir $outputDir_Primate

cd $outputDir_Primate

## Confirm that we are in primate_alignment_rates folder
correct_cur_Dir="primate_alignment_rates"
dir_check  $correct_cur_Dir

echo "Current Directory:"
pwd

echo "CHECKPOINT 2: RUN BOWTIE COMMAND"
## Do bowtie alignment to primate reference genome
bowtieAlignRate_Primate="primate_alignment_rate_"$left_read_file_base_name".txt"
outSam_Primate="primate_alignment_bowtie_"$left_read_file_base_name".sam"
bowtie2 -p 12 --no-mixed --no-discordant -x $bowtiePrimateIndex -1 $unalignedLeftFile -2 $unalignedRightFile 1>$outSam_Primate 2>$bowtieAlignRate_Primate

#From SAM file, view with 12 cores, select only unmapped reads (-f 4 option), select lines that are not
## part of the header (do not start with @), select only even or odd lines for read 1 and read 2 
## respectively. 

samtools view -@ 12 -f 4 $outSam_Primate | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstPrimate_"$left_read_file_base_name".fq"
samtools view -@ 12 -f 4 $outSam_Primate | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstPrimate_"$right_read_file_base_name".fq"

## Return to scripts folder at end of script
cd $starDir

echo "BOWTIE PRIMATE COMPLETE"