#!/bin/sh
#$ -N bowtieUnmaskedGenome
#$ -S /bin/bash
#$ -m n
#$ -pe threaded 12
#$ -l quick
#$ -cwd

module load ncurses/5.9-goolf-1.7.20
module load Bowtie2/2.3.3.1-py27pl5.22
module load FASTX-Toolkit/0.0.14-goolf-1.7.20
#module load BEDTools/2.25.0-goolf-1.7.20
module load SAMtools
#module load Java/1.8.0_92
module load Java/1.8.0_45
module load cURL/7.43.0-goolf-1.7.20
module load ncurses/5.9-goolf-1.7.20

#picard=/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/programs/picard/2.18.14/picard.jar

COUNTER=$SGE_TASK_ID
projectID=$1

# must be full paths:
left_read_file=$2
right_read_file=$3

left_read_file_base_name=$4
right_read_file_base_name=$5
outPathLeft=$6
outPathRight=$7

# debug
echo "out path right: "
echo $outPathRight

scriptsPath=$8

bowtieERCCIndex=$9
bowtieUnmaskedGenomeIndex=${10}

#debug
echo "code path: "
echo $scriptsPath
echo "bowtie ercc index: "
echo $bowtieERCCIndex
echo "bowtie unmasked genome index: "
echo $bowtieUnmaskedGenomeIndex

picard=${11}

# get starting directory so we can return to it
startDir=$(pwd)

## Source script for directory checking function
dos2unix $scriptsPath/dir_check.sh
source $scriptsPath/dir_check.sh 


echo "Line 28"
echo $left_read_file

rel_path_to_proj_dir="../"

#Todo: change all of the following variable name in the script
# because rel implies its relative, but it is not anymore
left_read_file_rel_path_from_script=$left_read_file
right_read_file_rel_path_from_script=$right_read_file



echo 'Line 38'
echo $left_read_file_rel_path_from_script


##########################################
#####

##########################################

generated_Data_dir=$outputPath"/Generated_Data/"
outputDir_unmaskedGenome=$outputPath"/unmasked_genome_alignment_rates/"
erccoutputDir=$outputPath"/ERCC_alignment_rates/"
erccAlignmentCountsDir=$outputPath"/ERCC_alignment_counts/"
insertSizeOutputDir=$outputPath"/insert_size_metrics"

outSam_ERCC="genome_alignment_ERCC_"$left_read_file_base_name".sam"
outSam_Unmasked_Genome="genome_alignment_unmasked_genome_"$left_read_file_base_name".sam"
#outSam="genome_alignment.sam"

cur_Dir=$(basename $(pwd))
echo $cur_Dir

#Confirm that sample folder exists
file_exist_check $outPathLeft

cd $outPathLeft

mkdir Generated_Data_1st_Bowtie_Alignment_ERCC

mkdir Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome
mkdir unmasked_genome_alignment_rates
mkdir insert_size_metrics
mkdir ERCC_alignment_rates
mkdir ERCC_alignment_counts

cd Generated_Data_1st_Bowtie_Alignment_ERCC/
#Confirm directory change was succesful
correct_cur_Dir="Generated_Data_1st_Bowtie_Alignment_ERCC"
dir_check  $correct_cur_Dir
echo $left_read_file
echo $right_read_file
echo "line 89"
##########################################
rm *.xml
##rm *.fq
rm *.sam
rm *.bam*
rm no_ERCC1.fq.gz
rm no_ERCC2.fq.gz
rm "unalignedRead1AgainstGenome_"${left_read_file_base_name}".fq"
rm "unalignedRead2AgainstGenome_"${right_read_file_base_name}".fq"
rm "unalignedRead1AgainstERCC_"$left_read_file_base_name".fq"
rm "unalignedRead2AgainstERCC_"$right_read_file_base_name".fq"
rm "readindex_index_"$index".fq.gz"
rm left*
rm right*
rm *.cxb
##rm *.txt
rm *.csv
##rm *.tab
##rm *tab.gz
rm *readcount
##########################################

pwd


bowtieAlignRate_ERCC="../ERCC_alignment_rates/ERCC_alignment_rate_"$left_read_file_base_name".txt"


# Todo: examine the following variable names
left_read_file_rel_path_from_ERCC_bowtie_folder=$left_read_file_rel_path_from_script
right_read_file_rel_path_from_ERCC_bowtie_folder=$right_read_file_rel_path_from_script

## Confirm path to read file works
file_exist_check $left_read_file_rel_path_from_ERCC_bowtie_folder
file_exist_check $right_read_file_rel_path_from_ERCC_bowtie_folder

# debug
echo "ALIGNMENT:"
echo "path to read file: "
echo $left_read_file_rel_path_from_ERCC_bowtie_folder
echo "output path: "
echo $outSam_ERCC
echo "alignment rate path: "
echo $bowtieAlignRate_ERCC
echo " "

# Perform first alignment step
bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieERCCIndex -1 $left_read_file_rel_path_from_ERCC_bowtie_folder -2 $right_read_file_rel_path_from_ERCC_bowtie_folder 1>$outSam_ERCC 2>$bowtieAlignRate_ERCC

bowtie_cmd_str="bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieERCCIndex -1 $left_read_file_rel_path_from_ERCC_bowtie_folder -2 $right_read_file_rel_path_from_ERCC_bowtie_folder 1>$outSam_ERCC 2>$bowtieAlignRate_ERCC"

echo "bowtie test command"
echo $bowtie_cmd_str
echo "Done printing bowtie test command"

outERCCbam="ERCCAlignments_"$left_read_file_base_name".bam"
#outERCCbam="ERCCAlignments.bam"
sortedERCCbam="sorted.ERCCAlignments_"$left_read_file_base_name".bam"
cur_Dir=$(basename $(pwd))

samtools_version=$(samtools --version)

# create indexes and human readable output for still unaligned reads
if [ -e $outSam_ERCC ]
then
    echo $cur_Dir >&1
    echo $outSam_ERCC >&1
    echo $outSam_ERCC >&1
    echo "Test line 245"
    echo $cur_Dir >&2
    echo $outSam_ERCC >&2
    echo $samtools_version 
    echo "Before samtools view call" >&2
    samtools view -@ 12 -b -F 4 $outSam_ERCC > $outERCCbam
    echo "After samtools view call" >&2

    # debug
    echo $sortedERCCbam
    samtools sort -@ 12 -m 7G $outERCCbam 1> $sortedERCCbam
    # debug
    echo $sortedERCCbam
    samtools index -@ 12 $sortedERCCbam
    # debug
    echo "current directory before calling samtools idxstats:"
    echo $cur_Dir
    echo $erccAlignmentCountsDir"ERCC_alignment_counts_rate_"$left_read_file_base_name".txt"
    samtools idxstats $sortedERCCbam 1> $erccAlignmentCountsDir"ERCC_alignment_counts_rate_"$left_read_file_base_name".txt"

    samtools view -@ 12 -f 4 $outSam_ERCC | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstERCC_"$left_read_file_base_name".fq"
    samtools view -@ 12 -f 4 $outSam_ERCC | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead2AgainstERCC_"$right_read_file_base_name".fq"

else
    printf '%s\n' "Sam file does not exist/not found " >&1
    printf '%s\n' "Should be in Generated_Data_1st_Bowtie_Alignment_ERCC " >&1
    printf '%s\n' "But instead is in  " >&1
    echo $cur_Dir >&1
    printf '%s\n' "Sam file does not exist/not found " >&2
    printf '%s\n' "Should be in Generated_Data_1st_Bowtie_Alignment_ERCC " >&2
    printf '%s\n' "But instead is in  " >&2
    echo $cur_Dir >&2
    exit 1

fi

# debug
echo "SECOND ALIGNMENT STEP: "

#Confirm directory change was succesful
correct_cur_Dir="Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome"
cd ../Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome

#bowtieIndex="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg19/unmasked/genome"
bowtieAlignRate_UnmaskedGenome="genome_alignment_rate_"$left_read_file_base_name".txt"

bowtieERCCDirPath="../Generated_Data_1st_Bowtie_Alignment_ERCC/"

ERRC_left_read_out_path_from_2nd_bowtie_folder=$bowtieERCCDirPath"unalignedRead1AgainstERCC_"$left_read_file_base_name".fq"
ERRC_right_read_out_path_from_2nd_bowtie_folder=$bowtieERCCDirPath"unalignedRead2AgainstERCC_"$right_read_file_base_name".fq"

## Confirm paths to ERCC ouptuts work (files exist)
file_exist_check $ERRC_left_read_out_path_from_2nd_bowtie_folder
file_exist_check $ERRC_left_read_out_path_from_2nd_bowtie_folder

# debug
echo $ERRC_left_read_out_path_from_2nd_bowtie_folder
echo $bowtieAlignRate_UnmaskedGenome

# 2nd alignment step
bowtie2 -p 12 --no-mixed --no-discordant -x $bowtieUnmaskedGenomeIndex -1 $ERRC_left_read_out_path_from_2nd_bowtie_folder -2 $ERRC_right_read_out_path_from_2nd_bowtie_folder 1>$outSam_Unmasked_Genome 2>$bowtieAlignRate_UnmaskedGenome
cp $bowtieAlignRate_UnmaskedGenome $outputDir_unmaskedGenome


samtools view -@ 12 -f 4 $outSam_Unmasked_Genome | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstGenome_"$left_read_file_base_name".fq"
samtools view -@ 12 -f 4 $outSam_Unmasked_Genome | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead2AgainstGenome_"$right_read_file_base_name".fq"

module load R/3.4.3-goolf-1.7.20
tempGenomeAlignments_filename="tempGenomeAlignments_"$left_read_file_base_name".bam"

samtools view -b $outSam_Unmasked_Genome > $tempGenomeAlignments_filename

tempSortedAlignments_filename="tempSortedAlignments_"$left_read_file_base_name".bam"
samtools sort -m 20G $tempGenomeAlignments_filename> $tempSortedAlignments_filename
echo "Reached line 145"
picard_stdout_file_name="picard_stdout_"$left_read_file_base_name".txt" 
picard_stderr_file_name="picard_stderr_"$left_read_file_base_name".txt" 
insert_size_metrics_file_name="insert_size_metrics_"$left_read_file_base_name".txt" 
insert_Hist_pdf_file_name="insert_size_histogram_"$left_read_file_base_name".pdf" 

# debug
echo "jar file name: "
echo $picard
echo "picard collectInsertSizeMetrics input: "
echo $tempSortedAlignments_filename

java -Xmx6G -jar $picard CollectInsertSizeMetrics I=$tempSortedAlignments_filename O=$insert_size_metrics_file_name H=$insert_Hist_pdf_file_name M=0.5 1>$picard_stdout_file_name 2>$picard_stderr_file_name
OUT=$?
echo "Reached line 148"

# debug
echo "insert_Hist_pdf_file_name exists: "
if [ -e $insert_Hist_pdf_file_name ]
then
    echo "true"
else
    echo "false"
fi


if [ -e $insert_Hist_pdf_file_name ]
then
	echo $left_read_file_base_name" Everything successful ("$OUT") and deleting intermediate files" >> "../finished_bowtieUnmaskedGenome.txt"
	rm $outSam_Unmasked_Genome
	rm $tempSortedAlignments_filename
	rm $tempGenomeAlignments_filename
	#rm unmapped.genome.bam
	#rm temp.unaligned.bam
	rm *ERCCAlignments*
	##rm unalignedRead1AgainstERCC.fq
	##rm unalignedRead2AgainstERCC.fq
	mv $insert_size_metrics_file_name "../insert_size_metrics/"$insert_size_metrics_file_name
	mv $insert_Hist_pdf_file_name "../insert_size_metrics/"$insert_Hist_pdf_file_name
	#cp $left_read_file_base_name"_insert_size_metrics.txt" $insertSizeOutputDir
	#cp $left_read_file_base_name"_insert_size_histogram.pdf" $insertSizeOutputDir
else
    echo $left_read_file_base_name" Something went wrong ("$OUT"), keeping first intermediate alignment file" >> "../finished_bowtieUnmaskedGenome.txt"
    rm $tempSortedAlignments_filename
	rm $tempGenomeAlignments_filename
   	#rm unmapped.genome.bam
   	#rm tempSortedAlignments.bam
fi

# return to starting directory
cd $startDir

echo "Reached end of bowtieUnmasked.sh"
