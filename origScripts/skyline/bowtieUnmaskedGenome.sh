#!/bin/sh
#SBATCH -J bowtieUnmaskedGenome
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

module load ncurses
module load bowtie2
module load fastx-toolkit
#module load BEDTools
module load samtools
module load javafx
module load curl

#picard=/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/programs/picard/2.18.14/picard.jar

COUNTER=$SGE_TASK_ID
projectID=$1
left_read_file=$2
right_read_file=$3
left_read_file_base_name=$4
right_read_file_base_name=$5
bowtieERCCIndex=$6
bowtieUnmaskedGenomeIndex=$7
picard=$8


## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 


echo "Line 28"
echo $left_read_file

rel_path_to_proj_dir="../"

left_read_file_rel_path_from_script=${rel_path_to_proj_dir}$left_read_file
right_read_file_rel_path_from_script=${rel_path_to_proj_dir}$right_read_file



echo 'Line 38'
echo $left_read_file_rel_path_from_script


##########################################
#####

##########################################

generated_Data_dir="../Generated_Data/"
outputDir_unmaskedGenome="../unmasked_genome_alignment_rates/"
erccoutputDir="../ERCC_alignment_rates/"
erccAlignmentCountsDir="../ERCC_alignment_counts/"
insertSizeOutputDir="../insert_size_metrics"

outSam_ERCC="genome_alignment_ERCC_"$left_read_file_base_name".sam"
outSam_Unmasked_Genome="genome_alignment_unmasked_genome_"$left_read_file_base_name".sam"
#outSam="genome_alignment.sam"

cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir



cd ../

cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir=$projectID
dir_check  $correct_cur_Dir

sample_folder_name="Sample_"$left_read_file_base_name
mkdir $sample_folder_name

#Confirm that sample folder exists
file_exist_check $sample_folder_name

cd $sample_folder_name

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

test_bowtie_cmd="bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieERCCIndex -1 $left_read_file_rel_path_from_script -2 $right_read_file_rel_path_from_script 1>$outSam_ERCC 2>$bowtieAlignRate_ERCC"

left_read_file_rel_path_from_ERCC_bowtie_folder="../"$left_read_file_rel_path_from_script
right_read_file_rel_path_from_ERCC_bowtie_folder="../"$right_read_file_rel_path_from_script

## Confirm path to read file works
file_exist_check $left_read_file_rel_path_from_ERCC_bowtie_folder
file_exist_check $right_read_file_rel_path_from_ERCC_bowtie_folder

echo "path to read file"
echo $left_read_file_rel_path_from_gen_data_folder
bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieERCCIndex -1 $left_read_file_rel_path_from_ERCC_bowtie_folder -2 $right_read_file_rel_path_from_ERCC_bowtie_folder 1>$outSam_ERCC 2>$bowtieAlignRate_ERCC

echo "bowtie test command"
echo $test_bowtie_cmd
echo "Done printing bowtie test command"
#cp $bowtieAlignRate_ERCC $erccoutputDir

outERCCbam="ERCCAlignments_"$left_read_file_base_name".bam"
#outERCCbam="ERCCAlignments.bam"
sortedERCCbam="sorted.ERCCAlignments_"$left_read_file_base_name".bam"
cur_Dir=$(basename $(pwd))

samtools_version=$(samtools --version)

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

    samtools sort -@ 12 -m 7G $outERCCbam 1> $sortedERCCbam
    samtools index -@ 12 $sortedERCCbam
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
java -Xmx6G -jar $picard CollectInsertSizeMetrics I=$tempSortedAlignments_filename O=$insert_size_metrics_file_name H=$insert_Hist_pdf_file_name M=0.5 1>$picard_stdout_file_name 2>$picard_stderr_file_name
OUT=$?
echo "Reached line 148"
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

echo "Reached line 172"

## Return to scripts folder at end of script
cd ../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir


