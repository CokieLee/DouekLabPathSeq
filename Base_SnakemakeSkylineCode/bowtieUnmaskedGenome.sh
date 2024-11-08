#!/bin/sh
#SBATCH -J bowtieUnmaskedGenome
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

echo "load open jdk"
module load openjdk/11.0.17_8-kzgd7of
echo "load ncurses"
module load ncurses/6.4-bhkjprr
echo "load bowtie2"
module load bowtie2/2.5.1-bfbrd2r
echo "module load fastx-toolkit"
module load fastx-toolkit/0.0.14-vzydz2l
echo "module load samtools"
module load samtools/1.17-pdu7bhd
echo "module load curl"
module load curl/7.88.1-cj5vz32

sampleNameLeft=$1
sampleNameRight=$2

# must be full paths:
left_read_file=$3
right_read_file=$4

outPath=$5
scriptsPath=$6

bowtieERCCIndex=$7
bowtieUnmaskedGenomeIndex=$8
picard=$9

##############################
# print input paths
echo "CHECKPOINT 1: BOWTIE UNMASKED INPUTS"
echo "1: sampleNameLeft: "
echo $sampleNameLeft
echo "sampleNameRight:"
echo $sampleNameRight

echo "2: left_read_file: "
echo $left_read_file
echo "3: right_read_file: "
echo $right_read_file

echo "5: outPath: "
echo $outPath
echo "6: code path: "
echo $scriptsPath

echo "7: bowtie ercc index: "
echo $bowtieERCCIndex
echo "8: bowtie unmasked genome index: "
echo $bowtieUnmaskedGenomeIndex
echo "9: picard: "
echo $picard
##############################

## Source script for directory checking function
dos2unix $scriptsPath/dir_check.sh
source $scriptsPath/dir_check.sh

$startDir=pwd
##########################################

outSam_ERCC="genome_alignment_ERCC_"$sampleNameLeft".sam"
outSam_Unmasked_Genome="genome_alignment_unmasked_genome_"$sampleNameLeft".sam"

sample_folder_name=$outPath"Sample_"$sampleNameLeft
mkdir $sample_folder_name

#Confirm that sample folder exists
file_exist_check $sample_folder_name

generatedDataFirstAlignDir=$sample_folder_name"/Generated_Data_1st_Bowtie_Alignment_ERCC/"
generatedDataSecondAlignDir=$sample_folder_name"/Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome/"
erccAlignmentRates=$sample_folder_name"/ERCC_alignment_rates/"
erccAlignmentCountsDir=$sample_folder_name"/ERCC_alignment_counts/"
insertSizeOutputDir=$sample_folder_name"/insert_size_metrics/"

mkdir $generatedDataFirstAlignDir
mkdir $generatedDataSecondAlignDir
mkdir $erccAlignmentRates
mkdir $erccAlignmentCountsDir
mkdir $insertSizeOutputDir

cd $generatedDataFirstAlignDir
#Confirm directory change was succesful
correct_cur_Dir=$generatedDataFirstAlignDir

echo "CHECKPOINT 2: Currently inside: "$correct_cur_Dir

##########################################
# rm *.xml
# ##rm *.fq
# rm *.sam
# rm *.bam*
# rm no_ERCC1.fq.gz
# rm no_ERCC2.fq.gz
# rm "unalignedRead1AgainstGenome_"${sampleNameLeft}".fq"
# rm "unalignedRead2AgainstGenome_"${sampleNameRight}".fq"
# rm "unalignedRead1AgainstERCC_"$sampleNameLeft".fq"
# rm "unalignedRead2AgainstERCC_"$sampleNameRight".fq"
# rm "readindex_index_"$index".fq.gz"
# rm left*
# rm right*
# rm *.cxb
# ##rm *.txt
# rm *.csv
# ##rm *.tab
# ##rm *tab.gz
# rm *readcount
##########################################

bowtieAlignRate_ERCC=$erccAlignmentRates"ERCC_alignment_rate_"$sampleNameLeft".txt"

bowtie_cmd="bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieERCCIndex -1 $left_read_file -2 $right_read_file 1>$outSam_ERCC 2>$bowtieAlignRate_ERCC"

echo "CHECKPOINT 3: first bowtie test command:"
echo $bowtie_cmd

echo "$bowtie_cmd" | bash

echo "FINISHED running first bowtie"

outERCCbam="ERCCAlignments_"$sampleNameLeft".bam"
#outERCCbam="ERCCAlignments.bam"
sortedERCCbam="sorted.ERCCAlignments_"$sampleNameLeft".bam"
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
    samtools idxstats $sortedERCCbam 1> $erccAlignmentCountsDir"ERCC_alignment_counts_rate_"$sampleNameLeft".txt"

    samtools view -@ 12 -f 4 $outSam_ERCC | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstERCC_"$sampleNameLeft".fq"
    samtools view -@ 12 -f 4 $outSam_ERCC | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead2AgainstERCC_"$sampleNameRight".fq"

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

echo "FINISHED extracting unaligned reads from first bowtie"

#Confirm directory change was succesful
correct_cur_Dir=$generatedDataSecondAlignDir
cd $generatedDataSecondAlignDir

#bowtieIndex="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg19/unmasked/genome"
bowtieAlignRate_UnmaskedGenome="genome_alignment_rate_"$sampleNameLeft".txt"

ERRC_left_unaligned=$generatedDataFirstAlignDir"unalignedRead1AgainstERCC_"$sampleNameLeft".fq"
ERRC_right_unaligned=$generatedDataFirstAlignDir"unalignedRead2AgainstERCC_"$sampleNameRight".fq"

## Confirm paths to ERCC ouptuts work (files exist)
file_exist_check $ERRC_left_unaligned
file_exist_check $ERRC_right_unaligned

bowtie2_cmd="bowtie2 -p 12 --no-mixed --no-discordant -x $bowtieUnmaskedGenomeIndex -1 $ERRC_left_unaligned -2 $ERRC_right_unaligned 1>$outSam_Unmasked_Genome 2>$bowtieAlignRate_UnmaskedGenome"

echo "CHECKPOINT 4: second bowtie test command:"
echo $bowtie2_cmd
echo "Current working directory: "
pwd

echo "$bowtie2_cmd" | bash

echo "FINISHED running second bowtie"
echo "alignments were sent to: "
echo $outSam_Unmasked_Genome

samtools view -@ 12 -f 4 $outSam_Unmasked_Genome | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead1AgainstGenome_"$sampleNameLeft".fq"
samtools view -@ 12 -f 4 $outSam_Unmasked_Genome | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > "unalignedRead2AgainstGenome_"$sampleNameRight".fq"

echo "FINISHED extracting unaligned reads from second bowtie"

module load r
tempGenomeAlignments_filename="tempGenomeAlignments_"$sampleNameLeft".bam"

samtools view -b $outSam_Unmasked_Genome > $tempGenomeAlignments_filename

tempSortedAlignments_filename="tempSortedAlignments_"$sampleNameLeft".bam"
samtools sort -m 20G $tempGenomeAlignments_filename> $tempSortedAlignments_filename

picard_stdout_file_name="picard_stdout_"$sampleNameLeft".txt"
picard_stderr_file_name="picard_stderr_"$sampleNameLeft".txt" 
insert_size_metrics_file_name="insert_size_metrics_"$sampleNameLeft".txt" 
insert_Hist_pdf_file_name="insert_size_histogram_"$sampleNameLeft".pdf"

echo "CHECKPOINT 5: run picard CollectInsertSizeMetrics command:"

sizeMetrics_cmd="java -Xmx6G -jar $picard CollectInsertSizeMetrics I=$tempSortedAlignments_filename O=$insert_size_metrics_file_name H=$insert_Hist_pdf_file_name M=0.5 1>$picard_stdout_file_name 2>$picard_stderr_file_name"

echo $sizeMetrics_cmd
echo "Current working directory: "
pwd

echo "$sizeMetrics_cmd" | bash
OUT=$?

echo "FINISHED running picard CollectInsertSizeMetrics command"

if [ -e $insert_Hist_pdf_file_name ]
then
	echo "$sampleNameLeft Everything successful ("${OUT}") and deleting intermediate files" >> "../finished_bowtieUnmaskedGenome.txt"
	# rm $outSam_Unmasked_Genome
	# rm $tempSortedAlignments_filename
	# rm $tempGenomeAlignments_filename
	# #rm unmapped.genome.bam
	# #rm temp.unaligned.bam
	# rm *ERCCAlignments*
	# ##rm unalignedRead1AgainstERCC.fq
	# ##rm unalignedRead2AgainstERCC.fq
	mv $insert_size_metrics_file_name "../insert_size_metrics/"$insert_size_metrics_file_name
	mv $insert_Hist_pdf_file_name "../insert_size_metrics/"$insert_Hist_pdf_file_name
	#cp $sampleNameLeft"_insert_size_metrics.txt" $insertSizeOutputDir
	#cp $sampleNameLeft"_insert_size_histogram.pdf" $insertSizeOutputDir
else
    echo $sampleNameLeft" Something went wrong ("${OUT}"), keeping first intermediate alignment file" >> "../finished_bowtieUnmaskedGenome.txt"
    # rm $tempSortedAlignments_filename
	# rm $tempGenomeAlignments_filename
   	# #rm unmapped.genome.bam
   	# #rm tempSortedAlignments.bam
fi

## Return to scripts folder at end of script
cd $startDir
echo "REACHED END OF BOWTIEUNMASKED"