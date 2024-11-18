#!/bin/sh
#SBTACH -J starAfterBowtie
#SBATCH --output=starAfterBowtie.out
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

module load star
module load samtools

COUNTER=$SGE_TASK_ID

projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
hg38_starDB=$4

## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

echo $left_read_file_base_name
echo $right_read_file_base_name
echo "line 25"


#hg38_starDB="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg38/Sequence/STAR/"
#hg_38_gtf="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg38/Annotation/Genes/genes.gtf"
#hg38_referenceFasta="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/hg38/Sequence/WholeGenomeFasta/genome.fa"

echo $hg38_starDB
#echo $hg_38_gtf
#echo $hg38_referenceFasta

demultiplexDir="../Input_Data/"
left="unalignedRead1AgainstGenome_"$left_read_file_base_name".fq"
right="unalignedRead2AgainstGenome_"$right_read_file_base_name".fq"

#outDir_Star="STAR"

##Confirm that we are  in scripts folder
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

##Change directories to project folder
cd "../"

##Confirm that we are in project folder
cur_Dir=$(basename $(pwd))
echo $cur_Dir
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

alignments="../alignment_stats/"
echo "alignments"
echo $alignments

mkdir "alignment_stats"

mkdir Generated_Data_Star_Alignment

cd Generated_Data_Star_Alignment/

## Confirm that we are in the folder Generated_Data_Star_Alignment
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="Generated_Data_Star_Alignment"
dir_check  $correct_cur_Dir


bowtieUnmaskedDir="../Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome/"


pwd

rel_path_to_bowtie_left_out_from_star_folder=$bowtieUnmaskedDir$left
rel_path_to_bowtie_right_out_from_star_folder=$bowtieUnmaskedDir$right

#Confirm that relative paths to unaligned reads work (files exist)
file_exist_check $rel_path_to_bowtie_left_out_from_star_folder
file_exist_check $rel_path_to_bowtie_right_out_from_star_folder

STAR --genomeDir $hg38_starDB --readFilesIn $rel_path_to_bowtie_left_out_from_star_folder $rel_path_to_bowtie_right_out_from_star_folder  --outFileNamePrefix ./$left_read_file_base_name --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt


bam_out_file_name=${left_read_file_base_name}"Aligned.sortedByCoord.out.bam"
log_out_file_name=${left_read_file_base_name}"Log.final.out"
echo "BAM out file name"
echo $bam_out_file_name
echo "Log.out file name"
echo $log_out_file_name
samtools index $bam_out_file_name
samtools idxstats $bam_out_file_name > ${left_read_file_base_name}"_star_alignment_stats.txt"
cp ${left_read_file_base_name}"_star_alignment_stats.txt" $alignments
mv $log_out_file_name ${left_read_file_base_name}"_star_align_summary.txt"
cp ${left_read_file_base_name}"_star_align_summary.txt" $alignments

mv $left_read_file_base_name"Unmapped.out.mate1" $left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq"
mv $left_read_file_base_name"Unmapped.out.mate2" $right_read_file_base_name"_unalignedRead2AgainstTranscriptome.fq"





rm $bowtieUnmaskedDir$left 
rm $bowtieUnmaskedDir$right

lineCount="$(wc -l $left_read_file_base_name"_unalignedRead1AgainstTranscriptome.fq" | cut -d' ' -f1)"
fastqCount=$(expr $lineCount / 4)
echo $left_read_file_base_name","$fastqCount >> "../finished_bowtie_star.csv"

## Change directories back to scripts folder
cd ../../scripts/

##Confirm that we are in scripts directory
echo "Line 140"
cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir
