#!/bin/sh
#SBATCH -J array_filter_host
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

module load blast-plus
module load openjdk

COUNTER=$SGE_TASK_ID
projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4  ## RNA, DNA or all
Prog_PathSeqRemoveHostForKaiju=$5
blastDB_Mammalia=$6



## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

## Confirm that we are in the scripts directory
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

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

## Confirm that the split trinity output folder exists
file_exist_check $trinitySplitDirectory

## Change directories into it
cd $trinitySplitDirectory

## Confirm that we are in the split Trinity directory
correct_cur_Dir=$trinitySplitDirectory
dir_check  $correct_cur_Dir

#trinityOutDirectory="myTrinity_Origin_"$origin"_Sample_"$left_read_file_base_name

#formatted_Trinity_fa_out_file_name="formatted_"$trinityOutDirectory".Trinity.fasta"


#trinityFile=$formatted_Trinity_fa_out_file_name
trinity_prefix="trinity_"$origin"_Sample_"$left_read_file_base_name"_split_"

## Recall that the suffixes for the split trinity output files have 4 digits. We are going to select one of the output files
## corresponding to the array job index (ranges from 1 to fileCount, which was the number of files-see start_filter script for more info).
## However, we can't just select file 1, 2, etc-need to add padding 0s, so it becomes 001, 002 etc. The code below adds that paddding.

##tempCounter is set to array job number (trinity output file number)
tempCounter=$COUNTER
echo $tempCounter
## Count the number of digits of temp counter (number of sigits in suffix)
num_digits_of_suffix_index=$(expr $tempCounter : '.*')
echo $num_digits_of_suffix_index
##If the number of digits of the array index is less than 4, need to pad by adding more 0s
while [[ $num_digits_of_suffix_index < 4 ]]
    	do
    		tempCounter="0$tempCounter" ## Add one 0 in front of number
    		num_digits_of_suffix_index=$(expr $tempCounter : '.*') ## Calculate number of digits (should increase by 1)
    		echo $tempCounter 
    		echo $num_digits_of_suffix_index
    	done
echo $tempCounter
echo $num_digits_of_suffix_index

##Now identify appropriate split trinity output file correspnding to that suffix
trinityFile=$trinity_prefix$tempCounter

## Confirm that it exists
file_exist_check $trinityFile

##Make name for .tab output file for blast search
trinity_contig_query_mammal_blast_HSPs=$trinity_prefix$tempCounter".tab"

## Run a command line blast query (for user manual see https://www.ncbi.nlm.nih.gov/books/NBK569856/,
## also see https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/)
## The -query argument is the file to be queried (in this case the split trinity file)
## The -db argument is the database to be used (in this case a mamillian reference database )
## The -outfmt argument specifies a custom output format-the value 6 produces tab separated rows and columns in a text file
## The e-value (0.000000000001) is defined (see above reference) as 
## the expected number of matches one might find by chance in a subject set of that size with the same score or higher. 
## The e-value arugment is the threshold to use (only high scoring pairs with e-values lower than the threshold should be reported)
## The perc_identity arugment is the percent identity cutoff, which is the perent of bases that are identical to the reference
## genome,  in this case 60% (all high scoring pairs reported must have at least 60% of bases identical to the reference genome.)
## The argument qcov_hsp_perc removes alignemtns below a specified query coverage, in this case 60%. The query coverage is the percent of the contig length that
## aligns to the NCBI hit. In our case, the NCBI hit must align to at least 60% of the split Trinity contig.
## The max_target_seqs argument is the number of aligned seqeuences to keep. In our case, we only keep 1 sequence (this makes sense,
## since any contigs that align to mamallian matches will be filtered, so even one alignement is sufficient since we don't care what specific genes 
## it aligned to as long as it aligned to a mamallian reference genome).
## The num_threads argument is the number of threads (CPUs) to use in the blast search. Here we use 4 threads.
## The output from the blast search for that split Trinity output is stored in a .tab format in the file nametab_file_out_name.
blastn -query $trinityFile -db $blastDB_Mammalia -outfmt 6 -evalue 0.000000000001 -perc_identity 60 -qcov_hsp_perc 60 -max_target_seqs 1 -num_threads 4 1> $trinity_contig_query_mammal_blast_HSPs

trinity_contigs_not_aligned_to_mammal_BLAST="unaligned_"$trinityFile

## This tab file is then fed into the program PathSeqRemoveHostForKaiju.jar along with the split Trininity output file.
## The program filters out matching host contigs, and outputs contigs that could not be aligned to the mammallian reference database. 
java -Xmx5G -jar $Prog_PathSeqRemoveHostForKaiju $trinity_contig_query_mammal_blast_HSPs $trinityFile 1> $trinity_contigs_not_aligned_to_mammal_BLAST

## Return to scripts folder at end of script
cd ../../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir