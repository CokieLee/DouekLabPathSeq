#!/bin/sh
#SBATCH -J array_filter_host
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

module load blast
module load java

codePath=$1
trinityContigFile=$2
outPath=$3
origin=$4  ## RNA, DNA or all
Prog_RemoveHostForKaiju=$5
blastDB_Mammalia=$6

## print input args
## check that all input args exist and are non-zero
#####################################################
echo "CHECKPOINT 1: ALL_TRINITY INPUTS:"

echo "1. codePath:"
echo $codePath
if [[ -d "$codePath" ]]; then
    echo "Directory exists: $codePath"
  else
    echo "Directory does not exist: $codePath. FAILED. QUITTING."
    exit 1
fi
## Source script for directory checking function
dos2unix $codePath"dir_check.sh"
source $codePath"/dir_check.sh"

echo "2. trinityContigFile:"
echo $trinityContigFile
file_exist_check $trinityContigFile

echo "3. outPath:"
echo $outPath
directory_exists $outPath
echo "4. origin:"
echo $origin
variable_is_empty $origin

echo "5. Prog_RemoveHostForKaiju"
variable_is_empty $Prog_RemoveHostForKaiju
echo $Prog_RemoveHostForKaiju
file_exist_check $Prog_RemoveHostForKaiju
echo "6. blastDB_Mammalia"
variable_is_empty $blastDB_Mammalia
echo $blastDB_Mammalia
file_exist_check $blastDB_Mammalia
#####################################################

## Source script for directory checking function
dos2unix $codePath/dir_check.sh
source $codePath/dir_check.sh

## Confirm that the split trinity output file exists
file_exist_check $trinityContigFile

##Make name for .tab output file for blast search
trinity_contig_query_mammal_blast_HSPs=$outPath"/trinity_"$origin"_Sample_"$left_read_file_base_name".tab"

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
blastn -query $trinityContigFile -db $blastDB_Mammalia -outfmt 6 -evalue 0.000000000001 -perc_identity 60 -qcov_hsp_perc 60 -max_target_seqs 1 -num_threads 4 1> $trinity_contig_query_mammal_blast_HSPs

trinity_contigs_not_aligned_to_mammal_BLAST=$outPath"/unaligned_"$trinityContigFile

## This tab file is then fed into the program PathSeqRemoveHostForKaiju.jar along with the split Trininity output file.
## The program filters out matching host contigs, and outputs contigs that could not be aligned to the mammallian reference database. 
java -Xmx5G -jar $Prog_PathSeqRemoveHostForKaiju $trinity_contig_query_mammal_blast_HSPs $trinityContigFile 1> $trinity_contigs_not_aligned_to_mammal_BLAST