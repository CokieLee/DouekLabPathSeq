#!/bin/sh
#SBATCH -J protein_kaiju
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## Load FASTX module
module load fastx-toolkit

left_read_file_base_name=$1
right_read_file_base_name=$2
nonHostContigs=$3
origin=$4  ## RNA, DNA or all            
program_prodigal=$5
program_kaiju=$6
kaiju_nodes=$7
kaiju_fmi=$8
codePath=$9
outPath=${10}

echo "CHECKPOINT 1: PROTEIN_KAIJU INPUTS:"
echo "1. left_read_file_base_name:"
echo $left_read_file_base_name
echo "2. right_read_file_base_name:"
echo $right_read_file_base_name

file_exist_check $nonHostContigs

echo "origin"
echo $origin

## Export paths to prodigal and kaiju programs
export PATH="$program_prodigal:$PATH"
export PATH="$program_kaiju:$PATH"


## Source script for directory checking function
dos2unix $codePath"dir_check.sh"
source $codePath"/dir_check.sh"

## Confirm that sample folder exists
file_exist_check $outPath

kaiju_output_folder=$outPath$origin"_kaiju_output"
mkdir $kaiju_output_folder
## Confirm that the trinity output folder exists
file_exist_check $kaiju_output_folder

## Change directories into that folder
cd $kaiju_output_folder

## Confirm that we are in that folder
correct_cur_Dir=$trinity_output_folder
dir_check  $correct_cur_Dir


## Prodigal is a program that predicts protein coding regions in bacterial and acrchael genomes,
## specifically identifying potential translation initiation sites with confidence scores and RBS
## motifs. 
## For summary of what prodigal does, see the wiki: https://github.com/hyattpd/prodigal/wiki/Introduction

## The -i argument specifies the input file, in this case non_host_formatted_Trinity.fa

## The -o arugment specifies the output file, in this case coords.gbk (genbank like format)
## The-a arugment means specify a protein translations output file,
## which in this case is the file non_host_proteins.faa. 
### The protein translation file lists all the proteins from all the sequences in multiple FASTA format.
### The header consists of: FASTAID_OridnalID #Leftmost coordinate of genome #Rightmost coordinate #strand 
### followed by all the information that would be found in a gene coordinates file
### separated by semi-colons (which 
### would list the location of each gene and scoring info such as GC content)
### For more info see 
### https://github.com/hyattpd/Prodigal/wiki/Understanding-the-Prodigal-Output#protein-translations
### and 
### https://github.com/hyattpd/Prodigal/wiki/Understanding-the-Prodigal-Output#gene-coordinates

### The -s arugment specifies the complete starts file, here protein_scores.txt
### The -p argument indicates the mode of operation to use, here "meta". Meta mode is the same as 
### anon mode, in which prodigal applies pre-calculated training files to the input seqeunce and 
### predicts genes based on the best results. The document 
### recomends that this mode be used on "metagenomes, low quality draft genomes, small viruses,
### and small plasmids". Excerpted from (and more info found on): 
### https://github.com/hyattpd/prodigal/wiki/gene-prediction-modes

### The -d argument specifies the nucleotide sequence file, which is non_host_proteins_nucleotide.fa
origin_sample_unique_id_tag=$origin"_Sample_"$left_read_file_base_name
non_host_proteins_translations_output_file_faa="non_host_proteins_"$origin_sample_unique_id_tag".faa"
protein_scores_starts_file_txt="protein_scores_"$origin_sample_unique_id_tag".txt"
non_host_proteins_nucleotide_sequences_file_fa="non_host_proteins_nucleotide_"$origin_sample_unique_id_tag".fa"
coords_output_genbank_like_format="coords_"$origin_sample_unique_id_tag".gbk"
formatted_non_host_proteins_translations_output_file_faa="formatted_non_host_proteins_"$origin_sample_unique_id_tag".faa"
formatted_non_host_proteins_nucleotide_sequences_file_fa="formatted_non_host_proteins_nucleotide_"$origin_sample_unique_id_tag".fa"

prodigal -a $non_host_proteins_translations_output_file_faa -o $coords_output_genbank_like_format -i $non_host_formatted_Trinity_file_name -s protein_scores.txt -p meta -d $non_host_proteins_nucleotide_sequences_file_fa

## changes the width of sequences line in a FASTA  or FASTAA fil    e
## (all nucleotide/ amino acid sequences appear on a single line)

fasta_formatter -i $non_host_proteins_translations_output_file_faa > $formatted_non_host_proteins_translations_output_file_faa
fasta_formatter -i $non_host_proteins_nucleotide_sequences_file_fa > $formatted_non_host_proteins_nucleotide_sequences_file_fa

protein_kaiju_output_file_tab="protein_kaiju_"$origin_sample_unique_id_tag".tab"

sorted_protein_kaiju_output_file_tab="sorted_protein_kaiju_"$origin_sample_unique_id_tag".tab"

#query="formatted_non_host_proteins.faa"
## Runs kaiju, which performs taxonomic classification.
## Kaiju command descriptions are excerpted from the readme (see link below)
## The -i arugment specifies the input file, in this case 
## $formatted_non_host_proteins_translations_output_file_faa
## The -t and -f arguments specify the files for the reference database
## According to the documentation, the -X option disables filtering of query sequences
## containing low complexity regions, which would otherwise be enabled by default. 
## We note that enabling it is usually recomended to avoid 
## "spurious matches due to simple repeat patterns or other sequencing noise".
## Source: https://github.com/bioinformatics-centre/kaiju/blob/master/README.md

## The -e argument specifies "the exponent of the suffix array checkpoint distances" which
## impacts the tradeoff between search speed and suffix array size. The default 
## value is 5, which is used here. 
## The -E argument specifies the E-value to be used when filtering (in addition to minimum 
## length and score). The default value is 0.01, which is also used here.
## 
## Recall (from earlier) that the e-value is the the expected number of matches 
## one might find by chance in a subject set of that size with the same score or higher.
## The -p argument prints the full taxon path, instead of the default setting 
## which prints only the taxon name.
## The -z option specifies the number of parallel threads that can be used. Here we use 4 threads.

## The -v option enables verbose output, which prints additional columns including if read
## is classified/unclassified, read name, NCBI taxon identifier of assigned taxon,
##  length/score of best match,
## taxon ids and ascension numbers of all database sequences with best match, 
## and matching fragment sequences.

## Kaiju output is stored in a tab file ($protein_kaiju_output_file_tab). 
kaiju -t $kaiju_nodes -f $kaiju_fmi -i $formatted_non_host_proteins_translations_output_file_faa -X -e 5 -E 0.01 -p -z 4 -v 1> $protein_kaiju_output_file_tab

## Sort kaiju output
## The -k option denotes sorting via a key, the 2 denotes sorting on field 2 (confirm)
sort -k 2 $protein_kaiju_output_file_tab > $sorted_protein_kaiju_output_file_tab