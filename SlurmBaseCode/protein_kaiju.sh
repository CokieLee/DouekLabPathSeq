#!/bin/sh
#SBATCH -J protein_kaiju
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# fastx-toolkit (fasta_formatter)

## indices ##
# kaiju indices (what reference?)
##########################

codePath=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
nonMammalContigs=$4
origin=$5  ## RNA, DNA or all            
program_prodigal=$6
program_kaiju=$7
kaiju_nodes=$8
kaiju_fmi=$9
outPath=${10}

## print input args
## check that all input args exist and are non-zero
#####################################################
echo "CHECKPOINT 1: PROTEIN_KAIJU INPUTS:"

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

echo "2. left_read_file_base_name:"
variable_is_empty $left_read_file_base_name
echo $left_read_file_base_name
echo "3. right_read_file_base_name:"
variable_is_empty $left_read_file_base_name
echo $right_read_file_base_name

echo "4. nonMammalContigs:"
variable_is_empty $nonMammalContigs
file_exist_check $nonMammalContigs
echo $nonMammalContigs
echo "5. origin:"
echo $origin

echo "6. program_prodigal:"
variable_is_empty $program_prodigal
file_exist_check $program_prodigal
echo $program_prodigal
echo "7. program_kaiju"
variable_is_empty $program_kaiju
file_exist_check $program_kaiju
echo $program_kaiju
## Export paths to prodigal and kaiju programs
# export PATH="$program_prodigal:$PATH"
# export PATH="$program_kaiju:$PATH"

echo "8. kaiju_nodes:"
variable_is_empty $kaiju_nodes
file_exist_check $kaiju_nodes
echo $kaiju_nodes
echo "9. kaiju_fmi:"
variable_is_empty $kaiju_fmi
file_exist_check $kaiju_fmi
echo $kaiju_fmi

echo "10. outPath:"
variable_is_empty $outPath
directory_exists $outPath
echo $outPath

kaiju_output_folder=$outPath/$origin"_kaiju_output/"
mkdir $kaiju_output_folder
directory_exists $kaiju_output_folder
echo "kaiju output folder: $kaiju_output_folder"
#####################################################


echo "CHECKPOINT 2: RUN PRODIGAL:"
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

echo "Prodigal output files:"
non_host_proteins_translations_output_file_faa=$kaiju_output_folder"non_host_proteins_"$origin".faa"
protein_scores_starts_file_txt=$kaiju_output_folder"protein_scores_"$origin".txt"
non_host_proteins_nucleotide_sequences_file_fa=$kaiju_output_folder"non_host_proteins_nucleotide_"$origin".fa"
coords_output_genbank_like_format=$kaiju_output_folder"coords_"$origin".gbk"

echo "non_host_proteins_translations_output_file: $non_host_proteins_translations_output_file_faa"
echo "protein_scores_starts_file_txt: $protein_scores_starts_file_txt"
echo "non_host_proteins_nucleotide_sequences_file_fa: $non_host_proteins_nucleotide_sequences_file_fa"
echo "coords_output_genbank_like_format: $coords_output_genbank_like_format"

prodigal -v

prodigal -a $non_host_proteins_translations_output_file_faa -o $coords_output_genbank_like_format -i $nonMammalContigs -s $protein_scores_starts_file_txt -p meta -d $non_host_proteins_nucleotide_sequences_file_fa
process_fail_check "Prodigal run FAILED. QUITTING"

echo "Prodigal run complete."


echo "CHECKPOINT 3: FORMAT INPUTS:"
formatted_non_host_proteins_translations_output_file_faa=$kaiju_output_folder"formatted_non_host_proteins_"$origin".faa"
formatted_non_host_proteins_nucleotide_sequences_file_fa=$kaiju_output_folder"formatted_non_host_proteins_nucleotide_"$origin".fa"
echo "formatted protein translations output file: $formatted_non_host_proteins_translations_output_file_faa"
echo "formatted protein nucleotide output file: $formatted_non_host_proteins_nucleotide_sequences_file_fa"

## changes the width of sequences line in a FASTA  or FASTAA file
## (all nucleotide/ amino acid sequences appear on a single line)
fasta_formatter -i $non_host_proteins_translations_output_file_faa > $formatted_non_host_proteins_translations_output_file_faa
process_fail_check "Formatting for protein translations output FAILED. QUITTING."

fasta_formatter -i $non_host_proteins_nucleotide_sequences_file_fa > $formatted_non_host_proteins_nucleotide_sequences_file_fa
process_fail_check "Formatting for protein nucleotide output FAILED. QUITTING."

echo "Formatting prodigal outputs complete."


echo "CHECKPOINT 4: RUN KAIJU"

protein_kaiju_output_file_tab=$kaiju_output_folder"protein_kaiju_"$origin".tab"
sorted_protein_kaiju_output_file_tab=$kaiju_output_folder"sorted_protein_kaiju_"$origin".tab"

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
kaiju -t $kaiju_nodes -f $kaiju_fmi -i $formatted_non_host_proteins_translations_output_file_faa -X -e 5 -E 0.01 -p -z 1 -v 1> $protein_kaiju_output_file_tab

## check if process worked
process_fail_check "kaiju run FAILED. QUITTING"

## Sort kaiju output
## The -k option denotes sorting via a key, the 2 denotes sorting on field 2 (confirm)
if [ -s $protein_kaiju_output_file_tab ]
then
  sort -k 2 $protein_kaiju_output_file_tab > $sorted_protein_kaiju_output_file_tab
else
  print "FAILED. $protein_kaiju_output_file_tab was not created by 'kaiju' command"
  exit 1
fi

echo "END OF PROTEIN_KAIJU"