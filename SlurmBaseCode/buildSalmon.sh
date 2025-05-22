#!/bin/sh
#SBATCH -J buildSalmon
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## REQUIREMENTS
##########################
## programs ##
# java
# salmon

## indices ##
# kaiju indices (what reference?)
##########################

codePath=$1
sorted_protein_kaiju_output_file_tab=$2
formatted_non_host_proteins_nucleotide_sequences_file_fa=$3
formatted_non_host_proteins_translations_output_file_faa=$4
PathSeqKaijuConcensusSplitter2_program=$5
origin=$6  ## RNA, DNA or all
NCBI_nt_kaiju_ref_taxonomy=$7
kaiju_nodes=$8
outPath=$9

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

echo "2. kaiju protein tab output:"
variable_is_empty $sorted_protein_kaiju_output_file_tab
file_exist_check $sorted_protein_kaiju_output_file_tab
echo $sorted_protein_kaiju_output_file_tab
echo "3. kaiju protein nucleotide output:"
variable_is_empty $formatted_non_host_proteins_nucleotide_sequences_file_fa
file_exist_check $formatted_non_host_proteins_nucleotide_sequences_file_fa
echo $formatted_non_host_proteins_nucleotide_sequences_file_fa
echo "4. kaiju protein translations output:"
variable_is_empty $formatted_non_host_proteins_translations_output_file_faa
file_exist_check $formatted_non_host_proteins_translations_output_file_faa
echo $formatted_non_host_proteins_translations_output_file_faa

echo "5. kaiju consensus splitter program:"
variable_is_empty $PathSeqKaijuConcensusSplitter2_program
file_exist_check $PathSeqKaijuConcensusSplitter2_program
echo $PathSeqKaijuConcensusSplitter2_program
echo "6. origin:"
variable_is_empty $origin
echo $origin

echo "7. kaiju taxonomy:"
variable_is_empty $NCBI_nt_kaiju_ref_taxonomy
file_exist_check $NCBI_nt_kaiju_ref_taxonomy
echo $NCBI_nt_kaiju_ref_taxonomy
echo "8. kaiju nodes:"
variable_is_empty $kaiju_nodes
file_exist_check $kaiju_nodes
echo $kaiju_nodes

echo "10. output path:"
variable_is_empty $outPath
directory_exists $outPath
echo $outPath
#####################################################

echo "CHECKPOINT 1: CREATE REFERENCE FASTA FILES FROM KAIJU NONHOST OUTPUT"

startDir=$(pwd)
echo "Starting Dir: $startDir"

salmonIndexOutDir=$outPath"$origin""_salmon_quant/salmon"
mkdir $salmonIndexOutDir
directory_exists $salmonIndexOutDir

echo "Changing directory to: $salmonIndexOutDir"
cd $salmonIndexOutDir

## concensus splitter 
concensusSplitterCmd="java -Xmx90G -jar $PathSeqKaijuConcensusSplitter2_program \
                        $sorted_protein_kaiju_output_file_tab \
                        $formatted_non_host_proteins_nucleotide_sequences_file_fa \
                        $formatted_non_host_proteins_translations_output_file_faa \
                        $NCBI_nt_kaiju_ref_taxonomy \
                        $kaiju_nodes \
                        2>stderr.txt"

echo "Consensus splitter command:"
echo $concensusSplitterCmd
eval "$concensusSplitterCmd"

process_fail_check "Kaiju concensus splitter FAILED. QUITTING"
echo "Finished kaiju concensus splitter."

## TODO: some kind of check so that its okay if some taxonomic levels were not found from kaiju


echo "CHECKPOINT 2: CREATE SALMON INDEX FROM FASTA FILES"

## Get indexes for salmon using the salmon indexer (no pre-defined reference, denoted by the -t flag)
salmon index -t kingdom_sequences.fa -i kingdom_salmon
process_fail_check "kingdom salmon index FAILED"
salmon index -t phylum_sequences.fa -i phylum_salmon
process_fail_check "phylum salmon index FAILED"
salmon index -t class_sequences.fa -i class_salmon
process_fail_check "class salmon index FAILED"
salmon index -t order_sequences.fa -i order_salmon
process_fail_check "order salmon index FAILED"
salmon index -t family_sequences.fa -i family_salmon
process_fail_check "family salmon index FAILED"
salmon index -t genus_sequences.fa -i genus_salmon
process_fail_check "genus salmon index FAILED"
salmon index -t species_sequences.fa -i species_salmon
process_fail_check "species salmon index FAILED"

cd $startDir
echo "returned to $(pwd)"

echo "END OF BUILD SALMON INDEX"