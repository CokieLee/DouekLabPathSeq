#!/bin/sh
#SBATCH -J prepForDiversity
#SBATCH --mem=25G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cokie.parker@nih.gov

## load Locus modules
module load javafx
module load r


projectID=$1
left_read_file_base_name=$2
right_read_file_base_name=$3
origin=$4 ## RNA, DNA or all
mytaxLevel=$5
PathSeqSplitOutputTableByTaxonomy_program=$6
rScriptDiv=$7

echo "left_read_file_base_name"
echo $left_read_file_base_name

echo "origin"
echo $origin

echo "mytaxLevel"
echo $mytaxLevel


echo "rScriptDiv"
echo $rScriptDiv


## Source script for directory checking function
dos2unix dir_check.sh
source ./dir_check.sh 

## Confirm that we are in the scripts directory
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir

file_exist_check $rScriptDiv

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

##Make diversity metric output folder
mkdir $origin"_diversity"

## Confirm that diversity folder exists
file_exist_check $origin"_diversity"

##Change into that directory
cd $origin"_diversity"

##Confirm we are in that directory
correct_cur_Dir=$origin"_diversity"
dir_check  $correct_cur_Dir

##Make directory for tax level
mkdir $mytaxLevel
file_exist_check $mytaxLevel

cd $mytaxLevel
##Confirm we are in that directory
correct_cur_Dir=$mytaxLevel
dir_check  $correct_cur_Dir

origin_sample_unique_id_tag=$origin"_Sample_"$left_read_file_base_name


cp "../../"$origin"_salmon_quant/"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag"_pseudocounts.csv" .
cp "../../"$origin"_salmon_quant/"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag"_tpm.csv" .

java -Xmx4G -jar $PathSeqSplitOutputTableByTaxonomy_program $origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag"_pseudocounts.csv" k
java -Xmx4G -jar $PathSeqSplitOutputTableByTaxonomy_program $origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag"_tpm.csv" k

full_path_to_rscript_from_div_tax_folder="../../../scripts/"$rScriptDiv

## Confirm we can get to the Rscript (we will be executing it from the current directory
## instead of returning to the scripts folder)
file_exist_check $full_path_to_rscript_from_div_tax_folder

function generateDiversity()
{
	local kingdom=$1
	origin_sample_unique_id_tag_for_func=$2
	file=$kingdom"_"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag_for_func"_pseudocounts.csv"
	diversityOut=$kingdom"_"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag_for_func"_diversity.csv"
	distMatrixOut=$kingdom"_"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag_for_func"_distance.csv"
	distTreeOut=$kingdom"_"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag_for_func"_distance_tree.pdf"
	distHeatMapOut=$kingdom"_"$origin"_"$mytaxLevel"_"$origin_sample_unique_id_tag_for_func"_distance_heatmap.pdf"

	if [ -f "$file" ]; then
		echo "$file exists."
		Rscript $full_path_to_rscript_from_div_tax_folder $file $diversityOut $distMatrixOut $distTreeOut $distHeatMapOut
	else 
		echo "$file does not exist."
	fi
}

kingdom=Eukaryota
generateDiversity $kingdom $origin_sample_unique_id_tag

kingdom=Bacteria
generateDiversity $kingdom $origin_sample_unique_id_tag

kingdom=Viruses
generateDiversity $kingdom $origin_sample_unique_id_tag

kingdom=Archaea
generateDiversity $kingdom $origin_sample_unique_id_tag

## Return to scripts folder at end of script
cd ../../../scripts/

## Confirm that we are in fact in the scripts folder
cur_Dir=$(basename $(pwd))
#echo $cur_Dir
correct_cur_Dir="scripts"
dir_check  $correct_cur_Dir