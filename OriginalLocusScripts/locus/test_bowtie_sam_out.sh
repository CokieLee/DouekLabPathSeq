module load ncurses/5.9-goolf-1.7.20
module load Bowtie2/2.3.3.1-py27pl5.22
module load FASTX-Toolkit/0.0.14-goolf-1.7.20
#module load BEDTools/2.25.0-goolf-1.7.20
module load SAMtools
#module load Java/1.8.0_92
module load Java/1.8.0_45
module load cURL/7.43.0-goolf-1.7.20


picard=/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/programs/picard/2.18.14/picard.jar

COUNTER=$SGE_TASK_ID
projectID=$1
left_read_file=$2
right_read_file=$3


echo "Line 28"
echo $left_read_file
echo $right_read_file

rel_path_to_proj_dir="../"

left_read_file_rel_path_from_script=${rel_path_to_proj_dir}$left_read_file
right_read_file_rel_path_from_script=${rel_path_to_proj_dir}$right_read_file


echo 'Line 38'
echo $left_read_file_rel_path_from_script

left_read_file_base_name=$(echo ${left_read_file} | cut -d'.' -f1| rev | cut -d'/' -f 1 | rev ) 
right_read_file_base_name=$(echo ${right_read_file} | cut -d'.' -f1| rev | cut -d'/' -f 1 | rev ) 
echo "Line 39"
echo $left_read_file
echo "Line 41"
echo $left_read_file_base_name
echo "Line 43"

## Check that input read files are either in .fastq or .fq.gz format
left_file_format=$(echo ${left_read_file} | cut -d'.' -f2)
if [ $left_file_format != "fastq" ]
then
    if [ $left_file_format != "fq" ]
    then
        printf '%s\n' "left file format invalid " >&1
        printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&1
        echo $left_file_format >&1
        echo $left_read_file_base_name >&1
        echo $left_read_file >&1

        printf '%s\n' "left file format invalid " >&2
        printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&2
        echo $left_file_format >&2
        echo $left_read_file_base_name >&2
        echo $left_read_file >&2

        exit 1
    fi
    gz_left_format=$(echo ${left_read_file} | cut -d'.' -f3)
    if [ $gz_left_format != "gz" ]
    then
       


        printf '%s\n' "left file format invalid: has .fq but no .gz " >&1
        printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&1
        echo $gz_left_format >&1
        echo $left_read_file_base_name >&1

        printf '%s\n' "left file format invalid: has .fq but no .gz " >&2
        printf '%s\n' "PathSeq only accepts .fq.gz and .fastq formats " >&2
        echo $gz_left_format >&2
        echo $left_read_file_base_name >&2

        exit 1

    fi
    echo "left file format is .fq.gz"
    echo $gz_left_format
    
    
fi
#right_read_split_array=split(/\./,$right_read_file)
right_read_file_base_name=$(echo ${right_read_file} | cut -d'.' -f2)
echo $left_read_file_base_name
echo $right_read_file_base_name
echo "line 25"
##########################################
#####

##########################################
inputDir=$4
generated_Data_dir="../Generated_Data/"
outputDir="../genome_alignment_rates/"
erccoutputDir="../ERCC_alignment_rates/"
erccAlignmentCountsDir="../ERCC_alignment_counts/"
insertSizeOutputDir="../insert_size_metrics"


outSam="genome_alignment.sam"

cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="scripts"
echo "correct_cur_Dir"
echo $correct_cur_Dir
if [ $cur_Dir != $correct_cur_Dir ]
then
    printf '%s\n' "In wrong directory should be in " >&1
    echo $correct_cur_Dir >&1
    printf '%s\n' "But instead is in  " >&1
    echo $cur_Dir >&1
    printf '%s\n' "In wrong directory should be in " >&2
    echo $correct_cur_Dir >&2
    printf '%s\n' "But instead is in  " >&2
    echo $cur_Dir >&2
    exit 1
fi


cd ../

cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir=$projectID
echo "correct_cur_Dir"
echo $correct_cur_Dir
if [ $cur_Dir != $correct_cur_Dir ]
then
    printf '%s\n' "In wrong directory should be in " >&1
    echo $correct_cur_Dir >&1
    printf '%s\n' "But instead is in  " >&1
    echo $cur_Dir >&1
    printf '%s\n' "In wrong directory should be in " >&2
    echo $correct_cur_Dir >&2
    printf '%s\n' "But instead is in  " >&2
    echo $cur_Dir >&2
    exit 1
fi

mkdir Generated_Data_1st_Bowtie_Alignment_Unmasked_Genome
mkdir genome_alignment_rates
mkdir insert_size_metrics
mkdir ERCC_alignment_rates
mkdir ERCC_alignment_counts

cd Generated_Data_1st_Bowtie_Alignment_Unmasked_Genome/

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
rm "unalignedRead2AgainstGenome_"$right_read_file_base_name".fq"
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
bowtieIndex=/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/ERCC/ERCC92
bowtieAlignRate="../ERCC_alignment_rates/"$left_read_file_base_name"_ERCC_alignment_rate.txt"


cmd="bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieIndex -1 $left_read_file_rel_path_from_script -2 $right_read_file_rel_path_from_script 1>$outSam 2>$bowtieAlignRate"
echo $cmd
bowtie2 --no-mixed --no-discordant -p 4 -x $bowtieIndex -1 $left_read_file_rel_path_from_script -2 $right_read_file_rel_path_from_script 1>$outSam 2>$bowtieAlignRate