#!/bin/sh
#$ -N convert_demux_names_bclfq_to_pathseq
#$ -S /bin/sh
#$ -M rahul.subramanian@nih.gov
#$ -m n
#$ -l h_vmem=8G
#$ -cwd






OLD_IFS=$IFS
echo "IFS test"
echo $IFS
IFS=' '

DemuxIndex=$SGE_TASK_ID

#Confirm that job is running and print working directory
echo "Job is running" 
full_dir_path=$(pwd)
echo $full_dir_path
echo $DemuxIndex
## Read in a single read file name
single_demux_row=$(ls -l ../Input_Data | awk -v DemuxIndex=$DemuxIndex 'NR==DemuxIndex')
echo $single_demux_row

FileDescArray=($single_demux_row)
Manning_read_file_name=${FileDescArray[8]}
echo $Manning_read_file_name
IFS=$OLD_IFS

##Extract the components of the read file name from the Manning naming conventino
## Anything to the left of _R is the sample number
## Anything to the right of _R is the read number
## We will make our own folders for the lane number (Manning_Lane) and flowcell ID "Manning_Flowcell"
## We assume that all samples are DNA 
IFS='_R'

Manning_read_file_name_base_name=$(basename "${Manning_read_file_name}" .fastq.gz)

Manning_format_file_name_array=($Manning_read_file_name_base_name)
origin="DNA"
sample_ID=${Manning_format_file_name_array[0]}
sample_index=${DemuxIndex}
lane="Manning_Lane"
read=${Manning_format_file_name_array[1]}
IFS=$OLD_IFS
phenotype="human"
adaptor="Unknown"
flow_cell_ID="Manning_Flowcell"

## Isolate all samples with this sample ID:
sample_ID_search_arg=$sample_ID"*"
ls -l $sample_ID_search_arg | wc -l 

## Identify the name of the file to be used as the sample sheet for PathSeq
sample_sheet_file_name="../PathSeq_Sample_Sheet_"$origin".csv"
echo "sample sheet file name"
echo $sample_sheet_file_name

## Create a row in the sample sheet corresponding to this sample
sample_sheet_output_line=$project_type","$project_id","$flow_cell_ID","$lane","$sample_ID","$species","$adaptor","$sample_index","$phenotype","$origin
echo $sample_sheet_output_line
echo $sample_sheet_output_line >> sample_sheet_file_name



##Read in PathSeq sample sheet information for read
#COUNTER=$SGE_TASK_ID
COUNTER=1
OLD_IFS=$IFS
IFS=','
#path_seq_matching_sample=$(awk 'BEGIN { FS = "," } ;{print NR,$5}' ../PathSeq_SampleSheet.csv) 
echo "print test"
echo $BCLFQ_format_file_name_array
echo $sample_ID
echo $lane

## The BCLFFQ lane format is of the form L001. We extract the last three characters. We assume that there are less than 100 lanes. 
lane_full_num=${lane:1:3}

echo "Lane full number"
echo $lane_full_num
##Remove the leading 0s from the lane number
lane_num=$(echo $lane_full_num | awk '{sub(/^0*/,"");}1')
#lane_full_num_int_test=$(awk '{sub(/^0*/,"");}1' $lane_full_num)
echo "lane_num"
echo $lane_num

## We check that there are not more than 4 characters in the Lane string (5 including eof). If not, we throw an error message and print it to both std out and st error and exit
## with error code
lane_word_length=$(echo $lane | wc -c)
echo $lane_word_length
if [ $lane_word_length != 5 ]
then
    printf '%s\n' "BCLFQ Lane variable does not have 5 characters-unable to read lane correctly" >&1
    printf '%s\n' "BCLFQ Lane variable does not have 5 characters-unable to read lane correctly" >&2
    exit 1
fi

echo "Sample index"
echo $sample_index
#Remove the leading S from the sample index in the BCLFQ format
sample_index_adj=$(echo $sample_index  | cut -c 2-)
echo "Sample index adjusted"
echo $sample_index_adj
dos2unix "../PathSeq_SampleSheet.csv"

# Find the path seq sample in the sample sheet that matches the name of the BCLFQ demultiplexed sample
path_seq_matching_sample=$(awk -v well_plate=$well_plate  -v lane_num=$lane_num -v sample_index_adj=$sample_index_adj -v origin=$origin -F',' '(($5 == well_plate) && ($8 == sample_index_adj) && ($9 == origin) && ($4 == lane_num)) {print NR}' ../PathSeq_SampleSheet.csv)

# Check for duplicate matches or if no matches were found, if so, throw error message and exit

num_match_samples=$(echo $path_seq_matching_sample | wc -w)
echo $num_match_samples


if [ $num_match_samples -gt 1 ]
then
    printf '%s\n' "Duplicate samples in samplesheet for single BCLFQ sample" >&1
    echo $num_match_samples >&1
    echo $path_seq_matching_sample >&1
    printf '%s\n' "Duplicate samples in samplesheet for single BCLFQ sample" >&2
    echo $num_match_samples >&2
    echo $path_seq_matching_sample >&2
    exit 1
fi

if [ $num_match_samples -lt 1 ]
then
    printf '%s\n' "No corresponding sample for BCLFQ found in Pathseq Sample Sheet" >&1
    echo $num_match_samples >&1
    printf '%s\n' "No corresponding sample for BCLFQ found in Pathseq Sample Sheet" >&2
    echo $num_match_samples >&2
    exit 1
fi


IFS=$OLD_IFS

file="../PathSeq_SampleSheet.csv"
sampleInformation="$(sed "${path_seq_matching_sample}q;d" $file)"
echo $path_seq_matching_sample
echo $sampleInformation

OLD_IFS=$IFS
IFS=","
newArray=($sampleInformation)

pathSeq_projectType=${newArray[0]}
pathSeq_projectID=${newArray[1]}
pathSeq_flowcellID=${newArray[2]}
pathSeq_lane=${newArray[3]}
pathSeq_index=${newArray[4]}
pathSeq_species=${newArray[5]}
pathSeq_adaptor=${newArray[6]}
pathSeq_sampleID=${newArray[7]}
pathSeq_phenotype=${newArray[8]}
pathSeq_chain=${newArray[9]}

echo $pathSeq_flowcellID
echo $pathSeq_sampleID

#demultiplexDir="/hpcdata/vrc/vrc1_data/douek_lab/Runs/"$flowcellID"/demultiplexed/"$lane"/"$index"/"

## Change directories to the Runs folder, where the demux files for pathseq will be stored
cd ../../../../Runs/

cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir="Runs"
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


##Confirm that we are in the Runs directory
echo "Current directory:"
full_dir_path=$(pwd)
echo $full_dir_path

mkdir $pathSeq_flowcellID
cd $pathSeq_flowcellID
echo $pathSeq_flowcellID

mkdir demultiplexed
cd demultiplexed
mkdir $pathSeq_lane
cd $pathSeq_lane
mkdir $pathSeq_index
cd $pathSeq_index

##Confirm that we are in the sample directory
echo "Current directory:"
full_dir_path=$(pwd)
echo $full_dir_path

cur_Dir=$(basename $(pwd))
echo $cur_Dir
correct_cur_Dir=$pathSeq_index
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

#Remove the R in front of the BCLFQ format read number
read_adj=$(echo $read  | cut -c 2-)
echo "Read adjusted"
echo $read_adj


#pathseq_read_file_name="paired_read1_index_"$index".fq.gz"
pathseq_read_file_name="paired_read"${read_adj}"_index_"${pathSeq_index}".fq.gz"

echo  $pathseq_read_file_name

cp ../../../../../projects/${pathSeq_projectType}/${pathSeq_projectID}/Demux_Samples_BCL_to_FQ_format/${BCLFQ_read_file_name} ${pathseq_read_file_name}

echo ../../../../../projects/${pathSeq_projectType}/${pathSeq_projectID}/Demux_Samples_BCL_to_FQ_format/${BCLFQ_read_file_name} ${pathseq_read_file_name}
#	head -n 40 $filename > ~/Downloads/unix_lesson/shell_review/${samp}_40_lines.fq
#date_today=$(date +%F)
date_today=$(date '+%A %W %Y %X')


output_string="demux_index_"$DemuxIndex", read_"${read_adj}", index_"${pathSeq_index}", BCLFQ_string_"${BCLFQ_read_file_name}", pathseq_file_name_"${pathseq_read_file_name}", date_"${date_today}

cd ../../../../../projects/${pathSeq_projectType}/${pathSeq_projectID}/
echo $output_string >> name_change_log.txt