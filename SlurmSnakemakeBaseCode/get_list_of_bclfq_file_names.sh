#!/bin/sh
#$ -N convert_demux_names_bclfq_to_pathseq
#$ -S /bin/bash
#$ -M rahul.subramanian@nih.gov
#$ -m n
#$ -cwd
#$ -l h_vmem=25G
OLD_IFS=$IFS
IFS="_"

#Confirm that job is running and print working directory
echo "Job is running" 
full_dir_path=$(pwd)
echo $full_dir_path

##Make folder to store file with name list and cd into it
cd ../
mkdir BCLFQ_file_info/
cd BCLFQ_file_info/

##Count the number of files in the folder, and select only the  





cd ../Demux_Samples_BCL_to_FQ_format


IFS=$OLD_IFS