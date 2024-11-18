#!/bin/sh
#$ -N name_conversion_submission_script
#$ -S /bin/sh
#$ -M rahul.subramanian@nih.gov
#$ -m n
#$ -l h_vmem=8G
#$ -cwd


OLD_IFS=$IFS
echo "IFS test"
echo $IFS
IFS=' '

#Confirm that job is running and print working directory
echo "Job is running" 
full_dir_path=$(pwd)
echo $full_dir_path
## Remove output and log files
rm convert_demux_names_bclfq_to_pathseq.o*
rm convert_demux_names_bclfq_to_pathseq.e*
rm ../name_change_log.txt
## Read in a single read file name
num_BCLFQ_reads=$(ls -l ../Demux_Samples_BCL_to_FQ_format | wc -l)

echo $num_BCLFQ_reads

conversion_script="convert_demux_names_bclfq_to_pathseq.sh"
echo $conversion_script
num_BCLFQ_reads_adj=`expr $num_BCLFQ_reads + 1`
qsub  -t 2-$num_BCLFQ_reads_adj  -tc 100 $conversion_script 
