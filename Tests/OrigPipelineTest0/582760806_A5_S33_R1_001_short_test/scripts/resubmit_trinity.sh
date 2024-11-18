#!/bin/sh
#$ -N resubmit_trinity
#$ -S /bin/bash
#$ -M rahul.subramanian@nih.gov
#$ -m be
#$ -pe threaded 12
#$ -l h_vmem=16.5G

export TMPDIR=/hpcdata/scratch/
export _JAVA_OPTIONS="-Djava.io.tmpdir=/hpcdata/scratch"

module load Trinity/2.13.2-goolf-1.7.20

projectID=$1
file=$2
origin=$3  ## RNA, DNA or all
MIN_CONTIG_LENGTH=$4

cd "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/"$origin"_trinity"

trinityOutDirectory=myTrinity
Trinity --FORCE --CPU 12 --output $trinityOutDirectory --min_contig_length $MIN_CONTIG_LENGTH --seqType fq --max_memory 175G  --min_kmer_cov 1 --left fastp_unmapped_left.fq --right fastp_unmapped_right.fq --no_version_check 1>"trinity_retry_out_"$MIN_CONTIG_LENGTH".txt" 2>trinity_err.txt

exitCode=$?

if [ $exitCode -eq 0 ]
	then
		echo "Trinity rerun was successful for $origin ("$exitCode")" >> "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/finished_Trinity.txt"
		module purge
		module load FASTX-Toolkit/0.0.14-goolf-1.7.20
		fasta_formatter -i myTrinity.Trinity.fasta > formatted_Trinity.fa
		
		module purge
		module load uge
		qsub "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/scripts/start_filter_host_kaiju.sh" $projectID $file $origin 20000
	else
		echo "Trinity failed for $origin ("$exitCode") retrying and STOPPING! You will need to resume the pipeline manually" > "/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/finished_Trinity.txt"

fi