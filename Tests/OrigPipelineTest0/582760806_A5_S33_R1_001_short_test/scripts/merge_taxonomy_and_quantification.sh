#!/bin/sh
#$ -N mergezaxonomy_and_quantification
#$ -S /bin/bash
#$ -M rahul.subramanian@nih.gov
#$ -m n
#$ -l h_vmem=25G

## load Locus modules
module load java
program=/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/programs/PathSeqMergeQIIME2TaxAndSalmon/dist/PathSeqMergeQIIME2TaxAndSalmon.jar

projectID=$1
file=$2
origin=$3 ## RNA, DNA or all

salmonQuantBase="/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/"$origin"_salmon_quant/"

## part1 is repetitive so let's define and use a function 
salmon_quantification() {
	
	taxLevel=$1
	cd $salmonQuantBase$taxLevel"/"
	taxFile="/hpcdata/vrc/vrc1_data/douek_lab/projects/PathSeq/"$projectID"/"$origin"_trinity/salmon/"$taxLevel"_table.csv"
	pwd
	sampleList=""
	fileList=""
	OLD_IFS=$IFS
	IFS=$'\n'
	for line in $(cat $file ); 
	do
		echo $line
		IFS=","
		newArray=($line)
		projectType=${newArray[0]}
		projectID=${newArray[1]}
		flowcellID=${newArray[2]}
		lane=${newArray[3]}
		index=${newArray[4]}
		phenotype=${newArray[7]}
		sample_origin=${newArray[8]}
		finalName=$phenotype"_"$sample_origin
		finalFile=$flowcellID"_"$lane"_"$index"_"$taxLevel"_quant.sf"
		len=`expr length "$sampleList"`
		echo $len
		if (( $len == 0 ));
		then
			sampleList=$finalName
			fileList=$finalFile
		else
			sampleList=$sampleList","$finalName
			fileList=$fileList","$finalFile
		fi
		echo $sampleList
		echo $fileList
		
	done
	IFS=$OLD_IFS
	echo $sampleList | sed "s/ /,/g"
	echo $fileList | sed "s/ /,/g"
	echo $sampleList
	echo $fileList
	java -Xmx20G -jar $program $fileList $sampleList $taxFile $origin"_"$taxLevel
	mv $origin"_"$taxLevel"_pseudocounts.csv" ../
	mv $origin"_"$taxLevel"_tpm.csv" ../
}

salmon_quantification kingdom
salmon_quantification phylum
salmon_quantification class
salmon_quantification order
salmon_quantification family
salmon_quantification genus
salmon_quantification species
