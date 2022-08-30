#!/bin/bash
#Made by Jonas N. Søndergaard, adapted by CS
#Made on 210715
#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 10:00:00
#SBATCH -J 210715_Trimmomatic
#SBATCH --output=210715_Trimmomatic.out
#SBATCH --error=210715_Trimmomatic.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load trimmomatic/0.36
module load FastQC
module load MultiQC

#file paths
FQ_PATH=
OUTPUT_PATH=

#loop to run Trimmomatic for 18 files
for i in {1..18}; do \
	FILE_NAME=`sed "${i}q;d" namelist`

	java -jar $TRIMMOMATIC_HOME/trimmomatic.jar \
		SE \
		-threads 6 \
		-phred33 \
		${FQ_PATH}/${FILE_NAME}.fastq.gz \
		${OUTPUT_PATH}/${FILE_NAME}_trim.fastq.gz \
		ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:1 \
		HEADCROP:4 \
		>>${OUTPUT_PATH}/${FILE_NAME}.trimmomatic.stdout.stderr.txt 2>&1

     fastqc ${OUTPUT_PATH}/${FILE_NAME}_trim.fastq.gz
done

multiqc *zip

#README
# LEADING AND TRAILING WERE REMOVED HERE
#PE: reads are paired end
#-phred33: the quality pipeline used
#ILLUMINACLIP: remove Illumina adapters
#CROP: cut away all bases after this base # from the end of the read
#HEADCROP: cut away the first bases corresponding to the #
#LEADING: remove leading low quality or N bases (below quality 3)
#TRAILING: Remove trailing low quality or N bases (below quality 3)
#SLIDINGWINDOW: Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
#MINLEN: minimum length of reads to keep.
#>> send all messages from Trimmomatic (including errors and warnings) into the specified file
