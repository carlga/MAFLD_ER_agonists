#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load trimmomatic/0.36
module load FastQC
module load MultiQC

#file paths. change paths to according to your directory structure

FQ_PATH= ~/FQ_PATH
OUTPUT_PATH= ~/OUTPUT_PATH

#loop to run Trimmomatic for 28 files
for i in {1..28}; do \
	FILE_NAME=`sed "${i}q;d" namelist`

	java -jar $TRIMMOMATIC_HOME/trimmomatic.jar \
		SE \
		-threads 6 \
		-phred33 \
		${FQ_PATH}/${FILE_NAME}.fastq.gz \
		${OUTPUT_PATH}/${FILE_NAME}_trim.fastq.gz \
		ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:1 \
		HEADCROP:4 \
		CROP:60 \
		MINLEN:30 \
		>>${OUTPUT_PATH}/${FILE_NAME}.trimmomatic.stdout.stderr.txt 2>&1

     fastqc ${OUTPUT_PATH}/${FILE_NAME}_trim.fastq.gz
done

multiqc ${OUTPUT_PATH}/*zip

#README
#SE: reads are single end
#-phred33: the quality pipeline used
#ILLUMINACLIP: remove Illumina adapters
#CROP: cut away all bases after this base # from the end of the read
#HEADCROP: cut away the first bases corresponding to the #
#LEADING: remove leading low quality or N bases (below quality 3)
#TRAILING: Remove trailing low quality or N bases (below quality 3)
#SLIDINGWINDOW: Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
#MINLEN: minimum length of reads to keep.
#>> send all messages from Trimmomatic (including errors and warnings) into the specified file
