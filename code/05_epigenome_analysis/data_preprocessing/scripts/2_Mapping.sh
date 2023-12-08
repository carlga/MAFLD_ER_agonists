
module load bioinfo-tools
module load bowtie2/2.3.5.1

#change paths to according to your directory structure
FILE_PATH= ~/FILE_PATH
REF_PATH= ~/REF_PATH
OUT_PATH= ~/OUT_PATH
###bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | b <bam>} -S [<sam>]
for i in {1..28}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

#align two sequences with bowtie2 for Chip-seq
bowtie2 -p 20 -x ${REF_PATH}/mm10_genome -5 3 -U ${FILE_PATH}/${FILE_NAME}_trim.fastq.gz -S ${OUT_PATH}/${FILE_NAME}_mapped.sam

done
