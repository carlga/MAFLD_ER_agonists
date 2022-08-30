#!/bin/bash
#ChristianSommerauer
#Date:2021-07-16
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 8:00:00
#SBATCH -J 210716_chip-seq_alignment_csommer
#SBATCH --output=210716_mapping.out
#SBATCH --error=210716_mapping.err


#Bowtie version module spider bowtie2/2.3.5.1
module load bioinfo-tools
module load bowtie2

FILE_PATH=
REF_PATH=
OUT_PATH=

###bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | b <bam>} -S [<sam>]

for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

#align two sequences with bowtie2 for Chip-seq
bowtie2 -p 20 -x ${REF_PATH}/mm10_genome -5 3 -U ${FILE_PATH}/${FILE_NAME}.fastq.gz -S ${OUT_PATH}/${FILE_NAME}_mapped.sam

done
