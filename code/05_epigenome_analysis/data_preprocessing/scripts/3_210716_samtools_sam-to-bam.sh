#!/bin/bash
#ChristianSommerauer
#Date:210716
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J 210716_sam_to_bam_sort
#SBATCH --output=210716_sam_to_bam_sort.out
#SBATCH --error=210716_sam_to_bam_sort.err

module load bioinfo-tools
module load samtools/1.12

#samtool version samtools/1.12

FILE_PATH=

for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

samtools view -S -b ${FILE_PATH}/${FILE_NAME}_mapped.sam > ${FILE_PATH}/210716_${FILE_NAME}_mapped.bam
samtools sort -o ${FILE_PATH}/210716_${FILE_NAME}_sort.bam -@ 4 ${FILE_PATH}/210716_${FILE_NAME}_mapped.bam

rm ${FILE_PATH}/${FILE_NAME}_mapped.sam
rm ${FILE_PATH}/210716_${FILE_NAME}_mapped.bam

done
