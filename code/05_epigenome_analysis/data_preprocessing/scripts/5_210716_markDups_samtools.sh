#!/bin/bash
#ChristianSommerauer
#Date:210716
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J 210716_mkdup
#SBATCH --output=210716_samtools_markdup.out
#SBATCH --error=210716_samtools_markdup.err

module load bioinfo-tools
module load samtools
#samtool version 1.12

FILE_PATH=

for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

#name sort the bam files. -o out.bam -n (namesorting) -@ threads input.bam
samtools sort -o ${FILE_PATH}/210716_${FILE_NAME}_nsort_BL.bam -n -@ 10 ${FILE_PATH}/210716_${FILE_NAME}_BL.bam

#fixmate scores the different read mate coordinates, insert size fields. requires name sorted bam files.
samtools fixmate -m ${FILE_PATH}/210716_${FILE_NAME}_nsort_BL.bam ${FILE_PATH}/210716_${FILE_NAME}_nsort_BL_fix.bam

#markdup only takes position-sorted files, thus do position sorting.
samtools sort -o ${FILE_PATH}/210716_${FILE_NAME}_psort_BL_fix.bam -@ 10 ${FILE_PATH}/210716_${FILE_NAME}_nsort_BL_fix.bam

#mark duplicates. -r removes the duplicates reads as well. -f writes statistics to the specified file
samtools markdup -r -f ${FILE_PATH}/210716_${FILE_NAME}_MARKDUP_STATS.txt ${FILE_PATH}/210716_${FILE_NAME}_psort_BL_fix.bam ${FILE_PATH}/210716_${FILE_NAME}_psort_BL_fix_MkDup.bam

samtools flagstat  ${FILE_PATH}/210716_${FILE_NAME}_psort_BL_fix_MkDup.bam > ${FILE_PATH}/210716_${FILE_NAME}_BEFORE_MkDup.flagstat
samtools flagstat  ${FILE_PATH}/210716_${FILE_NAME}_psort_BL_fix_MkDup.bam > ${FILE_PATH}/210716_${FILE_NAME}_AFTER_MkDup.flagstat


done
