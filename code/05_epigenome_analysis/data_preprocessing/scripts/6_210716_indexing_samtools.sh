#!/bin/bash
#ChristianSommerauer
#Date:210716
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -t 02:00:00
#SBATCH -J 210717_indices
#SBATCH --output=210717_bam_indices.out
#SBATCH --error=210717_bam_indices.err

module load bioinfo-tools
module load samtools
#samtool version 1.12

FILE_PATH=
for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

samtools index -b ${FILE_PATH}/210716_${FILE_NAME}_psort_BL_fix_MkDup.bam

done
