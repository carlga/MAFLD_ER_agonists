#!/bin/bash
#ChristianSommerauer
#Date:2021-07-16
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J 210718_bg_processing
#SBATCH --output=210718_bg_processing.out
#SBATCH --error=210718_bg_processing.err

OUT_PATH=/proj/snic2020-16-225/Chris/ChIPseq/210714_K27_me3_HFD_ESR/210713_CS_HFD_ESR_K4me3_K27ac/Coverage

for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

egrep '^(track|chr)' ${OUT_PATH}/210716_${FILE_NAME}_RPGC_norm.bedgraph  >  ${OUT_PATH}/210716_${FILE_NAME}_RPGC_norm2.bedgraph

echo "track type=bedGraph name=${FILE_NAME} description=${FILE_NAME}  \
visibility=2 color=118,12,96 windowingFunction=maximum alwaysZero=ON" | 
cat - ${OUT_PATH}/210716_${FILE_NAME}_RPGC_norm2.bedgraph > \
${OUT_PATH}/210716_${FILE_NAME}_RPGC_norm_UCSC.bedgraph


done
