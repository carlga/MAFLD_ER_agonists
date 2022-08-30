#!/bin/bash
#ChristianSommerauer
#Date:210717
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 03:00:00
#SBATCH -J 210717_bam_coverage
#SBATCH --output=210717_bamcoverage_bg.out
#SBATCH --error=210717_bamcoverage_bg.err


module load bioinfo-tools
module load deepTools/3.3.2

BAM_path=
OUT_path=


for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

bamCoverage --bam ${BAM_path}/210716_${FILE_NAME}_psort_BL_fix_MkDup.bam \
--outFileName ${OUT_path}/210716_${FILE_NAME}_RPGC_norm.bedgraph \
--normalizeUsing RPGC \
--effectiveGenomeSize 2407883318 \
--outFileFormat bedgraph \
--numberOfProcessors 6

done
