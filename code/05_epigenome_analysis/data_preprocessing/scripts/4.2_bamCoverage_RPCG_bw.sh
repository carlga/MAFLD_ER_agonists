module load bioinfo-tools
module load deepTools/3.3.2

BAM_path=~/PATH_TO_proc_BAMfiles
OUT_path=~/OUT_PATH


for i in {1..28}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

bamCoverage --bam ${BAM_path}/231127_${FILE_NAME}_psort_BL_fix_MkDup.bam \
--outFileName ${OUT_path}/231127_${FILE_NAME}_RPGC_norm.bw \
--normalizeUsing RPGC \
--effectiveGenomeSize 2407883318 \
--outFileFormat bigwig \
--numberOfProcessors 6

done
