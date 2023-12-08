
module load bioinfo-tools
module load samtools/1.12
module load NGSUtils/0.5.9


FILE_PATH= ~/FILE_PATH

for i in {1..28}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

samtools view -S -b ${FILE_PATH}/${FILE_NAME}_mapped.sam > ${FILE_PATH}/231127_${FILE_NAME}_mapped.bam
samtools sort -o ${FILE_PATH}/231127_${FILE_NAME}_sort.bam -@ 8 ${FILE_PATH}/231127_${FILE_NAME}_mapped.bam

rm ${FILE_PATH}/${FILE_NAME}_mapped.sam
rm ${FILE_PATH}/231127_${FILE_NAME}_mapped.bam


#bamutils script to exclude banlist
bamutils filter ${FILE_PATH}/231127_${FILE_NAME}_sort.bam ${FILE_PATH}/231127_${FILE_NAME}_BL.bam -excludebed mm10-blacklist.v2.bed  nostrand


#sort and mark duplicates
samtools sort -o ${FILE_PATH}/231127_${FILE_NAME}_nsort_BL.bam -n -@ 10 ${FILE_PATH}/231127_${FILE_NAME}_BL.bam

#fixmate scores the different read mate coordinates, insert size fields. requires name sorted bam files.
samtools fixmate -m ${FILE_PATH}/231127_${FILE_NAME}_nsort_BL.bam ${FILE_PATH}/231127_${FILE_NAME}_nsort_BL_fix.bam

#markdup only takes position-sorted files, thus do position sorting.
samtools sort -o ${FILE_PATH}/231127_${FILE_NAME}_psort_BL_fix.bam -@ 10 ${FILE_PATH}/231127_${FILE_NAME}_nsort_BL_fix.bam

#mark duplicates. -r removes the duplicates reads as well. -f writes statistics to the specified file
samtools markdup -r -f ${FILE_PATH}/231127_${FILE_NAME}_MARKDUP_STATS.txt ${FILE_PATH}/231127_${FILE_NAME}_psort_BL_fix.bam ${FILE_PATH}/231127_${FILE_NAME}_psort_BL_fix_MkDup.bam

samtools flagstat  ${FILE_PATH}/231127_${FILE_NAME}_psort_BL_fix_MkDup.bam > ${FILE_PATH}/231127_${FILE_NAME}_BEFORE_MkDup.flagstat
samtools flagstat  ${FILE_PATH}/231127_${FILE_NAME}_psort_BL_fix_MkDup.bam > ${FILE_PATH}/231127_${FILE_NAME}_AFTER_MkDup.flagstat

#Indexing
samtools index -b ${FILE_PATH}/231127_${FILE_NAME}_psort_BL_fix_MkDup.bam

done
