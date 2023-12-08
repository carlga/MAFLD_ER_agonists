#Enter the path to the coverage files
OUT_PATH=~/OUT_PATH
for i in {1..28}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

egrep '^(track|chr)' ${OUT_PATH}/231127_${FILE_NAME}_RPGC_norm.bedgraph  >  ${OUT_PATH}/231128_${FILE_NAME}_RPGC_norm2.bedgraph

echo "track type=bedGraph name=${FILE_NAME} description=${FILE_NAME}  \
visibility=2 color=118,12,96 windowingFunction=maximum alwaysZero=ON" |
cat - ${OUT_PATH}/231128_${FILE_NAME}_RPGC_norm2.bedgraph > \
${OUT_PATH}/231128_${FILE_NAME}_RPGC_norm_UCSC.bedgraph


done
