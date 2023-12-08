#!/usr/bin/env zsh

## Determine the closest 3 TSSs (1 left, 2 right) of the differentially acetylated enhancers.


file_path="../../results/Epigenome_analysis"
out_path="../../results/Epigenome_analysis"

sort-bed ${file_path}/promoters_genomewide_12958_anno.bed > ${file_path}/promoters_genomewide_12958_anno.sort.bed

## creates the closest TSS right
closest-features --no-overlaps ${file_path}/enhancers_allDB_2181_anno.bed \
 ${file_path}/promoters_genomewide_12958_anno.sort.bed | awk -F'|' '{print $3}' - > ${out_path}/closest_right_DB_enha_removeMe.bed

## creates the 2nd closest TSS right
closest-features --no-overlaps ${out_path}/closest_right_DB_enha_removeMe.bed \
  ${file_path}/promoters_genomewide_12958_anno.sort.bed | awk -F'|' '{print $3}' - >  ${out_path}/2ndclosest_right_DB_enha_removeMe.bed

## creates the closest TSS left
closest-features --no-overlaps ${file_path}/enhancers_allDB_2181_anno.bed \
 ${file_path}/promoters_genomewide_12958_anno.sort.bed | awk -F'|' '{print $2}' - > ${out_path}/closest_left_DB_enha_removeMe.bed


paste -d'|' ${file_path}/enhancers_allDB_2181_anno.bed \
    ${out_path}/closest_left_DB_enha_removeMe.bed \
    ${out_path}/closest_right_DB_enha_removeMe.bed ${out_path}/2ndclosest_right_DB_enha_removeMe.bed \
    > ${out_path}/origin_closest1_left_closest_1_2_right.bed

rm ${out_path}/*removeMe.bed
