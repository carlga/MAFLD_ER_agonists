## Determine the closest 3 TSSs (1 left, 2 right) of the differentially acetylated enhancers.

cd "/Users/christian.som/OneDrive - KI.SE/OneDrive - Karolinska Institutet/DATA_PhD/Estrogen_Receptor/ChIP_Seq/220124_Fig5_epigenome_new_analysis/220621_REANALYSE/"
file_path="2.1_Peak_distribution_anno"
out_path="5_Correlation_GRE_identification/5.1_bedops_closest"

sort-bed ${file_path}/220629_HFDup_promoters_13220_Mart-anno_unique_shortestDistance.bed > ${file_path}/220629_HFDup_promoters_13220_Mart-anno_unique_shortestDistance.sort.bed

## creates the closest TSS right
closest-features --no-overlaps ${file_path}/220629_enhancers_allDB_1816_Mart-anno.bed \
 ${file_path}/220629_HFDup_promoters_13220_Mart-anno_unique_shortestDistance.sort.bed | awk -F'|' '{print $3}' - > ${out_path}/220629_closest_right_DB_enha.bed

## creates the 2nd closest TSS right
closest-features --no-overlaps ${out_path}/220629_closest_right_DB_enha.bed \
  ${file_path}/220629_HFDup_promoters_13220_Mart-anno_unique_shortestDistance.sort.bed | awk -F'|' '{print $3}' - >  ${out_path}/220629_2ndclosest_right_DB_enha.bed

## creates the closest TSS left
closest-features --no-overlaps ${file_path}/220629_enhancers_allDB_1816_Mart-anno.bed \
 ${file_path}/220629_HFDup_promoters_13220_Mart-anno_unique_shortestDistance.sort.bed | awk -F'|' '{print $2}' - > ${out_path}/220629_closest_left_DB_enha.bed


paste -d'|' ${file_path}/220629_enhancers_allDB_1816_Mart-anno.bed \
    ${out_path}/220629_closest_left_DB_enha.bed \
    ${out_path}/220629_closest_right_DB_enha.bed ${out_path}/220629_2ndclosest_right_DB_enha.bed \
    > ${out_path}/220629_origin_closest1_left_closest_1_2_right.bed
