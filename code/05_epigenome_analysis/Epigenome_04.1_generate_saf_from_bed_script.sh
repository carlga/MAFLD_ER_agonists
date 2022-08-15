

cd "/Users/christian.som/OneDrive - KI.SE/OneDrive - Karolinska Institutet/DATA_PhD/Estrogen_Receptor/ChIP_Seq/220124_Fig5_epigenome_new_analysis/220621_REANALYSE/"
bed_path=2.1_Peak_distribution_anno
quant_path=3_Quantification_reads_in_peaks

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${bed_path}/220629_promoters_allDB_182_Mart-anno_unique_shortestDistance.bed > ${quant_path}/220629_allDB_promoters_182_shortestDistance.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${bed_path}/220629_enhancers_allDB_1816_Mart-anno.bed > ${quant_path}/220629_enhancers_allDB_1816_Mart-anno.saf
