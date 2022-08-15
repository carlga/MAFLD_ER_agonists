#220630: computematrices with DeepTools for all promoters and enhancers we identified. For promoters, each TSS has only one peak. For enhancers and promoters, the scaffolds were removed.


bigwig_path="/Users/christian.som/220124_Fig5_epigenome_new_analysis/bigwig_files"
cd "/Users/christian.som/OneDrive - KI.SE/OneDrive - Karolinska Institutet/DATA_PhD/Estrogen_Receptor/ChIP_Seq/220124_Fig5_epigenome_new_analysis/220621_REANALYSE/"
prom_enh_path=2.1_Peak_distribution_anno
out_path=4_Metaplots_deeptools
 #For the promoters
computeMatrix reference-point \
--referencePoint center -S \
${bigwig_path}/210716_CK0734_Input_CD9_RPGC_norm.bw \
${bigwig_path}/210716_CK0744_H3K27ac_CD2_RPGC_norm.bw \
${bigwig_path}/210716_CK0745_H3K27ac_CD9_RPGC_norm.bw \
${bigwig_path}/210716_CK0746_H3K27ac_HFD3_RPGC_norm.bw \
${bigwig_path}/210716_CK0747_H3K27ac_HFD4_RPGC_norm.bw \
${bigwig_path}/210716_CK0748_H3K27ac_DPN2_RPGC_norm.bw \
${bigwig_path}/210716_CK0749_H3K27ac_DPN3_RPGC_norm.bw \
${bigwig_path}/210716_CK0750_H3K27ac_E2_8_RPGC_norm.bw \
${bigwig_path}/210716_CK0751_H3K27ac_E2_9_RPGC_norm.bw \
 -R ${prom_enh_path}/220629_HFDdown_promoters_78_Mart-anno_unique_shortestDistance.bed \
${prom_enh_path}/220629_HFDup_promoters_104_Mart-anno_unique_shortestDistance.bed \
 -b 2000 -a 2000 \
 --outFileName ${out_path}/220630_DB_H3K27ac_broad_up_down_182prom_2K.gz -p 2

plotProfile -m ${out_path}/220630_DB_H3K27ac_broad_up_down_182prom_2K.gz \
--perGroup \
 -out ${out_path}/220630_DB_H3K27ac_broad_up_down_182prom_2K.pdf

 # For the enhancers
 computeMatrix reference-point \
 --referencePoint center -S \
${bigwig_path}/210716_CK0734_Input_CD9_RPGC_norm.bw \
${bigwig_path}/210716_CK0744_H3K27ac_CD2_RPGC_norm.bw \
${bigwig_path}/210716_CK0745_H3K27ac_CD9_RPGC_norm.bw \
${bigwig_path}/210716_CK0746_H3K27ac_HFD3_RPGC_norm.bw \
${bigwig_path}/210716_CK0747_H3K27ac_HFD4_RPGC_norm.bw \
${bigwig_path}/210716_CK0748_H3K27ac_DPN2_RPGC_norm.bw \
${bigwig_path}/210716_CK0749_H3K27ac_DPN3_RPGC_norm.bw \
${bigwig_path}/210716_CK0750_H3K27ac_E2_8_RPGC_norm.bw \
${bigwig_path}/210716_CK0751_H3K27ac_E2_9_RPGC_norm.bw \
  -R  ${prom_enh_path}/220629_enhancers_HFDdown_605_Mart-anno.bed \
${prom_enh_path}/220629_enhancers_HFDup_1211_Mart-anno.bed \
  -b 2000 -a 2000 \
  --outFileName ${out_path}/220630_DB_H3K27ac_broad_up_down_1816enh_2K.gz -p 2

 plotProfile -m ${out_path}/220630_DB_H3K27ac_broad_up_down_1816enh_2K.gz \
 --perGroup \
  -out ${out_path}/220630_DB_H3K27ac_broad_up_down_1816enh_2K.pdf
