#!/usr/bin/env zsh

#computematrices with DeepTools for all promoters and enhancers we identified. For promoters, each TSS has only one peak. For enhancers and promoters, the scaffolds were removed.
cd ../..
bigwig_path="data/files/Bigwig"
prom_enh_path="results/Epigenome_analysis/Bedfiles"
out_path="results/Epigenome_analysis"

  #For the promoters
 computeMatrix reference-point \
 --referencePoint center -S \
 ${bigwig_path}/231127_CK0734_Input_CD9_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0744_H3K27ac_CD2_RPGC_norm.bw \
 ${bigwig_path}/231127_CD6_H3K27ac_S2_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0745_H3K27ac_CD9_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0746_H3K27ac_HFD3_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0747_H3K27ac_HFD4_RPGC_norm.bw \
 ${bigwig_path}/231127_HFD6_H3K27ac_S8_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0748_H3K27ac_DPN2_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0749_H3K27ac_DPN3_RPGC_norm.bw \
 ${bigwig_path}/231127_DPN6_H3K27ac_S15_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_DIP3_H3K27ac_S9_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_DIP6_H3K27ac_S10_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_DIP10_H3K27ac_S11_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_E2_2_H3K27ac_S16_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0750_H3K27ac_E2_8_RPGC_norm.bw \
 ${bigwig_path}/231127_CK0751_H3K27ac_E2_9_RPGC_norm.bw \
 ${bigwig_path}/231127_PPT1_H3K27ac_S12_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_PPT2_H3K27ac_S13_R1_001_RPGC_norm.bw \
 ${bigwig_path}/231127_PPT3_H3K27ac_S14_R1_001_RPGC_norm.bw \
  -R ${prom_enh_path}/promoters_HFDdown_72_anno.bed \
 ${prom_enh_path}/promoters_HFDup_70_anno.bed \
  -b 2000 -a 2000 \
  --outFileName ${out_path}/DB_H3K27ac_up_down_142prom_2K_R.import.gz -p 2

#unzip and remove the header, we will later put it back, but modified
gunzip ${out_path}/DB_H3K27ac_up_down_142prom_2K_R.import.gz

  # For the enhancers
  computeMatrix reference-point \
  --referencePoint center -S \
  ${bigwig_path}/231127_CK0734_Input_CD9_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0744_H3K27ac_CD2_RPGC_norm.bw \
  ${bigwig_path}/231127_CD6_H3K27ac_S2_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0745_H3K27ac_CD9_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0746_H3K27ac_HFD3_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0747_H3K27ac_HFD4_RPGC_norm.bw \
  ${bigwig_path}/231127_HFD6_H3K27ac_S8_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0748_H3K27ac_DPN2_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0749_H3K27ac_DPN3_RPGC_norm.bw \
  ${bigwig_path}/231127_DPN6_H3K27ac_S15_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_DIP3_H3K27ac_S9_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_DIP6_H3K27ac_S10_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_DIP10_H3K27ac_S11_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_E2_2_H3K27ac_S16_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0750_H3K27ac_E2_8_RPGC_norm.bw \
  ${bigwig_path}/231127_CK0751_H3K27ac_E2_9_RPGC_norm.bw \
  ${bigwig_path}/231127_PPT1_H3K27ac_S12_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_PPT2_H3K27ac_S13_R1_001_RPGC_norm.bw \
  ${bigwig_path}/231127_PPT3_H3K27ac_S14_R1_001_RPGC_norm.bw \
   -R  ${prom_enh_path}/enhancers_HFDdown_681_anno.bed \
 ${prom_enh_path}/enhancers_HFDup_1500_anno.bed \
   -b 2000 -a 2000 \
   --outFileName ${out_path}/DB_H3K27ac_broad_up_down_2181enh_2K_R.import.gz -p 2

#unzip and remove the header, we will later put it back, but modified
   gunzip ${out_path}/DB_H3K27ac_broad_up_down_2181enh_2K_R.import.gz
