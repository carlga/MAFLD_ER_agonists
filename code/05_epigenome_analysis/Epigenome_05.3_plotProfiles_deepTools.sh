#!/usr/bin/env zsh

bigwig_path="bigwig_files"
prom_enh_path=2.1_Peak_distribution_anno
out_path=4_Metaplots_deeptools

gzip merged_promoters_replicates.mat
gzip merged_enhancers_replicates.mat

plotProfile -m ${out_path}/merged_promoters_replicates.mat.gz \
--perGroup --colors black #88CCEE #332288 #DDCC77 #CC6677 #AA4499  #882255 \
 --plotHeight 10 --plotWidth 13 \
 -out ${out_path}/DB_H3K27ac_up_down_142prom_2K_merged_reps.pdf

# For the enhancers

plotProfile -m ${out_path}/merged_enhancers_replicates.mat.gz \
--perGroup --colors black #88CCEE #332288 #DDCC77 #CC6677 #AA4499  #882255 \
 --plotHeight 10 --plotWidth 13 \
 -out ${out_path}/DB_H3K27ac_broad_up_down_2181enh_2K_merged_reps.pdf
