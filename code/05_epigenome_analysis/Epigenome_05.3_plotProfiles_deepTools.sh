#!/usr/bin/env zsh

cd ../..
matrix_path="results/Epigenome_analysis"

gzip ${matrix_path}/averaged_coverage_promoters.mat
gzip ${matrix_path}/averaged_coverage_enhancers.mat

plotProfile -m ${matrix_path}/averaged_coverage_promoters.mat.gz \
--perGroup --plotHeight 10 --plotWidth 12 \
-out ${matrix_path}/H3K27ac_avg_coverage_142_promoters.pdf \
--colors black "#88CCEE" "#332288" "#DDCC77" "#CC6677" "#AA4499" "#882255"


# For the enhancers
plotProfile -m ${matrix_path}/averaged_coverage_enhancers.mat.gz \
--perGroup --plotHeight 10 --plotWidth 12 \
-out ${matrix_path}/H3K27ac_avg_coverage_2181_enhancers.pdf \
--colors black "#88CCEE" "#332288" "#DDCC77" "#CC6677" "#AA4499" "#882255"

rm ${matrix_path}/DB_H3K27ac_broad_up_down_2181enh_2K_R.i*
rm ${matrix_path}/DB_H3K27ac_up_down_142prom_2K_R.i*

mkdir ${matrix_path}/Reads_in_peaks_analysis
mv ${matrix_path}/*saf ${matrix_path}/Reads_in_peaks_analysis
mv ${matrix_path}/*mat.gz ${matrix_path}/Reads_in_peaks_analysis
mv ${matrix_path}/H3K27ac*.pdf ${matrix_path}/Reads_in_peaks_analysis
mv ${matrix_path}/*readCount* ${matrix_path}/Reads_in_peaks_analysis
