#!/usr/bin/env zsh

cd ../../

CTCF_path="data"
analysis_path="results/Epigenome_analysis/Enhancer_Expression_correlation"


# this analysis finds diff. acetylated enhancer with CTCF closeby which are minus oriented of gene (left) - enhancer(right) pairs or
# plus oriented for gene (right) - enhancer (left) pairs. This orientation is considered canonical.

bedtools intersect -wa -a data/CTCF_consens_reps.bed \
-b ${analysis_path}/fimo_mm10_genome_CTCF_minus.bed \
> ${analysis_path}/CTCF_minus_peaks.bed

bedtools intersect -wa -a data/CTCF_consens_reps.bed \
-b ${analysis_path}/fimo_mm10_genome_CTCF_plus.bed \
> ${analysis_path}/CTCF_plus_peaks.bed

bedtools intersect -wa -a ${analysis_path}/H3K27ac_left_non_intersect.bed \
-b ${analysis_path}/CTCF_minus_peaks.bed \
> ${analysis_path}/H3K27ac_left_CTCF.intersect.canon.bed

bedtools intersect -wa -a ${analysis_path}/H3K27ac_right_non_intersect.bed \
-b ${analysis_path}/CTCF_plus_peaks.bed \
> ${analysis_path}/H3K27ac_right_CTCF.intersect.canon.bed

## Some entries have duplicate entries because they have multiple CTCF peaks in it. Deduplicate the rows with awk.
awk '!seen[$4]++'  ${analysis_path}/H3K27ac_left_CTCF.intersect.canon.bed > ${analysis_path}/H3K27ac_left_CTCF.intersect.canon.uniq.bed
awk '!seen[$4]++'  ${analysis_path}/H3K27ac_right_CTCF.intersect.canon.bed > ${analysis_path}/H3K27ac_right_CTCF.intersect.canon.uniq.bed



# this analysis finds diff. acetylated enhancer with CTCF closeby which are plus oriented of gene (left) - enhancer(right) pairs or
# minus oriented for gene (right) - enhancer (left) pairs. This orientation is considered non-canonical.

bedtools intersect -wa -a ${analysis_path}/H3K27ac_left_non_intersect.bed \
-b ${analysis_path}/CTCF_plus_peaks.bed \
> ${analysis_path}/H3K27ac_left_CTCF.intersect.noncanon.bed

bedtools intersect -wa -a ${analysis_path}/H3K27ac_right_non_intersect.bed \
-b ${analysis_path}/CTCF_minus_peaks.bed \
> ${analysis_path}/H3K27ac_right_CTCF.intersect.noncanon.bed

## Some entries have duplicate entries because they have multiple CTCF peaks in it. Deduplicate the rows with awk.
# Deduplicating based on the fourth column, so each original enhancer region can get one CTCF peak associated.

awk '!seen[$4]++'  ${analysis_path}/H3K27ac_left_CTCF.intersect.noncanon.bed > ${analysis_path}/H3K27ac_left_CTCF.intersect.noncanon.uniq.bed
awk '!seen[$4]++'  ${analysis_path}/H3K27ac_right_CTCF.intersect.noncanon.bed > ${analysis_path}/H3K27ac_right_CTCF.intersect.noncanon.uniq.bed

# Remove unused files.
rm ${analysis_path}/*canon.bed
rm ${analysis_path}/*sect.bed

echo "script has finished running"
