cd "/Users/christian.som/OneDrive - KI.SE/OneDrive - Karolinska Institutet/DATA_PhD/Estrogen_Receptor/ChIP_Seq/220124_Fig5_epigenome_new_analysis/220621_REANALYSE/"
CTCF_path="5_Correlation_GRE_identification/5.2_CTCF_intersection"



# this analysis finds diff. acetylated enhancer with CTCF closeby which are minus oriented of gene (left) - enhancer(right) pairs or
# plus oriented for gene (right) - enhancer (left) pairs

bedtools intersect -wa -a ${CTCF_path}/PEAK_CALLING/CTCF_consens_reps.bed \
-b ${CTCF_path}/220701_fimo_whole_mm10_genome_CTCF_minus.bed \
> ${CTCF_path}/220701_CTCF_minus_peaks.bed

bedtools intersect -wa -a ${CTCF_path}/PEAK_CALLING/CTCF_consens_reps.bed \
-b ${CTCF_path}/220701_fimo_whole_mm10_genome_CTCF_plus.bed \
> ${CTCF_path}/220701_CTCF_plus_peaks.bed

bedtools intersect -wa -a ${CTCF_path}/220701_8reps_H3K27ac_left_non_intersect.bed \
-b ${CTCF_path}/220701_CTCF_minus_peaks.bed \
> ${CTCF_path}/220701_8rep_H3K27ac_left_CTCF.intersect.canon.bed

bedtools intersect -wa -a ${CTCF_path}/220701_8reps_H3K27ac_right_non_intersect.bed \
-b ${CTCF_path}/220701_CTCF_plus_peaks.bed \
> ${CTCF_path}/220701_8rep_H3K27ac_right_CTCF.intersect.canon.bed

## Some entries have duplicate entries because they have multiple CTCF peaks in it. Deduplicate the rows with awk.
awk '!seen[$4]++'  ${CTCF_path}/220701_8rep_H3K27ac_left_CTCF.intersect.canon.bed > ${CTCF_path}/220701_8rep_H3K27ac_left_CTCF.intersect.canon.uniq.bed
awk '!seen[$4]++'  ${CTCF_path}/220701_8rep_H3K27ac_right_CTCF.intersect.canon.bed > ${CTCF_path}/220701_8rep_H3K27ac_right_CTCF.intersect.canon.uniq.bed



# this analysis finds diff. acetylated enhancer with CTCF closeby which are plus oriented of gene (left) - enhancer(right) pairs or
# minus oriented for gene (right) - enhancer (left) pairs. THIS IS LESS LIKELY TO OCCUR, BUT STILL HAPPENS.


bedtools intersect -wa -a ${CTCF_path}/220701_8reps_H3K27ac_left_non_intersect.bed \
-b ${CTCF_path}/220701_CTCF_plus_peaks.bed \
> ${CTCF_path}/220701_8rep_H3K27ac_left_CTCF.intersect.noncanon.bed

bedtools intersect -wa -a ${CTCF_path}/220701_8reps_H3K27ac_right_non_intersect.bed \
-b ${CTCF_path}/220701_CTCF_minus_peaks.bed \
> ${CTCF_path}/220701_8rep_H3K27ac_right_CTCF.intersect.noncanon.bed

## Some entries have duplicate entries because they have multiple CTCF peaks in it. Deduplicate the rows with awk.
# Deduplicating based on the fourth column, so each original enhancer region can get one CTCF peak associated.

awk '!seen[$4]++'  ${CTCF_path}/220701_8rep_H3K27ac_left_CTCF.intersect.noncanon.bed > ${CTCF_path}/220701_8rep_H3K27ac_left_CTCF.intersect.noncanon.uniq.bed
awk '!seen[$4]++'  ${CTCF_path}/220701_8rep_H3K27ac_right_CTCF.intersect.noncanon.bed > ${CTCF_path}/220701_8rep_H3K27ac_right_CTCF.intersect.noncanon.uniq.bed

# This seems to work! Can re-import the files to further analyze the pairs. The only thing that can be considered to be changed is that the intersect_outwards_orient
# CTCF peaks are not reported yet; but that can be changed in the bedtools intersect command.
