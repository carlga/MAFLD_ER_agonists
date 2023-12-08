#This script is for the Finger print plot, as a QC.
module load bioinfo-tools
module load deepTools/3.3.2


BAM_PATH=~/PATH_TO_proc_BAMfiles


plotFingerprint \
 -b  ${BAM_PATH}/231127_*H3K4me3*_psort_BL_fix_MkDup.bam \
 ${BAM_PATH}/231127_*Input*_psort_BL_fix_MkDup.bam \
--minMappingQuality 20 --skipZeros \
--numberOfSamples 50000 --smartLabels \
-T "Fingerprints of H3K4me3 in mouse liver vs input"  \
--plotFile ${BAM_PATH}/231128_H3K4me3_MkDup.pdf \
--outRawCounts  ${BAM_PATH}/231128_H3K4me3_MkDup.tab

plotFingerprint \
 -b  ${BAM_PATH}/231127_*H3K27*_psort_BL_fix_MkDup.bam \
 ${BAM_PATH}/231127_*Input*_psort_BL_fix_MkDup.bam \
--minMappingQuality 20 --skipZeros \
--numberOfSamples 50000 --smartLabels \
-T "Fingerprints of H3K27ac in mouse liver vs input"  \
--plotFile ${BAM_PATH}/231128_H3K27ac_MkDup.pdf \
--outRawCounts  ${BAM_PATH}/231128_H3K27ac_MkDup.tab

plotFingerprint \
 -b  ${BAM_PATH}/231127_CK*H3K27*_psort_BL_fix_MkDup.bam \
 ${BAM_PATH}/231127_CK*Input*_psort_BL_fix_MkDup.bam \
--minMappingQuality 20 --skipZeros \
--numberOfSamples 50000 --smartLabels \
-T "Fingerprints of H3K27ac 2021"  \
--plotFile ${BAM_PATH}/231128_H3K27ac_MkDup_2021experiment.pdf \
--outRawCounts  ${BAM_PATH}/231128_H3K27ac_MkDup_2021experiment.tab

plotFingerprint \
 -b  ${BAM_PATH}/231127_*H3K27ac_S*psort_BL_fix_MkDup.bam \
 ${BAM_PATH}/231127*Input*S*_psort_BL_fix_MkDup.bam \
--minMappingQuality 20 --skipZeros \
--numberOfSamples 50000 --smartLabels \
-T "Fingerprints of H3K27ac 2023"  \
--plotFile ${BAM_PATH}/231128_H3K27ac_MkDup_2023experiment.pdf \
--outRawCounts  ${BAM_PATH}/231128_H3K27ac_MkDup_2023experiment.tab
