
module load bioinfo-tools
module load MACS/2.2.6

FILE_PATH=~/PATH_TO_proc_BAMfiles
OUT_PATH=~/OUT_PATH


#THE CONTROL DIETS
macs2 callpeak -t \
${FILE_PATH}/231127_CK0736_H3K4me3_CD2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231127_CK0736_H3K4me3_CD2 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/231127_CK0737_H3K4me3_CD9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231127_CK0737_H3K4me3_CD9 -g 2407883318 -q 0.01

#THE HIGH-FAT DIETS
macs2 callpeak -t \
${FILE_PATH}/231127_CK0738_H3K4me3_HFD3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231127_CK0738_H3K4me3_HFD3 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/231127_CK0739_H3K4me3_HFD4_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231127_CK0739_H3K4me3_HFD4 -g 2407883318 -q 0.01 --to-large


# Samples from 2023
macs2 callpeak -t \
${FILE_PATH}/231127_CD6_H3K4me3_S1_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CD6_Input_S17_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231127_CD6_H3K4me3 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/231127_HFD6_H3K4me3_S7_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231127_HFD6_H3K4me3 -g 2407883318 -q 0.01 --to-large
