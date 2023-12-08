
module load bioinfo-tools
module load MACS/2.2.6

FILE_PATH=~/PATH_TO_proc_BAMfiles
OUT_PATH=~/OUT_PATH


#THE CONTROL DIETS

macs2 callpeak -t \
${FILE_PATH}/231127_CK0744_H3K27ac_CD2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0744_H3K27ac_CD2 -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_CK0745_H3K27ac_CD9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0745_H3K27ac_CD9 -g 2407883318 -q 0.001 --to-large --broad

#THE HIGH-FAT DIETS

macs2 callpeak -t \
${FILE_PATH}/231127_CK0746_H3K27ac_HFD3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0746_H3K27ac_HFD3 -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_CK0747_H3K27ac_HFD4_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0747_H3K27ac_HFD4 -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_CK0748_H3K27ac_DPN2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0748_H3K27ac_DPN2 -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_CK0749_H3K27ac_DPN3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0749_H3K27ac_DPN3 -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_CK0750_H3K27ac_E2_8_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0750_H3K27ac_E2_8 -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_CK0751_H3K27ac_E2_9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CK0751_H3K27ac_E2_9 -g 2407883318 -q 0.001 --to-large --broad

#2023 samples
macs2 callpeak -t \
${FILE_PATH}/231127_CD6_H3K27ac_S2_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_CD6_Input_S17_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_CD6_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_HFD6_H3K27ac_S8_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_HFD6_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

## DIP
macs2 callpeak -t \
${FILE_PATH}/231127_DIP3_H3K27ac_S9_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_DIP3_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_DIP6_H3K27ac_S10_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_DIP6_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_DIP10_H3K27ac_S11_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_DIP10_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

## PPT
macs2 callpeak -t \
${FILE_PATH}/231127_PPT1_H3K27ac_S12_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_PPT1_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_PPT2_H3K27ac_S13_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_PPT2_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_PPT3_H3K27ac_S14_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_PPT3_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

## Additional E2 and DPN samples
macs2 callpeak -t \
${FILE_PATH}/231127_E2_2_H3K27ac_S16_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_E2_2_H3K27ac -g 2407883318 -q 0.001 --to-large --broad

macs2 callpeak -t \
${FILE_PATH}/231127_DPN6_H3K27ac_S15_R1_001_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/231127_HFD6_Input_S18_R1_001_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/231128_DPN6_H3K27ac -g 2407883318 -q 0.001 --to-large --broad
