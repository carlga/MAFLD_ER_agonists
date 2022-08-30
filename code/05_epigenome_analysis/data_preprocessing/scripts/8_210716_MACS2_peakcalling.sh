#!/bin/bash
#ChristianSommerauer
#Date: 210717
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -t 04:00:00
#SBATCH -J 210717_peaks
#SBATCH --output=210717_MACS2callpeak.out
#SBATCH --error=210717_MACS2callpeak.err

module load bioinfo-tools
module load MACS/2.2.6

FILE_PATH=
OUT_PATH=


#THE CONTROL DIETS
macs2 callpeak -t \
${FILE_PATH}/210716_CK0736_H3K4me3_CD2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0736_H3K4me3_CD2 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0737_H3K4me3_CD9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0737_H3K4me3_CD9 -g 2407883318 -q 0.01

macs2 callpeak -t \
${FILE_PATH}/210716_CK0744_H3K27ac_CD2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0744_H3K27ac_CD2 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0745_H3K27ac_CD9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0734_Input_CD9_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0745_H3K27ac_CD9 -g 2407883318 -q 0.01 --to-large

#THE HIGH-FAT DIETS
macs2 callpeak -t \
${FILE_PATH}/210716_CK0738_H3K4me3_HFD3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0738_H3K4me3_HFD3 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0739_H3K4me3_HFD4_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0739_H3K4me3_HFD4 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0740_H3K4me3_DPN2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0740_H3K4me3_DPN2 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0741_H3K4me3_DPN3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0741_H3K4me3_DPN3 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0742_H3K4me3_E2_8_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0742_H3K4me3_E2_8 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0743_H3K4me3_E2_9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0743_H3K4me3_E2_9 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0746_H3K27ac_HFD3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0746_H3K27ac_HFD3 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0747_H3K27ac_HFD4_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0747_H3K27ac_HFD4 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0748_H3K27ac_DPN2_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0748_H3K27ac_DPN2 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0749_H3K27ac_DPN3_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0749_H3K27ac_DPN3 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0750_H3K27ac_E2_8_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0750_H3K27ac_E2_8 -g 2407883318 -q 0.01 --to-large

macs2 callpeak -t \
${FILE_PATH}/210716_CK0751_H3K27ac_E2_9_psort_BL_fix_MkDup.bam -c \
${FILE_PATH}/210716_CK0735_Input_HFD3_psort_BL_fix_MkDup.bam -n \
${OUT_PATH}/210716_CK0751_H3K27ac_E2_9 -g 2407883318 -q 0.01 --to-large
