#PROMOTERS
#PROMOTERS
#PROMOTERS
#PROMOTERS
#Determine all TSSs in CD and HFD - SPLIT FOR UP AND DOWN and all together

cd "/Users/christian.som/OneDrive - KI.SE/OneDrive - Karolinska Institutet/DATA_PhD/Estrogen_Receptor/ChIP_Seq/"

out_path="220124_Fig5_epigenome_new_analysis/220621_REANALYSE/2_Demarcating_enhancers_promoters"
peak_path="220124_Fig5_epigenome_new_analysis/220621_REANALYSE/PEAK_CALLING_H3K27ac_BROAD"
diffbind_path="220124_Fig5_epigenome_new_analysis/220621_REANALYSE/1_DiffBind_run_220621"


#All Peaks

#CD
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_all_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0736_H3K4me3_CD2_peaks.narrowPeak \
> ${out_path}/prom.1.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD2.bed
#Also do the second replicate
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_all_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0737_H3K4me3_CD9_peaks.narrowPeak \
> ${out_path}/prom.2.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD9.bed

#HFD
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_all_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0739_H3K4me3_HFD4_peaks.narrowPeak \
> ${out_path}/prom.3.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD4.bed
#Also do the second replicate
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_all_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0738_H3K4me3_HFD3_peaks.narrowPeak \
> ${out_path}/prom.4.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD3.bed

#Take the overlap peaks of the respective replicates
bedtools intersect -wa -a ${out_path}/prom.1.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD2.bed \
-b  ${out_path}/prom.2.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD9.bed \
> ${out_path}/prom.5.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.bed

bedtools intersect -wa -a ${out_path}/prom.4.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD3.bed \
-b  ${out_path}/prom.3.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD4.bed \
> ${out_path}/prom.6.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.bed

#Now, take the overlap of both CD and HFD.
module load BEDOPS
sort-bed ${out_path}/prom.5.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.bed > ${out_path}/prom.5.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.sort.bed
sort-bed ${out_path}/prom.6.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.bed > ${out_path}/prom.6.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.sort.bed

bedops --merge ${out_path}/prom.5.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.sort.bed \
${out_path}/prom.6.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.sort.bed \
> ${out_path}/prom.7.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed

#Lastly, merge within the same .bed file to prevent regions from being redundant.
bedtools merge -i ${out_path}/prom.7.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed \
> ${out_path}/prom.8.FINAL.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed

###Results in 184 total Promoter sequences with H3K27ac changes.

#HFDup peaks only. Starts from prom.9

#CD
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDup_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0736_H3K4me3_CD2_peaks.narrowPeak \
> ${out_path}/prom.9.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD2.bed
#Also do the second replicate
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDup_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0737_H3K4me3_CD9_peaks.narrowPeak \
> ${out_path}/prom.10.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD9.bed

#HFD
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDup_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0739_H3K4me3_HFD4_peaks.narrowPeak \
> ${out_path}/prom.11.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD4.bed
#Also do the second replicate
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDup_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0738_H3K4me3_HFD3_peaks.narrowPeak \
> ${out_path}/prom.12.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD3.bed

#Take the overlap peaks of the respective replicates
bedtools intersect -wa -a ${out_path}/prom.9.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD2.bed \
-b  ${out_path}/prom.10.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD9.bed \
> ${out_path}/prom.13.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.bed

bedtools intersect -wa -a ${out_path}/prom.12.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD3.bed \
-b  ${out_path}/prom.11.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD4.bed \
> ${out_path}/prom.14.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.bed

#Now, take the overlap of both CD and HFD.
module load BEDOPS
sort-bed ${out_path}/prom.13.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.bed > ${out_path}/prom.13.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.sort.bed
sort-bed ${out_path}/prom.14.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.bed > ${out_path}/prom.14.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.sort.bed

bedops --merge ${out_path}/prom.13.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.sort.bed \
${out_path}/prom.14.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.sort.bed \
> ${out_path}/prom.15.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed

#Lastly, merge within the same .bed file to prevent regions from being redundant.
bedtools merge -i ${out_path}/prom.15.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed \
> ${out_path}/prom.16.FINAL.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed

###Results in 105 Promoters with H3K27ac changes [ increasing in HFD ].


#HFDdown peaks only. Starts from prom.17

#CD
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDdown_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0736_H3K4me3_CD2_peaks.narrowPeak \
> ${out_path}/prom.17.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD2.bed
#Also do the second replicate
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDdown_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0737_H3K4me3_CD9_peaks.narrowPeak \
> ${out_path}/prom.18.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD9.bed

#HFD
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDdown_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0739_H3K4me3_HFD4_peaks.narrowPeak \
> ${out_path}/prom.19.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD4.bed
#Also do the second replicate
bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDdown_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b  ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0738_H3K4me3_HFD3_peaks.narrowPeak \
> ${out_path}/prom.20.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD3.bed

#Take the overlap peaks of the respective replicates
bedtools intersect -wa -a ${out_path}/prom.17.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD2.bed \
-b  ${out_path}/prom.18.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD9.bed \
> ${out_path}/prom.21.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.bed

bedtools intersect -wa -a ${out_path}/prom.20.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD3.bed \
-b  ${out_path}/prom.19.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD4.bed \
> ${out_path}/prom.22.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.bed

#Now, take the overlap of both CD and HFD.
module load BEDOPS
sort-bed ${out_path}/prom.21.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.bed > ${out_path}/prom.21.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.sort.bed
sort-bed ${out_path}/prom.22.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.bed > ${out_path}/prom.22.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.sort.bed

bedops --merge ${out_path}/prom.21.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_consens.sort.bed \
${out_path}/prom.22.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_HFD_consens.sort.bed \
> ${out_path}/prom.23.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed

#Lastly, merge within the same .bed file to prevent regions from being redundant.
bedtools merge -i ${out_path}/prom.23.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed \
> ${out_path}/prom.24.FINAL.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me3_CD_HFD.bed

###Results in 79 Promoters with H3K27ac changes [ decreasing in HFD ].

#79+105 = 184 total sites. All sites are retained when doing the split analysis.

#Re-create ALL PROMOTERS

#Combined H3K4me3 files
bedtools intersect -wa -a ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0736_H3K4me3_CD2_peaks.narrowPeak  \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0737_H3K4me3_CD9_peaks.narrowPeak \
| sort -k1,1 -k2,2n > ${out_path}/prom.1.genomewide.H3K4me3_CD_consens.bed

bedtools intersect -wa -a ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0739_H3K4me3_HFD4_peaks.narrowPeak  \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0738_H3K4me3_HFD3_peaks.narrowPeak \
| sort -k1,1 -k2,2n > ${out_path}/prom.2.genomewide.H3K4me3_HFD_consens.bed

#Merge CD and HFD
bedops --merge ${out_path}/prom.1.genomewide.H3K4me3_CD_consens.bed \
${out_path}/prom.2.genomewide.H3K4me3_HFD_consens.bed \
> ${out_path}/prom.3.genomewide.FINAL.H3K4me3_CD_HFD_merge.bed

# We have 15,892 total promoters.

### ENHANCERS
### ENHANCERS
### ENHANCERS
### ENHANCERS

#1) Subtract the H3K4me1 data with H3K4me3 data. We do not want H3K4m1 peaks at the TSS.
#remove all H3K4me3 sites from H3K4me1.
bedtools intersect -v -a ${peak_path}/K4me1_broad/H3K4me1_1_ERR2731772_peaks_peaks.broadPeak \
-b ${out_path}/prom.3.genomewide.FINAL.H3K4me3_CD_HFD_merge.bed \
> ${out_path}/enh.1.H3K4me1_1_minus_H3K4me3.bed

bedtools intersect -v -a ${peak_path}/K4me1_broad/H3K4me1_3_ERR2731776_peaks_peaks.broadPeak \
-b ${out_path}/prom.3.genomewide.FINAL.H3K4me3_CD_HFD_merge.bed \
> ${out_path}/enh.2.H3K4me1_3_minus_H3K4me3.bed

#Combine the replicates for H3K4me1 here, which are both reduced by H3K4me3 sites.
bedtools intersect -wa -a ${out_path}/enh.1.H3K4me1_1_minus_H3K4me3.bed \
-b ${out_path}/enh.2.H3K4me1_3_minus_H3K4me3.bed \
> ${out_path}/enh.3.H3K4me1_minus_H3K4me3_consens.bed

# Results in 97540 H3K4me1 regions which are NOT at TSS.
#Now we can take this file to identify H3K27ac sites which have H3K4me1, but without H4me3

#All DiffBound peaks

bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_all_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${out_path}/enh.3.H3K4me1_minus_H3K4me3_consens.bed \
> ${out_path}/enh.4.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.bed

sort -k1,1 -k2,2n ${out_path}/enh.4.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.bed  > ${out_path}/enh.4.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed
bedtools merge -i ${out_path}/enh.4.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed \
> ${out_path}/enh.5.FINAL.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed

#this results in 1816 enhancer sites.

#HFDup diffbound peaks. starts with enh.6.

bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDup_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${out_path}/enh.3.H3K4me1_minus_H3K4me3_consens.bed \
> ${out_path}/enh.6.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.bed

sort -k1,1 -k2,2n ${out_path}/enh.6.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.bed \
> ${out_path}/enh.6.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed

bedtools merge -i ${out_path}/enh.6.HFDup.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed \
> ${out_path}/enh.7.HFDup.FINAL.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed

#this results in 1211 enhancer sites.

#HFDdown diffbound peaks. starts with enh.8.

bedtools intersect -wa -a ${diffbind_path}/220629_diffbind_HFDdown_DB_regions_0585_K27_broad_CD_vs_HFD.bed \
-b ${out_path}/enh.3.H3K4me1_minus_H3K4me3_consens.bed \
> ${out_path}/enh.8.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.bed

sort -k1,1 -k2,2n ${out_path}/enh.8.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.bed \
> ${out_path}/enh.8.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed

bedtools merge -i ${out_path}/enh.8.HFDdown.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed \
> ${out_path}/enh.9.HFDdown.FINAL.H3K27ac_broad_CDHFD_DB_xxx_H3K4me1_minus_me3.sort.bed

#this results in 605 enhancer sites.
# Up and down also work nicely = 1211 + 605 = 1816.

## To have everything in the same script - also do all enhancers.

bedtools intersect -wa -a ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0744_H3K27ac_CD2_peaks.broadPeak \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0745_H3K27ac_CD9_peaks.broadPeak \
> ${out_path}/enh.1.genomewide.H3K27ac_broad_CD_consens.bed

bedtools intersect -wa -a ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0746_H3K27ac_HFD3_peaks.broadPeak \
-b ${peak_path}/Downloaded_Peaks_ArrayExpress/CK0747_H3K27ac_HFD4_peaks.broadPeak \
> ${out_path}/enh.2.genomewide.H3K27ac_broad_HFD_consens.bed

#Now merge them together to a combined H3K27ac file.
bedops --merge ${out_path}/enh.1.genomewide.H3K27ac_broad_CD_consens.bed \
${out_path}/enh.2.genomewide.H3K27ac_broad_HFD_consens.bed \
> ${out_path}/enh.3.genomewide.H3K27ac_all_CD_HFD_merge.bed
# We get 46231 total sites for H3K27ac.

#Consider H3K4me1 (reduced by H3K4me3)
bedtools intersect -wa -a ${out_path}/enh.3.genomewide.H3K27ac_all_CD_HFD_merge.bed \
-b ${out_path}/enh.3.H3K4me1_minus_H3K4me3_consens.bed \
> ${out_path}/enh.4.genomewide.H3K27ac_all_xxx_H3K4me1_minus_K4me3.bed

sort -k1,1 -k2,2n ${out_path}/enh.4.genomewide.H3K27ac_all_xxx_H3K4me1_minus_K4me3.bed  > ${out_path}/enh.4.genomewide.H3K27ac_all_xxx_H3K4me1_minus_K4me3.sort.bed
bedtools merge -i ${out_path}/enh.4.genomewide.H3K27ac_all_xxx_H3K4me1_minus_K4me3.sort.bed  \
> ${out_path}/enh.5.FINAL.genomewide.H3K27ac_all_xxx_H3K4me1_minus_K4me3.bed

wc -l ${out_path}/*bed > ${out_path}/220630_summary_file_rownumbers.txt
awk '{gsub("220124_Fig5_epigenome_new_analysis/220621_REANALYSE/2_Demarcating_enhancers_promoters/", "");print}' \
${out_path}/220630_summary_file_rownumbers.txt > ${out_path}/220630_summary_file_rownumbers.txt

mkdir ${out_path}/Intermediate_enhancer_files
mkdir ${out_path}/Intermediate_promoter_files
mkdir ${out_path}/Final_enhancer_files
mkdir ${out_path}/Final_promoter_files

mv ${out_path}/enh*FINAL* ${out_path}/Final_enhancer_files
mv ${out_path}/prom*FINAL* ${out_path}/Final_promoter_files
mv ${out_path}/prom* ${out_path}/Intermediate_promoter_files
mv ${out_path}/enh* ${out_path}/Intermediate_enhancer_files
