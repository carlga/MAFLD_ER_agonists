#!/bin/bash
#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A snic2021-22-960
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J 220630_FeatureCounts
#SBATCH --output=220630_FeatureCounts_enh.out
#SBATCH --error=220630_FeatureCounts_enh.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load subread/2.0.0

#file paths
BAM_PATH=/proj/snic2020-16-225/Chris/ChIPseq/210714_K27_me3_HFD_ESR/210713_CS_HFD_ESR_K4me3_K27ac/Mapping
OUTPUT_PATH=/proj/snic2020-16-225/Chris/ChIPseq/210714_K27_me3_HFD_ESR/210713_CS_HFD_ESR_K4me3_K27ac/220622_reads_in_peaks_reanalyze
SAF_PATH=/proj/snic2020-16-225/Chris/ChIPseq/210714_K27_me3_HFD_ESR/210713_CS_HFD_ESR_K4me3_K27ac/220622_reads_in_peaks_reanalyze

#run featureCounts
featureCounts \
        -g gene_id \
        -s 2 \
        -T 2 \
        -M \
	-F SAF \
        -O \
        -C \
        -a ${SAF_PATH}/220629_enhancers_allDB_1816_Mart-anno.saf  \
        -o ${OUTPUT_PATH}/220630_diffbind_HFDall_enhancers1816_H3K27ac.readCount \
        ${BAM_PATH}/*H3K27ac*MkDup.bam \
        &> ${OUTPUT_PATH}/220630_diffbind_HFDall_enhancers1816_H3K27ac.readCount.log


#Readme
#-t: Specify feature type in GTF annotation. Here I chose exon, because on transcript level we might not get the lncRNAs
#-g: Specify attribute type in GTF annotation. Here we could chose e.g. transcript ID or gene ID.
#I chose gene ID, because I want to do DE analysis on gene level.
#-s: use '-s 2' if reversely stranded (as is the case for the Illumina Truseq library prep protocol, also NEB NextSeq and generally all UTP methods)
#-T: Number of computational cores/threads used for the analysis
#-p: The experiment is paired end
#-M: Multi-mapping reads will also be counted. Each alignment will have 1 count or a fractional count if --fraction is specified
#-O: Allow reads that overlaps multiple features to be counted
#-C: If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be counted.
#-a: Name of the annotation file. Here it's gencode v27 chr1-22
#-o: output file name


