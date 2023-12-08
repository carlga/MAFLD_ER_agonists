
#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load subread/2.0.0

#file paths
BAM_PATH="../old_new_fastq_combined/Mapping/H3K27ac"
OUTPUT_PATH="."
SAF_PATH="."

#run featureCounts
featureCounts \
        -g gene_id \
        -s 2 \
        -T 2 \
        -M \
	-F SAF \
        -O \
        -C \
        -a ${SAF_PATH}/enhancers_allDB_2181_anno.saf  \
        -o ${OUTPUT_PATH}/DAc_enhancers_2181_H3K27ac.readCount \
        ${BAM_PATH}/*H3K27ac*MkDup.bam \
        &> ${OUTPUT_PATH}/DAc_enhancers_2181_H3K27ac.readCount.log


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