#!/usr/bin/env zsh

#Made by Jonas N. Søndergaard
#Used to clean up the ugly header from Featurecount output
cd ../../results/Epigenome_analysis
#grab the header
awk 'NR==2' DAc_promoters_142_H3K27ac.readCount > Ugly.header

#substitute unnecessary text. Please adapt this for the individual path.
awk '{gsub("\\../old_new_fastq_combined/Mapping/H3K27ac/231127_", "");print}' Ugly.header > Ugly.header.2
awk '{gsub("_psort_BL_fix_MkDup.bam", "");print}' Ugly.header.2 > Nice.header

#grab data without header
awk 'NR!=1 && NR!=2' DAc_promoters_142_H3K27ac.readCount > readCount.copy

#add nice header
cat Nice.header readCount.copy > readCount.copy2

#remove columns that will not be used
cut -f -1,7- readCount.copy2 > DAc_promoters_142_H3K27ac.clean.readCount

#remove files
rm *.copy*
rm *.header*



#grab the header
awk 'NR==2' DAc_enhancers_2181_H3K27ac.readCount > Ugly.header

#substitute unnecessary text
awk '{gsub("\\../old_new_fastq_combined/Mapping/H3K27ac/231127_", "");print}' Ugly.header > Ugly.header.2
awk '{gsub("_psort_BL_fix_MkDup.bam", "");print}' Ugly.header.2 > Nice.header

#grab data without header
awk 'NR!=1 && NR!=2' DAc_enhancers_2181_H3K27ac.readCount > readCount.copy

#add nice header
cat Nice.header readCount.copy > readCount.copy2

#remove columns that will not be used
cut -f -1,7- readCount.copy2 > DAc_enhancers_2181_H3K27ac.clean.readCount

#remove files
rm *.copy*
rm *.header*