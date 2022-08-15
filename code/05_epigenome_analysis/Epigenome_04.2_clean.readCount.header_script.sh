#Made by Jonas N. Søndergaard
#Made on 190204
#Used to clean up the ugly header from Featurecount output

#grab the header
awk 'NR==2' 220630_diffbind_HFDall_promoters182_H3K27ac.readCount > Ugly.header

#substitute unnecessary text
awk '{gsub("/proj/snic2020-16-225/Chris/ChIPseq/210714_K27_me3_HFD_ESR/210713_CS_HFD_ESR_K4me3_K27ac/Mapping/210716_", "");print}' Ugly.header > Ugly.header.2
awk '{gsub("_psort_BL_fix_MkDup.bam", "");print}' Ugly.header.2 > Nice.header

#grab data without header
awk 'NR!=1 && NR!=2' 220630_diffbind_HFDall_promoters182_H3K27ac.readCount > readCount.copy

#add nice header
cat Nice.header readCount.copy > readCount.copy2

#remove columns that will not be used
cut -f -1,7- readCount.copy2 > 220630_diffbind_HFDall_promoters182_H3K27ac.clean.readCount

#remove files
rm *.copy*
rm *.header*







#grab the header
awk 'NR==2' 220630_diffbind_HFDall_enhancers1816_H3K27ac.readCount > Ugly.header

#substitute unnecessary text
awk '{gsub("/proj/snic2020-16-225/Chris/ChIPseq/210714_K27_me3_HFD_ESR/210713_CS_HFD_ESR_K4me3_K27ac/Mapping/210716_", "");print}' Ugly.header > Ugly.header.2
awk '{gsub("_psort_BL_fix_MkDup.bam", "");print}' Ugly.header.2 > Nice.header

#grab data without header
awk 'NR!=1 && NR!=2' 220630_diffbind_HFDall_enhancers1816_H3K27ac.readCount > readCount.copy

#add nice header
cat Nice.header readCount.copy > readCount.copy2

#remove columns that will not be used
cut -f -1,7- readCount.copy2 > 220630_diffbind_HFDall_enhancers1816_H3K27ac.clean.readCount

#remove files
rm *.copy*
rm *.header*
