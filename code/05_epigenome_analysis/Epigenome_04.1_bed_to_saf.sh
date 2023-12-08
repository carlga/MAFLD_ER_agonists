#!/usr/bin/env zsh


# Copy this into your terminal (in ~GitHub/MAFLD_ER_agonists/code/05_epigenome_analysis folder)
# or run sh ~GitHub/MAFLD_ER_agonists/code/05_epigenome_analysis/Epigenome_04.1_generate_saf_from_bed_script.sh


bed_path="../../results/Epigenome_analysis"

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${bed_path}/promoters_allDB_142_anno.bed > ${bed_path}/promoters_allDB_142_anno.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${bed_path}/enhancers_allDB_2181_anno.bed > ${bed_path}/enhancers_allDB_2181_anno.saf
