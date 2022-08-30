#!/bin/bash
#ChristianSommerauer
#Date:210311
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -t 10:00:00
#SBATCH -J 210311_massiveGzip
#SBATCH --output=210311_gunzip.out
#SBATCH --error=210311_gunzip.err


cd /proj/snic2020-16-225/ALLNEW/Wade_Hnf4a_Cebpa_H3K27ac_GSE124463/scripts

gzip *.bam
gzip *fastq
