#!/bin/bash
#ChristianSommerauer
#Date:210716
#UPPMAX_commands
#SBATCH -A snic2020-15-291
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J 210617_blacklist
#SBATCH --output=210716_blocklist.out
#SBATCH --error=210716_blocklist.err

module load bioinfo-tools
module load NGSUtils/0.5.9

FILE_PATH=

for i in {1..18}; do \
        FILE_NAME=`sed "${i}q;d" namelist`

bamutils filter ${FILE_PATH}/210716_${FILE_NAME}_sort.bam ${FILE_PATH}/210716_${FILE_NAME}_BL.bam -excludebed mm10-blacklist.v2.bed  nostrand

done


# -excludebed file.bed {nostrand}
                              # Remove reads that are in any of the regions
                              # from the given BED file. If 'nostrand' is given,
                              # strand information from the BED file is ignored.
