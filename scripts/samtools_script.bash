#!/usr/bin/env sh

# https://www.biostars.org/p/49306/ -- use the -M option for much faster

# Create a BED file for each chromosome, with each regions being around 50bp
# apart.

BAM_INDEX=$1
BAM=$2
GENOME=$3

# NOTE: If it was in multiple files, use this
# for i in {1..22} X Y
# do
#     echo "Processing chromosome $i"
#     bedtools makewindows -g $genome -w 500 > $i.bed
# done

# Split genome into window sizes
bedtools makewindows -g $GENOME -w 500 > windows.bed

# Use samtools bedcov to get coverage at each window/region
# -c is the number of reads that are within a region, not just starting inside it
# -x OPTIONAL, if the index file is not in the same folder
samtools bedcov -c windows.bed $bam > $i.cov
