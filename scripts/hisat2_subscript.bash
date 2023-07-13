#!/usr/bin/env bash


source ./utils.bash
REFERENCE=$1
REF_IDX_NAME="refidx"
READS=$2
NAME_BASE=$(strip_extension ${READS})

# This script creates an index for the reference sequence, and then creates an RNA-seq alignment from the given bam file

# Build reference index
hisat2-build $REFERENCE $REF_IDX_NAME

# Creating first RNA-seq alignment/mapping
hisat2 -x refidx -U $READS -S ${NAME_BASE}.sam

# Sort it by coordinate, useful (and sometimes required) by downstream processes
samtools sort ${NAME_BASE}.sam -o ${NAME_BASE}_coord_sorted.bam
rm ${NAME_BASE}.sam
