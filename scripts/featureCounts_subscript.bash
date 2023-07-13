#!/usr/bin/env bash

# This is one way to get a quantification of transcripts/reads using a given BAM, GFF file

BAM_DIR=$1
GFF=$2
FEATURE=$3     # What feature to use in the GTF/GFF file i.e exon, gene, transcript
ATTRIBUTE=$4   # In Prokaryotic dataset it is ID
OUTPUT_DIR=$5  

mkdir -p "${OUTPUT_DIR}/featureCounts_output"

for bam_file in "$BAM_DIR"/*.bam; do
  if [ -f "$bam_file" ]; then

    output_file="${bam_file##*/}"
    output_file="${output_file%.bam}_counts.txt"

    featureCounts -F GTF        \
                  -t $FEATURE   \
                  -g $ATTRIBUTE \
                  -T 4          \
                  -a "$GFF"     \
                  -o "$OUTPUT_DIR/featureCounts_output/$output_file" \
                  "$bam_file"
  fi
done
