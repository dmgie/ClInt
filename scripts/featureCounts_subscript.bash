#!/usr/bin/env bash

# This is one way to get a quantification of transcripts/reads using a given BAM, GFF file

BAM_DIR=$1
GFF=$2
FEATURE=$3     # What feature to use in the GTF/GFF file i.e exon, gene, transcript
ATTRIBUTE=$4   # In Prokaryotic dataset it is ID
OUTPUT_DIR=$5  

mkdir -p "${OUTPUT_DIR}/featureCounts_output"

bam_files=$(find "$BAM_DIR" -name "*.bam" -type f -maxdepth 1 -printf "%p ")

featureCounts -F GTF        \
              -t $FEATURE   \
              -g $ATTRIBUTE \
              -T 4          \
              -a "$GFF"     \
              -o "$OUTPUT_DIR/featureCounts_output/counts.txt" \
              $bam_files

# Generate DESeq2 metadata file
metadata_file="${OUTPUT_DIR}/metadata.txt"
echo -e "\condition\ttype" > "$metadata_file"

sample_index=1
for bam_file in $bam_files; do
  sample_name=$(basename "$bam_file" .bam)
  echo -e "$sample_name\ttreated\tsingle-read" >> "$metadata_file"
  ((sample_index++))
done
