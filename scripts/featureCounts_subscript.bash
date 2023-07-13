#!/usr/bin/env bash

# This is one way to get a quantification of transcripts/reads using a given BAM, GFF file

BAM=$1
GFF=$2
FEATURE=$3 # What feature to use in the GTF/GFF file i.e exon, gene, transcript

featureCounts -F GTF -t $FEATURE -a $GFF -o ${BAM}_counts.txt $BAM
