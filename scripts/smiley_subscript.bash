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

fasta_to_bed() {
    # Convert a genomic fasta file into a bed file of a certain window size.
    local input=$1
    local size=$2
    #awk '/^>/ {if (seqlen){print substr(seqname, 2) "\t" seqlen}; seqname=$1; seqlen=0; next; } { seqlen += length($0)}END{print substr(seqname, 2) "\t" seqlen}' $1 > genome.txt
    # awk '/^>/ {if (seqlen){print seqname "\t" seqlen}; printf "%s\t", substr($0, 2); seqlen=0; next; } { seqlen += length($0)}END{print seqname "\t" seqlen}' $1 > genome.txt

    # Alternatively, since bedtools doesnt like spaces in the name, replace them all with hyphens
    awk '/^>/ {if (seqlen){print seqname, seqlen}; seqname = substr($0, 2); gsub(" ", "-", seqname); seqlen=0; next; } { seqlen += length($0)}END{print seqname, seqlen}' $1

    bedtools makewindows -g genome.txt -w $2 > genome_windows.bed
}

check_bam_index() {
    # Check if the bam index exists, and if not, create it
    local bam_index=$1
    if [ ! -f $bam_index ]; then
        echo "Bam index not found, creating it now"
        samtools index -b --threads $NUM_CORES $bam
    fi
}

# Split genome into window sizes
fasta_to_bed $GENOME 5000

# Use samtools bedcov to get coverage at each window/region
# -c is the number of reads that are within a region, not just starting inside it
# -x OPTIONAL, if the index file is not in the same folder
check_bam_index $BAM_INDEX
samtools bedcov -c windows.bed $bam > $i.cov



