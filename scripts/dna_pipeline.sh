#!/bin/bash

input_dir="$1"
output_dir="$2"
human_reference="$3"

# Create a directory for the output
output_dir_fastqc="${output_dir}/fastqc_output/"
output_dir_fastp="${output_dir}/fastp_output/"
output_dir_spades="${output_dir}/spades_output"
output_dir_minimap2="${output_dir}/minimap2_output"

mkdir -p "$output_dir_fastqc"
mkdir -p "$output_dir_fastp"
mkdir -p "$output_dir_spades"
mkdir -p "$output_dir_minimap2"

# Loop over all .fastq files in the input directory
for file in "$input_dir"/*.fastq; do

    # Extract the filename without extension
    base_name=$(basename "$file" .fastq)

    # Perform quality checking with FastQC
    fastqc "$file" -o "$output_dir_fastqc"

    # Perform quality checking and read trimming wth fastp
    ./scripts/fastp -i "$file" -o "${output_dir_fastp}/${base_name}_trimmed.fq"

    # Assembly using SPAdes
    #spades.py -s "${output_dir_fastp}/${base_name}_trimmed.fq" -o "$output_dir_spades"

    # Assembly using minimap2
    echo "Start assembly"
    ./scripts/minimap2/minimap2 -a "$human_reference" "${output_dir_fastp}/${base_name}_trimmed.fq" > "${output_dir_minimap2}/${base_name}_alignment.sam"

    
    echo "The file $file has been processed."
done

echo "Processing is complete."

