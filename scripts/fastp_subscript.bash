#!/usr/bin/env bash

input_dir="$1"
output_dir="$2"

output_dir_fastp="${output_dir}/fastp_output/"
mkdir -p "$output_dir_fastp"

## fastp
echo "Running fastp"

for file in "$input_dir"*.fastq; do

   # Extract the filename without extension
   base_name=$(basename "$file" .fastq)

   # Perform quality checking / preprocessing with fastp
   ../programs/fastp -i "$file" \
                   -o "${output_dir_fastp}/${base_name}_trimmed.fq" \
                   --json /dev/null \
                   --html /dev/null

   echo "The file $file has been processed."
done

echo "fastp processing is complete."

