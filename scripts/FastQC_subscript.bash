#!/usr/bin/env bash

input_dir="$1"
output_dir="$2"

output_dir_fastqc="${output_dir}/fastqc_output/"
mkdir -p "$output_dir_fastqc"

## FastQC

for file in $input_dir/*.f*q*; do

   # Perform quality checking with FastQC
    fastqc "$file" -o "$output_dir_fastqc"

  echo "The file $file has been processed."
done

echo "FastQC processing is complete."
