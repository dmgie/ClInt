#!/usr/bin/env bash

input_dir="$1"
output_dir="$2"
human_reference="$3"

output_dir_trinity="${output_dir}/trinity_output/"
mkdir -p "$output_dir_trinity"

## Trinity DeNovo Assembly
echo "Running Trinity"

for file in "${input_dir}/"*.fastq *.fq; do

  Trinity --seqType fq \
          --max_memory 50G \
          --CPU 6 \
          --single file  

  echo "The file $file has been processed."
done
