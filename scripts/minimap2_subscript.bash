#!/usr/bin/env bash

human_reference="$1"
input_dir="$2"
output_dir="$3"

output_dir_minimap2="${output_dir}/minimap2_output"
mkdir -p "$output_dir_minimap2"

echo "Start assembly using minimap2"

for file in ../output/fastp_output*.fq; do

  # Extract the filename without extension
  base_name=$(basename "$file" _trimmed.fastq)

    
  ./scripts/minimap2/minimap2 -a "$human_reference"\
                                 "$file" >\ 
                                 "${output_dir_minimap2}/${base_name}_alignment.sam"

  echo "The file $file has been processed."
done


