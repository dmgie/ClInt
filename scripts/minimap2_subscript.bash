#!/usr/bin/env bash

human_reference="$1"
input_dir="$2"
output_dir="$3"

output_dir_minimap2="${output_dir}/minimap2_output"
mkdir -p "$output_dir_minimap2"

# Run minimap2

for file in "${input_dir}/"*.fastq *fq; do

  # Extract the filename without extension
  base_name=$(basename "$file" .fq)

  $MINIMAP_PATH -a "$human_reference"\
                                  "$file" > \
                                  "${output_dir_minimap2}/${base_name}.sam"

    echo "The file $file has been processed."
done


