#!/usr/bin/env bash
source ./utils.bash
input_dir="$1"
output_dir="$2"

output_dir_fastp="${output_dir}/fastp_output/"
mkdir -p "$output_dir_fastp"

## fastp
for file in $input_dir/*.f*q*; do

   echo $file
   # Extract the filename without extension
   base_name=$(strip_extension $file)

   # Perform quality checking / preprocessing with fastp
   echo "Running fastp"
   $FASTP_PATH -i "$file" \
                   -o "${output_dir_fastp}/${base_name}_trimmed.fastq" \
                   --json /dev/null \
                   --html /dev/null

   echo "The file $file has been processed."
done

echo "fastp processing is complete."