#!/usr/bin/env bash

input_dir="$2"
MODE="$1"
human_reference="$3"

output_dir_trinity="${output_dir}/trinity_output/"
mkdir -p "$output_dir_trinity"

## Trinity DeNovo Assembly
echo "Running Trinity"

case $MODE in
    
    DeNovo)
    echo "Trinity DeNovo Assembly Mode"
    ;;

    Ref_Based)
    echo "Trinity Genome-Guided Assembly Mode"
    ;;

esac

for file in "../../resources/fq_files/"*.fq; do

echo "################$file"

  Trinity --seqType fq \
          --max_memory 50G \
          --CPU 6 \
          --single $file  

  echo "The file $file has been processed."
done
