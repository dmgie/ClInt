#!/usr/bin/env bash

MODE="$1"
input_dir="$2"
human_reference="$3"
output_dir="$4"

output_dir_trinity="${output_dir}/trinity_output/"
mkdir -p "$output_dir_trinity"

## Trinity DeNovo / Genome Guided Assembly

case $MODE in

DeNovo)
  echo "Trinity DeNovo Assembly Mode"
  for file in "../../resources/fq_files/"*.fq; do

    echo "################$file"

    Trinity --seqType fq \
      --max_memory 50G \
      --CPU 6 \
      --single $file \
      --output $output_dir

    echo "The file $file has been processed."
  done
  ;;

Ref_Based)
  echo "Trinity Genome-Guided Assembly Mode"
  ;;

esac
