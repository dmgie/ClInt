#!/usr/bin/env bash

MODE="$1"
input_dir="$2"
human_reference="$3"
output_dir="$4"
num_cores="$5"
max_ram="$6"
max_intron="$7"

## Trinity DeNovo / Genome Guided Assembly

case $MODE in

DeNovo)
  echo "Trinity DeNovo Assembly Mode"
  for file in "${input_dir}/"*fastq; do

    Trinity --seqType fq \
      --max_memory $max_ram \
      --CPU $num_cores \
      --single $file \
      --output ${output_dir}/trinity_output

    echo "The file $file has been processed."
  done
  ;;

Ref_Based)
  echo "Trinity Genome-Guided Assembly Mode"
  for file in "${input_dir}/"*fastq; do
   Trinity --genome_guided_bam $input_dir \
           --genome_guided_max_intron $max_intron \
           --max_memory $max_ram  \
           --CPU $num_cores 
  done
  ;;

*)
  echo "Invalid MODE: $MODE use DeNovo or Ref_Based instead."
  exit 1

esac
