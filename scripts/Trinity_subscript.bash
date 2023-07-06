#!/usr/bin/env bash

MODE="$1"
human_reference="$3"
output_dir="$4"
num_cores="$5"
max_ram="$6"
max_intron="$7"

source $WORKING_DIR/temp_variables.sh

IFS="§"

for program in ${rna_categories["Quality Control"]}; do
  if [ "${rna_programs[$program]}" -eq 1 ]; then
    echo "Aktueller Wert #####: $program"

    ## Trinity DeNovo / Genome Guided Assembly

    input_dir="${output_dir}/${program}_output/"

    echo $input_dir

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

    Guided)
      echo "Trinity Genome-Guided Assembly Mode"
      for file in "${input_dir}/"*fastq; do
        Trinity --genome_guided_bam $input_dir \
          --genome_guided_max_intron $max_intron \
          --max_memory $max_ram \
          --CPU $num_cores
      done
      ;;

    *)
      echo "Invalid MODE: $MODE use DeNovo or Ref_Based instead."
      exit 1
      ;;

    esac

  fi
done
