#!/usr/bin/env bash

source ./utils.bash

reference=$1
input_fastq=$2
output_dir=$3
output_name=$(strip_extension "$input_fastq" ".fastq")
index_dir=${output_dir}/hisat2_output/index/

mkdir -p "$index_dir"

echo "############# $output_dir $output_name"

printf "\n### Start mapping reads to reference using HISAT2 ###\n"
hisat2-build $reference $index_dir

hisat2 -x $index_dir \
       -U $input_fastq \
       | samtools view -bS - > "${output_dir}/hisat2_output/${output_name}.bam"