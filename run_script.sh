#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --job-name=clint-nextflow
#SBATCH --ntasks=1
#SBATCH --nodes=1


# Full run
nextflow run ./main.nf --reference_file ../TestData/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --input_dir ../RawData/RNA/NGS05/ --gff_file ../TestData/Homo_sapiens.GRCh38.110.gtf \
    --output_dir proper_runthrough --genome_index star_index/ --paired -profile singularity \
    -resume

