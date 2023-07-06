#!/usr/bin/env bash

# Resource: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

GENOME=$1 # The human genome assembly
OUTPUT_DIR=$3 # The output directory

# Source arrays from temp_variables.sh
source $WORKING_DIR/temp_variables.sh

# TODO: Split the string then loop over it
assembly_programs=$(${rna_categories["Assembly"]} | tr '§' ' ')
for program in $assembly_programs; do
    if [[ rna_programs["${program}"] -eq 1 ]]; then
        echo $program
    fi
done

# Preprocessing required to make formatting for RNA-aligned output similar to DNA aligned output
# which is needed for HaploType caller later
echo "Adjusting RNA bam output to allow for Haplotype caller"
gatk SplitNCigarReads \
    -R $GENOME \
    -I $BAM \
    -O ${OUTPUT_DIR}/${BAM}-split.bam

# OPTIONAL: If known polymorphic sites (.vcf file) are there, we can skip over them to exclude from recalibration
# If no known sites, use "--unsorted" (use all sites for recalibration) or "--interval" (specific interval of fasta file used for calibration)
echo "Recalibrating bases"

# Actual Haplotype calling
echo "Running Haplotype caller to analyse variants"
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R $GENOME \
    -I $BAM \
    -O $BAM-output.g.vcf.gz \
    -ERC GVCF
# Alternatively
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R $GENOME \
    -I $BAM \
    -O ${OUTPUT_DIR}/${BAM}-output.g.vcf.gz \
    -bamout $BAM-bamout.bam

# Variant Filtering
# As there is not much truth data for training for the CNN/VQSR methods, use hard filtering
echo "Running variant filtration for..."
gatk VariantFiltration \
    -R $GENOME \
    -V ${OUTPUT_DIR}/${BAM}-output.g.vcf.gz \
    -O output.vcf.gz \
    --filter-name "my_filter1" \
    --filter-expression "AB < 0.2" \
    --filter-name "my_filter2" \
    --filter-expression "MQ0 > 50"
