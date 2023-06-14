#!/usr/bin/env bash

## This is the main entry of the workflow bashscript file. We can call some sub-scripts here if 
## we want to make it a bit more modular.

## Import shared functions
source ./utils.bash

## Set up various (environment) variables that we will use throughout the workflow
NUM_CORES=8
NUM_THREADS=8
MAX_RAM_GB=32

## Argument Parsing: Long and short form
_setArgs(){
  while [ "${1:-}" != "" ]; do
    case "$1" in
      "-i" | "--input-genome")
        shift
        INPUT_GENOME_PATH=$1
        ;;
      "-r" | "--reference")
        shift
        HUMAN_REFERENCE_PATH=$1
        ;;
      "-g" | "--reference-gff")
        shift
        HUMAN_REFERENCE_GFF_PATH=$1
        ;;
      "-o" | "--output-dir")
        shift
        OUTPUT_DIR=$1
        # If output is a file, then we should exit with an error
        if [ -f "$OUTPUT_DIR" ]; then
          echo "Output directory is a file. Please specify a directory."
          exit 1
        fi
        # if it doesn't exist, then we should create it
        if [ ! -d "$OUTPUT_DIR" ]; then
          mkdir -p "$OUTPUT_DIR"
        fi
        ;;
    esac
    shift
  done
}

_arg_log(){
    # Script argument logging function
    echo -e "\e[1m    $1:\t\e[0m $2 "
}
_log(){
    # Log the script arguments to see what was used
    echo "Ran script at: [$(date +'%Y-%m-%dT%H:%M:%S%z')]: $*"
    echo -e "\e[32mUsed Arguments:\e[0m"
    # Bold escape sequence for the input arguments until the colon
    _arg_log "Input Genome" "$INPUT_GENOME_PATH"
    _arg_log "Human Reference" "$HUMAN_REFERENCE_PATH"
    _arg_log "Human Reference GFF" "$HUMAN_REFERENCE_GFF_PATH"
    _arg_log "Output Directory" "$OUTPUT_DIR"
}


## Set up the paths/program names to the various tools
# TODO: Make it read from a template file instead of hardcoding it here
# e.g "SPADES_PATH=" and left empty for them to fill in. This is then read in using a function
# TODO: If a command accepts cores/threads/ram, then we should add it to the command
SPADES_PATH=spades.py
FASTQC_PATH=fastqc
FASTQC_COMMAND="fastqc -t $NUM_CORES -o $OUTPUT_DIR"
BWA_PATH=bwa
SAMTOOLS_PATH=samtools
BCFTOOLS_PATH=bcftools
FREEBAYES_PATH=freebayes
VCFUTILS_PATH=vcfutils.pl
SNPEFF_PATH=snpeff.jar
SNPEFF_CONFIG_PATH=snpEff.config
# .... a lot more needed


## Argument Parsing and initial logging
_setArgs "$@"
_log

## Set up the (default) paths to the various input files
RNA_READS_DIR=reads
# READS_PATH=reads.fastq.gz
# HUMAN_REFERENCE_PATH=human_g1k_v38.fasta
# HUMAN_REFERENCE_GFF_PATH=human_g1k_v38.gff

## Save all arguments in this script to a file, if wanting it to make it more modular
# echo "$0 $*" > "$OUTPUT_DIR"/command.txt

######################################################### Workflow ########################################################

# TODO: Do we do all the question/selection beforehand, or as we go along?

## Assosciative array to store the various options for the user to select from
## The order of the options is also stored in an array, so that we can iterate over it later in order (otherwise it is random)
declare -A categories;                           declare -a order;
categories["Quality Control"]="FastQC MultiQC"; order+=("Quality Control")
categories["Assembly"]="SPAdes Trinity STRING";                order+=("Assembly")
categories["Mapping"]="BWA HISAT STAR";                    order+=("Mapping")
categories["Variant Calling"]="FreeBayes BCFTools"; order+=("Variant Calling")
categories["Annotation"]="SnpEff";              order+=("Annotation")


declare -A programs;                            declare -a order_programs;
programs["FastQC"]=0;                           order_programs+=("FastQC")
programs["MultiQC"]=0;                          order_programs+=("MultiQC")
programs["SPAdes"]=0;                           order_programs+=("SPAdes")
programs["Trinity"]=0;                           order_programs+=("Trinity")
programs["STRING"]=0;                           order_programs+=("STRING")
programs["HISAT"]=0;                            order_programs+=("HISAT")
programs["STAR"]=0;                             order_programs+=("STAR")
programs["BWA"]=0;                              order_programs+=("BWA")
programs["FreeBayes"]=0;                        order_programs+=("FreeBayes")
programs["BCFTools"]=0;                         order_programs+=("BCFTools")
programs["SnpEff"]=0;                           order_programs+=("SnpEff")

## Call display_menu for each category values, and change the value of the selected program to 1. Multiple programs can be selected for each category


# for category in "${order[@]}"; do
#   selected=$(display_menu "Which $category should be used?" ${categories[$category]})
#   programs[$selected]=1
# done

## Echo the selected programs, which will be used in the workflow
echo -e "\e[1mSelected Programs: \e[0m"
for program in "${order_programs[@]}"; do
  if [ "${programs[$program]}" -eq 1 ]; then
    echo "  $program"
  fi
done


echo -e "\e[1mRunning Workflow: \e[0m"

