#!/usr/bin/env bash

## This is the main entry of the workflow bashscript file. We can call some sub-scripts here if 
## we want to make it a bit more modular.

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
  echo -e "\e[1m    $1:\t\e[0m $2 "
}
_log(){
  echo "Ran script at: [$(date +'%Y-%m-%dT%H:%M:%S%z')]: $*"
  echo -e "\e[32mUsed Arguments:\e[0m"
  # Bold escape sequence for the input arguments until the colon
  _arg_log "Input Genome" "$INPUT_GENOME_PATH"
  _arg_log "Human Reference" "$HUMAN_REFERENCE_PATH"
  _arg_log "Human Reference GFF" "$HUMAN_REFERENCE_GFF_PATH"
  _arg_log "Output Directory" "$OUTPUT_DIR"
}

## Argument Parsing and initial logging
_setArgs "$@"
_log


## Set up various (environment) variables that we will use throughout the workflow
NUM_CORES=8
NUM_THREADS=8
MAX_RAM_GB=32

## Set up the paths to the various tools
SPADES_PATH=spades.py
BWA_PATH=bwa
SAMTOOLS_PATH=samtools
BCFTOOLS_PATH=bcftools
FREEBAYES_PATH=freebayes
VCFUTILS_PATH=vcfutils.pl
SNPEFF_PATH=snpeff.jar
SNPSIFT_PATH=SnpSift.jar
SNPEFF_CONFIG_PATH=snpEff.config

############### .... a lot more needed

## Set up the (default) paths to the various input files
# READS_PATH=reads.fastq.gz
# HUMAN_REFERENCE_PATH=human_g1k_v38.fasta
# HUMAN_REFERENCE_GFF_PATH=human_g1k_v38.gff



## Save all arguments in this script to a file
echo "$0 $*" > "$OUTPUT_DIR"/command.txt

## Workflow

