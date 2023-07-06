#!/usr/bin/env bash

## This is the main entry of the workflow bashscript file. We can call some sub-scripts here if
## we want to make it a bit more modular.

## Import some shared functions
source ./utils.bash

## XML config file
export XML_FILE="./arguments.xml"

## Argument Parsing: Long and short form
## TODO: Do we even want this? Maybe prompt for user input?
## It might get unwieldy if we have a lot of arguments for every program

export INPUT_READS_PATH="INPUT_READS_PATH_EMPTY"
export HUMAN_REFERENCE_PATH="HUMAN_REFERENCE_PATH_EMPTY"
export HUMAN_REFERENCE_GFF_PATH="HUMAN_REFERENCE_GFF_PATH_EMPTY"
export WORKING_DIR=$(pwd)
export OUTPUT_DIR="OUTPUT_DIR_EMPTY" # TODO: Maybe default this to "$WORKING_DIR/output"
export NUM_CORES="NUM_CORES_EMPTY"
export MAX_RAM="MAX_RAM_EMPTY"


_setArgs() {
  while [ "${1:-}" != "" ]; do
    case "$1" in
    "-i" | "--input-reads")
      shift
      INPUT_READS_PATH=$1
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

_log() {
  # Log the script arguments to see what was used
  echo "Ran script at: [$(date +'%Y-%m-%dT%H:%M:%S%z')]: $*"
  echo -e "\e[32mUsed Arguments:\e[0m"
  # Bold escape sequence for the input arguments until the colon
  _log_format "Input Genome" "$INPUT_READS_PATH"
  _log_format "Human Reference" "$HUMAN_REFERENCE_PATH"
  _log_format "Human Reference GFF" "$HUMAN_REFERENCE_GFF_PATH"
  _log_format "Output Directory" "$OUTPUT_DIR"
}
_log_format() {
  # Script argument logging function
  echo -e "\e[1m    $1:\t\e[0m $2 "
}

## Set up the paths/program names to the various tools
# TODO: Make it read from a template file instead of hardcoding it here
# e.g "SPADES_PATH=" and left empty for them to fill in. This is then read in using a function
# TODO: If a command accepts cores/threads/ram, then we should add it to the command
SPADES_PATH=spades.py
FASTQC_PATH=fastqc
BWA_PATH=bwa
SAMTOOLS_PATH=samtools
BCFTOOLS_PATH=bcftools
FREEBAYES_PATH=freebayes
VCFUTILS_PATH=vcfutils.pl
SNPEFF_PATH=snpeff.jar
SNPEFF_CONFIG_PATH=snpEff.config
FASTP=fastp
MINIMAP2=minimap2
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

## Ask for CORE and RAM
export NUM_CORES=$(get_core_count)
export MAX_RAM=$(get_available_ram)

# Exchange environment variable placeholders for input variables
# Create temporary arguments file
envsubst < arguments.xml \
         > temp_arguments.xml

# Load programs, arguments, descriptions etc. from .xml config file

declare -A rna_categories
declare -a rna_categories_order
declare -A rna_programs
declare -a rna_programs_order

get_config "$XML_FILE" rna_categories rna_categories_order rna_programs rna_programs_order 

# Export the variables so that they can be used in the subscripts
declare -p rna_categories rna_categories_order rna_programs rna_programs_order > temp_variables.sh

######################################################### Workflow ########################################################

## Assosciative array to store the various options for the user to select from
## The order of the options is also stored in an array, so that we can iterate over it later in order (otherwise it is random)
## which would not be ideal when presenting the final option selection

echo -e "----------------\e[1mRunning Workflow: \e[0m-------------------"

# Iterate through program list
for program in "${!rna_programs[@]}"; do

  # Execute when program entry is set to 1
  if [ "${rna_programs[$program]}" -eq 1 ]; then
    # Get arguments from xml config file, turn into array
    arguments=$(get_arguments "program")
    arguments_array=($arguments)

    # Only run program when arguments complete returns 0
    # Execute program subscript
    echo "Running $program using the arguments: $arguments"
    if arguments_complete "${arguments_array[@]}"; then 
      ./scripts/${program}_subscript.bash $arguments
    else
      echo "Arguments missing - unable to execute $program"
    fi
  fi
done

## Delete temporary argument file
rm temp_arguments.xml
rm temp_variables.sh


## TODO: Pause after FastQC, since we need to determine how much we want to trim, so we can ask whether to continue
## TODO: Is this needed? Things like trimgalore and fastp do the trimming for you, so you don't need to pause
##echo -e "\e[1m Read analysis complete, check quality and \e[0m"
##read -p "Press enter to continue..."


# TODO: Adapt towards single end or double end read


