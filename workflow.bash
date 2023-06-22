#!/usr/bin/env bash

## This is the main entry of the workflow bashscript file. We can call some sub-scripts here if 
## we want to make it a bit more modular.

## Import some shared functions
source ./utils.bash

## Set up various (environment) variables that we will use throughout the workflow, these are default?
## NOTE: Can be skipped, we do prompt-based input instead (below)
# NUM_CORES=8
# NUM_THREADS=8
# MAX_RAM_GB=32

## Argument Parsing: Long and short form
## TODO: Do we even want this? Maybe prompt for user input? 
## It might get unwieldy if we have a lot of arguments for every program
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

_log(){
    # Log the script arguments to see what was used
    echo "Ran script at: [$(date +'%Y-%m-%dT%H:%M:%S%z')]: $*"
    echo -e "\e[32mUsed Arguments:\e[0m"
    # Bold escape sequence for the input arguments until the colon
    _log_format "Input Genome" "$INPUT_GENOME_PATH"
    _log_format "Human Reference" "$HUMAN_REFERENCE_PATH"
    _log_format "Human Reference GFF" "$HUMAN_REFERENCE_GFF_PATH"
    _log_format "Output Directory" "$OUTPUT_DIR"
}
_log_format(){
    # Script argument logging function
    echo -e "\e[1m    $1:\t\e[0m $2 "
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
FASTP=./scripts/fastp
MINIMAP2=./scripts/minimap2/minimap2
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

## Assosciative array to store the various options for the user to select from
## The order of the options is also stored in an array, so that we can iterate over it later in order (otherwise it is random)
## which would not be ideal when presenting the final option selection
declare -A "rna_programs";         
declare -a rna_programs_order;
declare -A rna_categories;                                      
declare -a rna_categories_order;

RNA_Questionnaire() {
    # RNA-based questionnaire
    
    rna_categories["Quality Control"]="FastQC,\
MultiQC,\
fastp (All-in-one),\
TrimGalore";                                                     rna_categories_order+=("Quality Control")
    rna_categories["Assembly"]="SPAdes,Trinity,STRING,minimap2"; rna_categories_order+=("Assembly")
    rna_categories["Mapping"]="BWA,HISAT,STAR";                  rna_categories_order+=("Mapping")
    rna_categories["Variant Calling"]="FreeBayes,BCFTools";      rna_categories_order+=("Variant Calling")
    rna_categories["Annotation"]="SnpEff";                       rna_categories_order+=("Annotation")


    rna_programs["FastQC"]=0;                   rna_programs_order+=("FastQC")
    rna_programs["MultiQC"]=0;                  rna_programs_order+=("MultiQC")
    rna_programs["fastp"]=0;                    rna_programs_order+=("fastp")
    rna_programs["TrimGalore"]=0;               rna_programs_order+=("TrimGalore")
    rna_programs["SPAdes"]=0;                   rna_programs_order+=("SPAdes")
    rna_programs["Trinity"]=0;                  rna_programs_order+=("Trinity")
    rna_programs["STRING"]=0;                   rna_programs_order+=("STRING")
    rna_programs["minimap2"]=0;                 rna_programs_order+=("minimap2")
    rna_programs["HISAT"]=0;                    rna_programs_order+=("HISAT")
    rna_programs["STAR"]=0;                     rna_programs_order+=("STAR")
    rna_programs["BWA"]=0;                      rna_programs_order+=("BWA")
    rna_programs["FreeBayes"]=0;                rna_programs_order+=("FreeBayes")
    rna_programs["BCFTools"]=0;                 rna_programs_order+=("BCFTools")
    rna_programs["SnpEff"]=0;                   rna_programs_order+=("SnpEff")

    echo "State before choosing: ${rna_programs["FastQC"]}, $rna_programs[]"

    # TODO: For each program, add its function/script to run as an arra i.e
    # declare -A rna_program_script;                            declare -a rna_script_order;
    # rna_program_script["FastQC"]="./scripts/fastqc.bash";  rna_script_order+=("FastQC")

    # From ./utils.bash
    category_chooser rna_categories rna_categories_order rna_programs rna_programs_order
}


## Ask for CORE and RAM 
NUM_CORES=$(get_core_count)
MAX_RAM=$(get_available_ram)

RNA_Questionnaire

echo -e "----------------\e[1mRunning Workflow: \e[0m-------------------"

#### Quality Control ####

PROCESSED_READS_DIR="$INPUT_GENOME_PATH"

## FastQC
if [ "${rna_programs["FastQC"]}" -eq 1 ]; then
  echo "FASTQC works"
   ./scripts/fastqc_subscript.bash "$INPUT_GENOME_PATH"\
                                   "$OUTPUT_DIR"
fi

## fastp
if [ "${rna_programs[fastp]}" -eq 1 ]; then
   ./scripts/fastp_subscript.bash "$INPUT_GENOME_PATH"\
                                  "$OUTPUT_DIR"

    PROCESSED_READS_DIR="${OUTPUT_DIR}/fastp_output"
fi

## TODO: Pause after FastQC, since we need to determine how much we want to trim, so we can ask whether to continue
## TODO: Is this needed? Things like trimgalore and fastp do the trimming for you, so you don't need to pause
echo -e "\e[1m Read analysis complete, check quality and \e[0m"
read -p "Press enter to continue..."



## Trimming reads



#### Assembly ####
if [ "${rna_programs[minimap2]}" -eq 1 ]; then
  #if [ "${programs[fastp]}" -eq 1 ]; then
    ./scripts/minimap2_subscript.bash "$HUMAN_REFERENCE_PATH"\
                                      "$PROCESSED_READS_DIR"\
		                                  "$OUTPUT_DIR"
  #else
    #echo "When using an assembler, you should use fastp first to ensure
    #      only high quality reads being used."
  #fi
fi



## DNA Pipeline
#echo "DNA pipeline starting"
#./scripts/dna_pipeline.sh "$INPUT_GENOME_PATH" "$OUTPUT_DIR" "$HUMAN_REFERENCE_PATH"


