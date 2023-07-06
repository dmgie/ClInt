

INPUT_FILE=$1
SPADES_OUTPUT="${OUTPUT_DIR}/spades_output"


# rnaspades.py/ spades.py --rna
## --ss-fr = first read in stranded pair corresponsds to the actual gene strand <--- used for single-end
## --ss-rf = first read in stranded pair corresponds to the reverse strand
$SPADES_PATH --rna --ss-fr -o $SPADES_OUTPUT -s $INPUT_FILE -t $NUM_CORES -m $MAX_RAM

# The output is transcripts.fasta
