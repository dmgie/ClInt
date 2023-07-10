#!/usr/bin/env bash

REF=$1
READS=$2
ANNOTATIONS=$2
READ_LENGTH=$3
OUTPUT_DIR=${WORKING_DIR}/STAR_output
SINGLE_OR_PAIRED=$4 # 0 for single, 1 for paired
mkdir $OUTPUT_DIR

# Make sure its empty
rm -rf $OUTPUT_DIR/*

# Create indexes
# Check if ANNOTATIONS is gtf or gff
if [[ $ANNOTATIONS == *.gtf ]]
then
    echo "GTF file detected"
    STAR --runThreadN $NUM_CORES \
        --runMode genomeGenerate \
        --genomeDir $OUTPUT_DIR \
        --genomeFastaFiles $REF \
        --sjdbGTFfile $ANNOTATIONS \
        --sjdbOverhang $((READ_LENGTH - 1))
elif [[ $ANNOTATIONS == *.gff ]]
then
    echo "GFF file detected"
    STAR --runThreadN $NUM_CORES \
        --runMode genomeGenerate \
        --genomeDir $OUTPUT_DIR \
        --genomeFastaFiles $REF \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbOverhang $((READ_LENGTH - 1))
else
    echo "ERROR: Annotations file must be either gtf or gff"
    exit 1
fi


# Mapping
if [[ $SINGLE_OR_PAIRED -eq 0 ]]; then
    # Make ${READS}/*.fastq.gz into an array
    READS_ARRAY=(${READS}/*.f*q.gz)

    # Create a comma separated list of the reads (filepaths)
    READS=$(echo ${READS_ARRAY[@]} | sed 's/ /,/g')


    # According to the docs, this can be used for the SAM files READ_GROUP-related
    # tags and improtant for GATK later
    # Create another " , " separated list of the reads, prepending ID: to each element
    READS_ARRAY_BASENAME=(${READS_ARRAY[@]##*/}) # Remove path, only filename
    READS_ARRAY_BASENAME=(${READS_ARRAY_BASENAME[@]%%.*}) # Remove extensions
    READS_ARRAY_BASENAME=(${READS_ARRAY_BASENAME[@]/#/ID:}) # Prepend ID:

    # Add $OTHER_ATTRIBUTES to each element, (use ";" for separation for now
    # so its easier to replace later with "sed", otherwise all spaces turn into
    # commas)
    OTHER_ATTRIBUTES="PL:ILLUMINA;PU:UNKNOWN;LB:UNKNOWN;SM:UNKNOWN"
    for ((i=0; i<${#READS_ARRAY_BASENAME[@]}; i++)); do
        READS_ARRAY_BASENAME[$i]="${READS_ARRAY_BASENAME[$i]};$OTHER_ATTRIBUTES"
    done

    SAM_ATTRIBUTE_LINE=$(echo ${READS_ARRAY_BASENAME[@]} | sed 's/ / , /g' | sed 's/;/ /g') # Make it " , "

    # NOTE: The above ensures  that the SAM ReadGroup attribute line is in the correct format
    # Where the filename is the ID, and PL/PU/LB/SM are also given
    # This is formatted according to the STAR manual in a way so that an entire folder
    # of (single-end) reads can be given, and each of them have their RG line

    STAR --runThreadN $NUM_CORES \
        --genomeDir $OUTPUT_DIR \
        --readFilesIn $READS \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUTPUT_DIR}/ \
        --outSAMattrRGline $SAM_ATTRIBUTE_LINE \
else
    # TO BE IMPLEMENTED -> Look at the STAR manual on the github page
    # STAR --runThreadN $NUM_CORES \
        #     --genomeDir $OUTPUT_DIR \
        #     --readFilesIn ${READS}/*.fastq.gz \
        #     --readFilesCommand zcat \
        #     --outSAMtype BAM SortedByCoordinate \
        #     --outFileNamePrefix ${OUTPUT_DIR}/
    echo "Wow"
fi
