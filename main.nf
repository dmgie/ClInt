#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QUALITYCONTROL } from './preprocess'
include { MAPPING } from './mapping'
include { ASSEMBLY } from './assembly'
include { VARIANT_CALLING } from './variant_calling'

// Set up some params

REFERENCE = file(params.reference_file); // Require as file() can't be empty
ANNOTATION = file(params.gff_file); // Require as file() can't be empty

def docker_current_dir() {
    REF_FULLPATH = "realpath ${params.reference_file}".execute().text.trim()
    "docker run -v $PWD:$PWD -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
}

def CHECKPARAMS() {
    println "Checking parameters..."
    if (params.input_reads_path == '') {
        error "ERROR: input_reads_path is not set"
    } else if (params.output_dir == '') {
        error "ERROR: output_dir is not set"
    } else if (params.reference_file == '') {
        error "ERROR: reference_file is not set"
    } else if (params.gff_file == '') {
        error "ERROR: gff_file is not set"
    } else {
        println "All parameters are set"
        println "   User ID:\t ${params.USER_ID}"
        println "   Group ID:\t ${params.GROUP_ID}"
    }
}

workflow {
    CHECKPARAMS()
 
    READS = QUALITYCONTROL(Channel.fromPath(params.input_reads))
    MAPPING(REFERENCE, READS) 
    ASSEMBLY(REFERENCE, READS, MAPPING.out) 
    VARIANT_CALLING(MAPPING.out.concat(ASSEMBLY.out), REFERENCE) // Places files in output folder
    // VARIANT_CALLING.out.view()
}