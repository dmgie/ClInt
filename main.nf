#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QUALITYCONTROL } from './preprocess'
include { MAPPING } from './mapping'
include { ASSEMBLY } from './assembly'
include { VARIANT_CALLING } from './variant_calling'

// Set up some params

REFERENCE = file(params.reference_file); // Require as file() can't be empty
ANNOTATION = file(params.gff_file); // Require as file() can't be empty


def CHECKPARAMS() {
    println "Checking parameters..."
    if (params.input_dir == '') {
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

    INPUT_READS = Channel.fromPath("${params.input_dir}/*.f*q.gz"); // Matches regex
    if (paired) {
        INPUT_READS = Channel.fromPath(${params.metadata}) | splitCsv(header:true) \
            | map { row -> tuple(row.sample_name,path(row.r1_path),path(row.r2_path) }
            | view { it }
    }



    // QC_READS = QUALITYCONTROL(INPUT_READS)
    // MAPPING(REFERENCE, QC_READS)
    // VARIANT_CALLING(MAPPING.out, REFERENCE) // Places files in output folder

    // If assembly
    // ASSEMBLY(REFERENCE, READS, MAPPING.out) 
    // VARIANT_CALLING(MAPPING.out.concat(ASSEMBLY.out), REFERENCE) // Places files in output folder

    VARIANT_CALLING(MAPPING.out, REFERENCE) // Places files in output folder

}
