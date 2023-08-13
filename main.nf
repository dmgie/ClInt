#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QUALITYCONTROL } from './preprocess'
include { MAPPING } from './mapping'
include { ASSEMBLY } from './assembly'
include { VARIANT_CALLING } from './variant_calling'

// Set up some params

REFERENCE = file(params.reference_file); // Require as file() can't be empty
ANNOTATION = file(params.gff_file); // Require as file() can't be empty
INPUT_DIR = file(params.input_dir);


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

    // TODO: Call CSV_CREATION depending on if it was paired data or not
    // TODO: If metadata already given, don't run CSV creation, this allows for adding metadata
    // TODO: Somehow enable giving strandedness as an option and adjusting the methods based on that
    //      Could do: If during CSV read = strandedness column empty, we read from nextflow.config (i.e variable with default RF or FR)
    // TODO: Should input_dir fastq.gz search be recursive to look at all subfolders?
    // TODO: Build up commands inside the program processes i.e "BASE_COMMANDS", "READ_COMMANDS" etc all depending if paired or unpaired
    //       Also use ternary's if necessary; Example: https://stackoverflow.com/questions/68442177/nextflow-change-part-of-the-script-basing-on-a-parameter
    // TODO: Maybe just use `Channel.fromFilePairs` for READS_{1,2}
    // METADATA = CSV_CREATION(INPUT_DIR)
    // Setup paired/unpaired reads Channel
    // INPUT_READS = METADATA.splitCsv(header: true)
    //     .map { row -> tuple(row.sample_name,file(row.r1_path),file(row.r2_path)) }
    //     .view { it }
    // METADATA.std.view()

    // Input Reads; {gz,bz2} needed since sometimes naming is bad i.e .gz.normalised.vcf != read file
    if (params.paired) {
        // TODO: placing the "." inside i.e  [.gz|.bz2] causes it to not function?
        INPUT_READS = Channel.fromFilePairs("${params.input_dir}/*{${params.r1_pattern},${params.r2_pattern}}*.f*q.[gz|bz2]?",
                                            type: 'file',
                                            maxDepth: 5)
    } else {
        INPUT_READS = Channel.fromPath("${params.input_dir}/*.f*q.[gz|bz2]?", type: 'file', maxDepth: 5)
    }

    INPUT_READS.view()


    // QC_READS = QUALITYCONTROL(INPUT_READS)
    // MAPPING(QC_READS, REFERENCE)
    // VARIANT_CALLING(MAPPING.out, REFERENCE) // Places files in output folder

    // If assembly wanted
    // ASSEMBLY(INPUT_READS, MAPPING.out, REFERENCE) 
    // VARIANT_CALLING(MAPPING.out.concat(ASSEMBLY.out), REFERENCE) // Places files in output folder

}
