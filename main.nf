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

    // Input Reads; {gz,bz2} needed since sometimes naming is bad i.e .gz.normalised.vcf != read file
    // TODO: Either format unpaired to be in the same format as paired i.e tuple val(x), path(read)
    //       OR let each process accept a val(x) and then on each process determine whether its a tuple
    if (params.paired) {
        // TODO: placing the "." inside i.e  [.gz|.bz2] causes it to not function?
        // FIXME: the input_dir/**/files_here makes it so it has to be in sub-dirs. Make it also take ones that are directly
        // in the dir i.e input_dir/files_here
        INPUT_READS = Channel.fromFilePairs("${params.input_dir}/**/*{${params.r1_pattern},${params.r2_pattern}}*.f*q.[gz|bz2]?",
                                            type: 'file',
                                            maxDepth: 5)
    } else {
        // TODO: What to put as the second element in the list for the tuple i.e [name_id, [read, _]]? Or we can just leave it fully alone?
        INPUT_READS = Channel.fromPath("${params.input_dir}/*.f*q.[gz|bz2]?", type: 'file', maxDepth: 5)
                    .map(read -> tuple(read.simpleName, read))
    }

    // Temporary filter, move them out of dir when permissions allow, it removes all DNA-fastq's
    // INPUT_READS.filter { !it[0].startsWith("FO")  }
    // INPUT_READS.view()


    QC_READS = QUALITYCONTROL(INPUT_READS)
    MAPPING(QC_READS, REFERENCE)
    VARIANT_CALLING(MAPPING.out, REFERENCE) // Places files in output folder

    // If assembly wanted
    // ASSEMBLY(INPUT_READS, MAPPING.out, REFERENCE) 
    // VARIANT_CALLING(MAPPING.out.concat(ASSEMBLY.out), REFERENCE) // Places files in output folder

}
