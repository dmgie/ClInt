#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Set up some params
params.input_reads = file("${params.input_reads_path}/*.f*q.gz");
params.output_dir = file("./workflow_output"); // Default to subdir in current dir
REFERENCE = file(params.reference_file); // Require as file() can't be empty
ANNOTATION = file(params.gff_file); // Require as file() can't be empty

def docker_current_dir() {
    // Simple function to return the docker command for the current process (as the PWD is different for each process)

    // This is a bit of a hack, since the docker container would need to get the reference file from the host,
    // but inputs into processes are symlinked to their original location, so the docker container can't see them

    // That is the reason why we are mounting REF_FULLPATH into the container, so that it can see the reference file
    // from the symlink that is given as input
    REF_FULLPATH = "realpath ${params.reference_file}".execute().text.trim()
    // "docker run -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
    "docker run -v $PWD:$PWD -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
}

def CHECKPARAMS() {
    // Simple function to check if all parameters are set

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
        // println "   Reference file:\t ${REF_FULLPATH}"
        // println "   Docker command:\t ${DOCKER_COMMAND}"
    }
}

process QUALITYCONTROL {
    maxForks 5
    // if ($PUBLISH_DIRECTORIES == 1) {
    //     publishDir "${params.output_dir}/QC" // mode: 'copy',
    // }
    
    input:
        path read

    output:
        path "trimmed_*"
        // stdout emit: verb


    script:
    """
    echo "Working on ${read}"
    fastp -i ${read} -o trimmed_${read} -j /dev/null -h /dev/null
    """
}


workflow {
    CHECKPARAMS()
    READS = QUALITYCONTROL(Channel.fromPath(params.input_reads))
    

    MAPPING(REFERENCE, READS) // Will output bams
    ASSEMBLY(READS, MAPPING.out) // Will output fasta

    // FIXME: Give GATK both the MAPPING.out and ASSEMBLY.out (both bam files)
    GATK(MAPPING.out.concat(ASSEMBLY.out), ref_fai, ref_dict, REFERENCE)
    GATK.out.view()
}

/*
    * Notes on nextflow:
        * If you wanna see stdout, its "PROCESS.out.view()"
        * If you want to define a variable within a process, use "def" after the input/output section
        * the 'val' keyword can be used an ulimited number of times, in comparison to channles which can only be used once 
            (and are consumed)
            Actually, it seems that DSL2 allows for re-use of channels
        * For passing i.e ".ht2" files all at once, we can do a collect() on the output of the process which produces these index files

*/
