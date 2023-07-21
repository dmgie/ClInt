#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Set up params
params.input_reads_path = '.'; // Default to the current directory
params.input_reads = file("${params.input_reads_path}/*.fastq.gz");
params.output_path = file("./workflow_output"); // Default to subdir in current dir
params.reference_file = ''; 
params.gff_file = '';
params.num_threads = 4;
params.max_memory = '8.GB';

REFERENCE = file(params.reference_file); // Require as file() can't be empty

def docker_current_dir(inp) {
    // Simple function to return the docker command for the current process (as the PWD is different for each process)

    // This is a bit of a hack, since the docker container would need to get the reference file from the host,
    // but inputs into processes are symlinked to their original location, so the docker container can't see them

    // That is the reason why we are mounting REF_FULLPATH into the container, so that it can see the reference file
    // from the symlink that is given as input
    REF_FULLPATH = "realpath ${params.reference_file}".execute().text.trim()
    "docker run -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
}

def CHECKPARAMS() {
    // Simple function to check if all parameters are set

    println "Checking parameters..."
    if (params.input_reads_path == '') {
        error "ERROR: input_reads_path is not set"
    } else if (params.output_path == '') {
        error "ERROR: output_path is not set"
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

// --------------------Program-specific Processes-------------------------

process HISAT_BUILD {
    input:
        path ref_file
    output:
        path "ref_idx*.ht2"

    script:
    def NUM_THREADS = 4
    """
    # hisat2-build -p $NUM_THREADS  ${ref_file} ref_idx
    echo "Running hisat2-build"
    touch ref_idx.1.ht2
    touch ref_idx.2.ht2
    touch ref_idx.3.ht2
    ls -lah
    """
}

process HISAT2 {
    input:
        path ref_idx
        path reads

    output:
        path "aligned_${read}.bam"

    script:
    def NUM_THREADS = 4
    def index_name = "ref_idx"
    """
    # hisat2 -p $NUM_THREADS --rg-id aligned_${reads} --rg SM:None --rg LB:None --rg PL:Illumina -x ${index_name} -x ${ref_idx} -U ${reads} -S aligned_${reads}.sam
    # samtools sort -@ $NUM_THREADS -o aligned_${reads}.bam aligned_${reads}.sam
    ls -lah
    """
}

process STAR_BUILD {
    input:
        path ref_file
        path annotation
    output:
        path "*"

    script:
    def NUM_THREADS = 4
    def extension = annotation.extension
    def feature = "gene"
    """
    if [[ ${extension} == *.gtf ]]
        then
            echo "GTF file detected"
            echo "STAR --runThreadN $NUM_THREADS \
            --runMode genomeGenerate \
            --genomeDir . \
            --genomeFastaFiles $ref_file \
            --sjdbGTFfile $annotation \
            --sjdbOverhang \$((READ_LENGTH - 1))"
        elif [[ ${extension} == *.gff ]]
            then
            echo "GFF file detected"
            echo "STAR --runThreadN $NUM_THREADS \
            --runMode genomeGenerate \
            --genomeDir . \
            --genomeFastaFiles $ref_file \
            --sjdbGTFfile $annotation \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbGTFfeatureExon $feature \
            --sjdbOverhang \$((READ_LENGTH - 1))"
    else
        echo "ERROR: Annotations file must be either gtf or gff"
            exit 1
    fi
    """
}





// -----------------------------------------------------------------------

process QUALITYCONTROL {
    // if $PUBLISH_DIRECTORIES == 1 {
    //     publishDir "${params.output_path}/QC" // mode: 'copy',
    // }
    maxForks 2
    
    input:
        path read

    output:
        // path "trimmed_*"
        stdout emit: verb


    //fastp -i ${read} -o trimmed_${read} -j /dev/null -h /dev/null
    script:
    """
    echo "Working on ${read}"
    """
}


process ASSEMBLY {
    // if $PUBLISH_DIRECTORIES == 1 {
    //     publishDir "${params.output_path}/ASSEMBLY" // mode: 'copy',
    // }

    input:
        each method
        path qc_reads

    output:
        stdout
        path "${method}_*"

    def docker = docker_current_dir("${PWD}")
    

    script:
    """
    echo ${method} 
    """
}

process REFERENCE_HELP_FILES {
    // This process mainly relates to the necessary side-files
    // such as the .fai and .dict files for the reference genome
    // as these are required in some of the processes

    input: 
        path ref_file

    output:
        val "${ref_file.simpleName}.fai" // val used to allow for multiple consumption
        val "${ref_file.simpleName}.dict" // val used to allow for multiple consumption
        stdout emit: verbo


    def docker = docker_current_dir("ref_file.fasta")
    script:
    """
    echo "Running samtools faidx and docker"
    #samtools faidx ${ref_file} && mv ${ref_file}.fai ${ref_file.simpleName}.fai
    #${docker} broadinstitute/gatk gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file.simpleName}.dict
    ls -lah
    """
}

workflow RNASEQ_QUANT_PIPE {
    take:
      transcriptome_ch
      read_pairs_ch
    main:
      transcriptome_ch = channel.fromPath(params.transcriptome)
      INDEX(transcriptome_ch)
      QUANT(INDEX.out,read_pairs_ch)
}

workflow {
    CHECKPARAMS()
    
    // (ref_fai, ref_dict) = REFERENCE_HELP_FILES(Channel.fromPath(REFERENCE))
    // REFERENCE_HELP_FILES.out.verbo.view()

    // QC_ch = QUALITYCONTROL(Channel.fromPath(params.input_reads))
    // QUALITYCONTROL.out.verb.view()
    // ASSEMBLY_ch = ASSEMBLY(ASSEMBLY_METHODS,QC_ch)

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
