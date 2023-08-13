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

process CSV_CREATION {
    publishDir "${baseDir}", mode: 'copy', overwrite: true

    // Create a CSV file detailing the sample_name and the file paths to the read(s)
    // This is depending on whether the "paired" flag was enabled or disabled when launching the command

    input: 
        path read_dir
    output:
        path "*.csv"
        stdout emit: std

    shell:
    r1 = "_R1" 
    r2 = "_R2"
    depth = 1

    shell:
    if (params.paired) {
        '''#!/usr/bin/env bash
        # Get the input arguments
        input_path=$(realpath !{read_dir})
        input_files=($(find $input_path -name "*.f*q.*"))
        r1_pattern=!{r1}
        r2_pattern=!{r2}

        # Remove non-paired files i.e files which do not have both *r1* and *r2* patterns
        for file in ${input_files[@]}; do
            if [[ $file =~ $r1_pattern ]]; then
                # check if the corresponding r2 file exists
                r2_file=${file/$r1_pattern/$r2_pattern}
                if [[ -f $r2_file ]]; then
                    sample_name=${file/$r1_pattern/}
                    sample_name=${sample_name%_}
                    sample_list+=($sample_name)
                    r1_files+=($file)
                    r2_files+=($r2_file)
                fi
            fi
        done

        # Remove paths and extension from the sample names
        for i in ${!sample_list[@]}; do
            sample_list[$i]=${sample_list[$i]##*/}
            sample_list[$i]=${sample_list[$i]%%.*}
        done

        # Create the csv file with the sample names and paths
        echo "sample_name,r1_path,r2_path" > clint_metadata.csv # Header
        for i in ${!sample_list[@]}; do
            echo "${sample_list[$i]},${r1_files[$i]},${r2_files[$i]}" >> clint_metadata.csv
        done

        cat clint_metadata.csv
        '''
    } else 
        '''
        input_path=$(realpath !{read_dir})
        input_files=($(find $input_path -maxdepth !{depth} -type f -regex '.*f*q\\(.gz\\|.bz2\\)?'))
        files_path=($(find $input_path -maxdepth !{depth} -type f -regex '.*f*q\\(.gz\\|.bz2\\)?'))

        echo ${input_files[@]}

        for i in ${!input_files[@]}; do
            sample_name[$i]=${input_files[$i]##*/}
            sample_name[$i]=${sample_name[$i]%%.*}
        done
        echo "sample_name,r1_path,r2_path" > clint_metadata.csv # Header
        for i in ${!sample_name[@]}; do
            echo "${sample_name[$i]}, ${files_path[$i]}" >> clint_metadata.csv
        done

        cat clint_metadata.csv
        '''
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
