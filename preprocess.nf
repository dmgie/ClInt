process FASTP {
    label 'preprocess'
    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("*.f*q.*")

    script:
    def (read1,read2) = [reads[0], reads[1]]
    def arguments = params.paired ? "--in1 ${read1} --in2 ${read2} --out1 trimmed_${read1} --out2 trimmed_${read2}" : 
                                    "-i ${read1} -o trimmed_${read1}"
    """
    echo "Working on ${reads} using ${arguments}"
    fastp ${arguments} --threads ${task.cpus} -j /dev/null -h /dev/null
    """

    stub:
    """
    touch trimmed_${reads}
    """
}

workflow QUALITYCONTROL {
    take:
        reads
    main:
        trimmed_reads = FASTP(reads)
    emit:
        trimmed_reads
}
