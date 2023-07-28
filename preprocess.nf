process FASTP {
    maxForks 5
    input:
        path reads
    output:
        path "trimmed_*"

    script:
    """
    echo "Working on ${reads}"
    fastp -i ${reads} -o trimmed_${reads} -j /dev/null -h /dev/null
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
