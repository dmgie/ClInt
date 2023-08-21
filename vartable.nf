#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VarTable {
    input:
        path(vcfs)

    output:
        path("*.tsv")

    script:
    """
    python3 -script- -arguments-
    """
}

workflow {
    take:
        vcfs

    main:
        tables = VarTabl(vcfs)

    emit:
        tables
}