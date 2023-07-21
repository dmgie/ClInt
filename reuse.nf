workflow MAPPING {
    take: 
        ref_file
        reads
    main:
        mapping_method = params.mapping.toLowerCase()
        if (mapping_method == "hisat2") {
            ref_idx = HISAT_BUILD(ref_file).collect()
            aligned_bams = HISAT2(ref_idx, reads)
            sorted_bams = SAMTOOLS_SORT(aligned_bams)
        } else if (mapping_method == "star") {
            ref_idx = STAR_BUILD(ref_file, ANNOTATION).collect()
            sorted_bams = STAR(ref_idx, reads) // The command itself aligns the bams
        } else {
            println "ERROR: Mapping method not recognised"
        }
    emit:
        sorted_bams // output bam files (sorted)
}
