

process LearnReadOrientationModel {
    // This requires all the f1r2 files from scattered analysis
    input:
    val(all_f1r2)

    output:
    path("*.tar.gz")

    script:

    def input_args = ""
    for (f1r2 in all_f1r2) {
        input_args += "-I ${f1r2} "
    }
    """
    gatk LearnReadOrientationModel -I ${input_args} -O read-orientation-model.tar.gz
    """
}

process MergeMutectStats {
    input:
    val(all_stats)

    output:
    path("*.stats")

    script:

    def input_args = ""
    for (stats in all_stats) {
        input_args += "-stats ${stats} "
    }
    """
    gatk MergeMutectStats \
        ${input_args}
        -O merged.stats
    """
}





process VariantFiltration {
    label 'variant_calling'
    publishDir "${params.output_dir}/vcf/filtered/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    publishDir "${params.output_dir}/vcf/filtered/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
    publishDir "${params.output_dir}/vcf/filtered/normal", mode: 'copy', overwrite: true, pattern: "*snc*.vcf"

    input:
        tuple val(sample_id), path(vcf)
        path ref
        path ref_fai
        path ref_dict

    output:
        path "*.vcf"

    script:
    // Layout: [Filter Expression, Filtername]
    def filter_options = [
        ["FS > 20", "FS20"]
        ["QUAL > 20", "FS20"]
        ]
    def filtering_args = ""
    filter_options.each { expr, name ->
        filtering_args += "--genotype-filter-expression \"${expr}\" --genotype-filter-name \"${name}\" "
    }
    println ${filtering_args}
    // TODO: Integrate the filtering args into the command block
    """
    gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4' VariantFiltration \
        -R \$PWD/${ref} \
        -I \$PWD/${vcf} \
        -O \$PWD/${vcf.simpleName}_filtered.vcf \
        ${filtering_args}
    """
}

workflow VARIANT_PROCESSING {
    take:
        vcfs // Sample ID + vcfs
        ref_auxillary
    main:
        haplotype_vcf = VariantFiltering(bam_split_n,
                                ref_auxillary.out.fai,
                                ref_auxillary.out.dict,
                                ref_auxillary.out.ref)
    emit:
        haplotype_vcf
}
