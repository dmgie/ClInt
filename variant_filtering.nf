process VariantFiltration {
    label 'variant_calling'
    publishDir "${params.output_dir}/vcf/filtered/", mode: 'copy', overwrite: true, pattern: "*snc*.vcf"
    // publishDir "${params.output_dir}/vcf/filtered/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    // publishDir "${params.output_dir}/vcf/filtered/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
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
