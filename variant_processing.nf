

process MergeOrientationModel {
    label 'variant_calling'
    label 'falliable'
    // This requires all the f1r2 files from scattered analysis
    input:
    tuple val(sample_id), path(all_f1r2) // this takes all paths to f1r2 for one sample ID

    output:
    tuple val(sample_id), path("*.tar.gz")

    script:

    def input_args = ""
    for (f1r2 in all_f1r2) {
        input_args += "-I ${f1r2} "
    }

    """
    gatk LearnReadOrientationModel ${input_args} -O ${sample_id}_read-orientation-model.tar.gz
    """
}

process MergeMutectStats {
    label 'variant_calling'
    input:
    tuple val(sample_id), path(all_stats)

    output:
    tuple val(sample_id), path("*.stats")

    script:

    def input_args = ""
    for (stats in all_stats) {
        input_args += "-stats ${stats} "
    }
    """
    gatk MergeMutectStats \
        ${input_args} \
        -O ${sample_id}_merged.stats
    """
}

process FilterMutect {
    label 'variant_calling'
    label 'falliable'
    publishDir "${params.output_dir}/vcf/filtered/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.vcf"

    input:
    tuple val(sample_id), path(req_files)
    path fai
    path dict
    path ref

    output:
    tuple val(sample_id), path("*.vcf")

    script:
    def (vcf, read_orient, stats) = req_files
    """
    gatk FilterMutectCalls \
        -R ${ref} \
        -V ${vcf} \
        --stats ${stats} \
        --ob-priors ${read_orient} \
        -O filtered_${sample_id}.vcf
    """
}

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
