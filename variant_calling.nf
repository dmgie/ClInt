process MarkDuplicates {
    // NOTE: We an add --REMOVE_DUPLICATES=true to remove duplicates from the final BAM file
    //       intead of just switching the flag for that read
    label 'variant_calling'
    publishDir "${params.output_dir}/deduped_bam/", mode: 'copy', overwrite: true
    input:
        path aligned_bam

    output:
        path "dedup_${aligned_bam}"


    script:
    """
    echo "Working on ${aligned_bam}"
    gatk MarkDuplicates -I \$PWD/${aligned_bam} -O \$PWD/dedup_${aligned_bam} -M \$PWD/dedup_${aligned_bam}.metrics
    """

    stub:
    """
    touch dedup_${aligned_bam}
    """
}

process SplitNCigarReads {
    label 'variant_calling'
    input:
        path aligned_bam
        path ref_fai
        path ref_dict
        path ref

    output:
        path "snc_${aligned_bam}"
        // stdout emit: temp
    

    script:
    """
    echo "Working on ${aligned_bam}"
    gatk SplitNCigarReads -R \$PWD/${ref} -I \$PWD/${aligned_bam} -O \$PWD/snc_${aligned_bam}
    """

    stub:
    """
    touch snc_${aligned_bam}
    """
}

process VariantFiltering {
    label 'variant_calling'
    publishDir "${params.output_dir}/filtered_vcf/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    publishDir "${params.output_dir}/filtered_vcf/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
    publishDir "${params.output_dir}/filtered_vcf/normal", mode: 'copy', overwrite: true, pattern: "*snc_trimmed*.vcf"

    input:
        path vcf
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
    // println ${filtering_args}
    // TODO: Integrate the filtering args into the command block
    """
    gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4' \
        -R \$PWD/${ref} \
        -I \$PWD/${vcf} \
        -O \$PWD/${vcf.simpleName}_filtered.vcf \
        ${filtering_args}
    """
}

process HaplotypeCaller {
    label 'variant_calling'
    publishDir "${params.output_dir}/haplotype_vcf/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    publishDir "${params.output_dir}/haplotype_vcf/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
    // FIXME: Maybe fix this so that non-assembled ones have their own name? But currently based upon that
    //        only the non-assembled ones don't have a method between "snc" and "trimmed"
    publishDir "${params.output_dir}/haplotype_vcf/normal", mode: 'copy', overwrite: true, pattern: "*snc_trimmed*.vcf"
    input:
        path split_bam
        path ref_fai
        path ref_dict
        path ref

    output:
        path "haplotype_*.vcf"


    // TODO: Split by chromosome. Either have a for loop in the command block (each of them uses chromosome interval)
    // or then do it on the process-level
    script:
    """
    echo "Working on ${split_bam}"
    samtools index ${split_bam}
    gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
    --pair-hmm-implementation FASTEST_AVAILABLE \
    --smith-waterman FASTEST_AVAILABLE \
    -R \$PWD/${ref} -I \$PWD/${split_bam} -O \$PWD/haplotype_${split_bam}.vcf
    #touch haplotype_${split_bam.simpleName}.vcf
    ls -lah
    """

    stub:
    """
    touch haplotype_${split_bam.simpleName}.vcf
    """
}


process REFERENCE_HELP_FILES {
    // This process mainly relates to the necessary side-files
    // such as the .fai and .dict files for the reference genome
    // as these are required in some of the processes

    input: 
        path ref_file

    output:
        path "${ref_file}.fai" 
        path "${ref_file.baseName}.dict" 
        // stdout emit: verbo


    script:
    """
    echo "Running samtools faidx and docker"
    samtools faidx ${ref_file}
    gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file.baseName}.dict
    ls -lah
    """

    stub:
    """
    touch ${ref_file}.fai
    touch ${ref_file.baseName}.dict
    """
}

workflow VARIANT_CALLING {
    // Run SplitNCigarReads and HaplotypeCaller

    take:
        bam
        ref
    main:
        // TODO: Combine ref_fai, ref_dict and ref into one thing
        REFERENCE_HELP_FILES(ref)
        (ref_fai, ref_dict) = REFERENCE_HELP_FILES.out
        split_bam = SplitNCigarReads(bam, ref_fai, ref_dict, ref) | MarkDuplicates
        // SplitNCigarReads.out.view()
        haplotype_vcf = HaplotypeCaller(split_bam, ref_fai, ref_dict, ref)
    emit:
        // split_bam
        haplotype_vcf
}

