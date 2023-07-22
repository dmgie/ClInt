def docker_current_dir() {
    REF_FULLPATH = "realpath ${params.reference_file}".execute().text.trim()
    "docker run -v $PWD:$PWD -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
}

process MarkDuplicates {
    // NOTE: We an add --REMOVE_DUPLICATES=true to remove duplicates from the final BAM file
    //       intead of just switching the flag for that read
    maxForks 5
    publishDir "${params.output_dir}/deduped_bam/", mode: 'copy', overwrite: true
    input:
        path aligned_bam

    output:
        path "dedup_${aligned_bam}"

    def docker = docker_current_dir()

    script:
    """
    echo "Working on ${aligned_bam}"
    ${docker} broadinstitute/gatk gatk MarkDuplicates -I \$PWD/${aligned_bam} -O \$PWD/dedup_${aligned_bam} -M \$PWD/dedup_${aligned_bam}.metrics
    """

    stub:
    """
    touch dedup_${aligned_bam}
    """
}

process SplitNCigarReads {
    // TODO: Change aligned_bam to dedup_bam, and change the output to dedup_${aligned_bam}, because we wanna remove duplicates first
    maxForks 8
    input:
        path aligned_bam
        path ref_fai
        path ref_dict
        path ref

    output:
        path "snc_${aligned_bam}"
        // stdout emit: temp
    
    def docker = docker_current_dir()

    script:
    def NUM_THREADS = 4
    """
    echo "Working on ${aligned_bam}"

    ${docker} broadinstitute/gatk gatk SplitNCigarReads -R \$PWD/${ref} -I \$PWD/${aligned_bam} -O \$PWD/snc_${aligned_bam}

    # ${docker} broadinstitute/gatk bash -c "cd \$PWD/; ls -lah ../../"
    # ${docker} broadinstitute/gatk gatk SplitNCigarReads -R \$PWD/${ref} -I \$PWD/${aligned_bam} -O \$PWD/snc_${aligned_bam}
    # echo "gatk SplitNCigarReads -R ${ref} -I ${aligned_bam} -O split_${aligned_bam}"
    # touch snc_${aligned_bam}
    # ls -lah
    """

    stub:
    """
    touch snc_${aligned_bam}
    """
}

process HaplotypeCaller {
    maxForks 8
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

    def docker = docker_current_dir()

    script:
    def NUM_THREADS = 4
    """
    echo "Working on ${split_bam}"
    samtools index ${split_bam}
    $docker broadinstitute/gatk gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
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


    def docker = docker_current_dir()
    script:
    """
    echo "Running samtools faidx and docker"
    samtools faidx ${ref_file}
    ${docker} broadinstitute/gatk gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file.baseName}.dict
    #${docker} broadinstitute/gatk gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file}.dict
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
        REFERENCE_HELP_FILES(ref)
        (ref_fai, ref_dict) = REFERENCE_HELP_FILES.out
        split_bam = SplitNCigarReads(bam, ref_fai, ref_dict, ref) | MarkDuplicates
        // SplitNCigarReads.out.view()
        haplotype_vcf = HaplotypeCaller(split_bam, ref_fai, ref_dict, ref)
    emit:
        // split_bam
        haplotype_vcf
}
