include { SAMTOOLS_INDEX } from './mapping'
include { SAMTOOLS_SORT } from './mapping'

process REF_AUXILLARY {
    // This process mainly relates to the necessary side-files
    // such as the .fai and .dict files for the reference genome
    // as these are required in some of the processes

    input:
        path ref_file

    output:
        path "${ref_file}.fai", emit: fai
        path "${ref_file.baseName}.dict", emit: dict
        path ref_file, emit: ref
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

process MarkDuplicates {
    // NOTE: We an add --REMOVE_DUPLICATES=true to remove duplicates from the final BAM file
    //       intead of just switching the flag for that read
    label 'variant_calling'
    publishDir "${params.output_dir}/deduped_bam/", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), path(aligned_bam)

    output:
        tuple val(sample_id), path("dedup_*.bam")

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
        tuple val(sample_id), path(bam), path(bai)
        path ref_fai
        path ref_dict
        path ref
        each chr_interval // [1,2], use "var" if using entire chromosome list [1..21,X,Y]

    output:
        tuple val(sample_id), path("snc_*.bam")
        // stdout emit: temp

    // TODO: Parallelise using interval list, in pairs of 2
    // TODO :This is currently done as a single process/job submission, maybe split it into multiple jobs? Then collect and converge
    // FIXME: This might just be overwriting at each interval, so do the looping somewhere else
    script:
    def name = "snc_${bam.simpleName}."
    def interval_args = ""
    for (chr in chr_interval) {
        interval_args += " -L ${chr}"
        name += "${chr}_"
    }
    println "Processing SplitNCigarReads for ${interval_args} for ${sample_id}"

    """
    echo "Working on ${bam}"
    gatk SplitNCigarReads -R ${ref} -I ${bam} -O ${name}.bam ${interval_args}
    """

    stub:
    """
    touch snc_${bam_bai}
    """
}

process Mutect2 {
    label 'variant_calling'
    publishDir "${params.output_dir}/vcf/intermediate/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    publishDir "${params.output_dir}/vcf/intermediate/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
    // FIXME: Maybe fix this so that non-assembled ones have their own name? But currently based upon that
    //        only the non-assembled ones don't have a method between "snc" and "trimmed"
    publishDir "${params.output_dir}/vcf/intermediate/normal", mode: 'copy', overwrite: true, pattern: "*snc*.vcf"
    input:
        tuple val(sample_id), path(split_bam), path(bai)
        path ref_fai
        path ref_dict
        path ref
        each chr_interval

    output:
        tuple val(sample_id), path("*.vcf"), emit: vcfs
        tuple val(sample_id), path("*.tar.gz"), emit: f1r2
        tuple val(sample_id), path("*.stats"), emit: stats

    script:
    def name = "haplotype_${split_bam.simpleName}."
    def interval_args = ""
    for (chr in chr_interval) {
        interval_args += " -L ${chr}"
        name += "${chr}_"
    }

    println "Processing mutect2 in interval ${interval_args}"

    """
    echo "Working on ${split_bam}"
    gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}' Mutect2 \
        --pair-hmm-implementation FASTEST_AVAILABLE \
        --native-pair-hmm-threads ${task.cpus} \
        --smith-waterman FASTEST_AVAILABLE \
        -R ${ref} \
        -I ${split_bam} \
        --f1r2-tar-gz f1r2.tar.gz \
        ${interval_args} \
        -O ${name}.vcf \
    # touch haplotype_${name}.vcf
    # ls -lah
    """

    stub:
    """
    touch haplotype_${split_bam.simpleName}.vcf
    """
}



process MergeBams {
    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("*.bam")

    script:
    def fname = bams[0].simpleName
    def allBams = ""
    for (bam in bams) {
        allBams += "-I ${bam} " // "${bams} " if using samtools
    }

    """
    gatk GatherBamFiles ${allBams} -O ${fname}.bam
    """
}

process MergeVcfs {
    // publishDir "${params.output_dir}/vcf/intermediate/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    // publishDir "${params.output_dir}/vcf/intermediate/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
    publishDir "${params.output_dir}/vcf/unfiltered/normal", mode: 'copy', overwrite: true, pattern: "*snc*.vcf"
    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("*.vcf")

    script:
    def fname = vcfs[0].simpleName
    def allVCFs = ""
    for (vcf in vcfs) {
        allVCFs += "-I ${vcf} " // "${bams} " if using samtools
    }
    """
    gatk MergeVcfs ${allVCFs} -O ${fname}.vcf
    """
}

workflow VARIANT_CALLING {
    take:
        sorted_index_bam // Sample ID + BAM/BAI
        ref
    main:

        // NOTE: Each sample pair gets processed by:
        // 1. Defining chromosomal intervals
        // 2. Pair up the intervals
        // 3. Create SNC_BAM files for each interval, collect
        // 4. Merge collected intervals

        def chromosomes = (1..21) + ['X', 'Y']
        def groupedPairs = [] // i.e [[1,2], [3,4], [5,6]]
        for (int i = 0; i < chromosomes.size(); i += 2) {
            def pair = chromosomes.subList(i, Math.min(i + 2, chromosomes.size()))
            groupedPairs.add(pair)
        }
        // println chromosomes
        // println groupedPairs
        // This would launch 8 (Processes|Groups) * 2 (Chromosomes at a time) * 6 (Cores per process) ~=144 cores
        REF_AUXILLARY(ref)
        bam_split_n = SplitNCigarReads(sorted_index_bam,
                                       REF_AUXILLARY.out.fai,
                                       REF_AUXILLARY.out.dict,
                                       REF_AUXILLARY.out.ref,
                                       groupedPairs) // or [chromosomes]
            .groupTuple() | MergeBams | SAMTOOLS_SORT | MarkDuplicates | SAMTOOLS_INDEX

        haplotype_vcf = Mutect2(bam_split_n,
                                REF_AUXILLARY.out.fai,
                                REF_AUXILLARY.out.dict,
                                REF_AUXILLARY.out.ref,
                                groupedPairs) .groupTuple() | MergeVcfs
        // haplotype_vcf
    emit:
        haplotype_vcf
}

// TODO: Avoid sending tuple from sorted_index_bam to all the others (since some of them require )
// use named outputs?
