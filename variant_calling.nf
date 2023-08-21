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
    label 'variant_calling'
    publishDir "${params.output_dir}/bams/deduped/", mode: 'symlink', overwrite: true, pattern: "*.bam"
    publishDir "${params.output_dir}/bams/deduped/", mode: 'symlink', overwrite: true

    input:
        tuple val(sample_id), path(aligned_bam)

    output:
        tuple val(sample_id), path("*.bam")

    script:
    // NOTE: We an add --REMOVE_DUPLICATES=true to remove duplicates from the final BAM file
    //       intead of just switching the flag for that read
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
        each chr_interval // i.e [[1,2,..], [3,4,..]], does for each [x,y,..]  set

    output:
        tuple val(sample_id), path("snc_*.bam")

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
    def name = "${split_bam.simpleName}"
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
        --f1r2-tar-gz f1r2_${name}.tar.gz \
        ${interval_args} \
        -O ${name}.vcf \
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
    publishDir "${params.output_dir}/vcf/unfiltered/", mode: 'copy', overwrite: true, pattern: "*.vcf"
    // publishDir "${params.output_dir}/vcf/unfiltered/rna_spades", mode: 'copy', overwrite: true, pattern: "*spades_*.vcf"
    // publishDir "${params.output_dir}/vcf/unfiltered/Trinity-GG", mode: 'copy', overwrite: true, pattern: "*Trinity-GG_*.vcf"
    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("merged_${sample_id}.vcf")

    script:
    def fname = vcfs[0].simpleName
    def allVCFs = ""
    for (vcf in vcfs) {
        allVCFs += "-I ${vcf} " // "${bams} " if using samtools
    }
    """
    gatk MergeVcfs ${allVCFs} -O merged_${sample_id}.vcf
    """
}

workflow VARIANT_CALLING {
    take:
        sorted_index_bam // Sample ID + BAM/BAI
        ref
    main:


        // NOTE: Each sample_id (i.e paired-end pair) gets processed by:
        // 1. Defining chromosomal intervals (sets of size i.e 3)
        // 3. BAM -> Merged Bam & VCF files for each interval, collect
        // 4. Merge collected intervals
        def group_size = 5 // How many intervals each GATK command should take
        def chromosomes = (1..21) + ['X', 'Y']

        // NOTE: This (along with the "size: num_lists") parameter allows nextflow to know
        // how many elements to expect for each sample_id. This allows it know when it can start
        // the next process much faster rather than waiting on the current one for all to finish
        // If we have 23 chromosomes, and pair them in 2's, we have at the end 12 intervals, so we supply "12" to the
        // "size" parameter in "groupTuple"
        def num_lists = ((chromosomes.size() / group_size) + (chromosomes.size() % group_size > 0 ? 1 : 0)) as int
        println num_lists

        groups = Channel.fromList(chromosomes).collate(group_size)
        groups.view()

        // This would launch 8 (Processes|Groups) * 2 (Chromosomes at a time) * 6 (Cores per process) ~=144 cores
        REF_AUXILLARY(ref)
        bam_split_n = SplitNCigarReads(sorted_index_bam,
                                       REF_AUXILLARY.out.fai,
                                       REF_AUXILLARY.out.dict,
                                       REF_AUXILLARY.out.ref,
                                       groups) // or [chromosomes]
        .groupTuple(size: num_lists) | MergeBams | SAMTOOLS_SORT | MarkDuplicates | SAMTOOLS_INDEX


        Mutect2(bam_split_n,
                REF_AUXILLARY.out.fai,
                REF_AUXILLARY.out.dict,
                REF_AUXILLARY.out.ref,
                groups)

    // Collect for each sample ID (i.e paired end read set) the (per-chromosome) scattered
    // vcfs & f1r2, to be merged
        vcfs = Mutect2.out.vcfs.groupTuple()
        f1r2 = Mutect2.out.f1r2.groupTuple()


    // haplotype_vcf
    // emit:
    //     Mutect2.out
}

// TODO: Avoid sending tuple from sorted_index_bam to all the others (since some of them require )
// use named outputs?
