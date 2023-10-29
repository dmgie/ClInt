include { SAMTOOLS_INDEX } from './mapping'
include { SAMTOOLS_SORT } from './mapping'
include { VARIANT_PREPROCESSING } from './variant_preprocessing'

workflow VARIANT_CALLING {
    take:
        sorted_index_bam // Sample ID + BAM/BAI
        ref
    main:
        def group_size = 5 // How many intervals each GATK command should take
        def chromosomes = (1..21) + ['X', 'Y']
        def num_lists = ((chromosomes.size() / group_size) + (chromosomes.size() % group_size > 0 ? 1 : 0)) as int
        groups = Channel.fromList(chromosomes).collate(group_size)

        // Create .dict and .fai files for reference fasta file
        REF_AUXILLARY(ref)

        bam_split_n = SplitNCigarReads(sorted_index_bam,
                                       REF_AUXILLARY.out.fai,
                                       REF_AUXILLARY.out.dict,
                                       ref,
                                       groups)
        .groupTuple(size: num_lists) | MergeBams | SAMTOOLS_SORT | MarkDuplicates | SAMTOOLS_INDEX


        // bam_split_n = BAM_PREPROCESSING(sorted_index_bam, REF_AUXILLARY, [groups, num_lists])

        // NOTE: The bams can either be the same or recalibrated, we have it uncalibrated
        // recalibrated =
        recalibrated = bam_split_n // OR VARIANT_PREPROCESSING(bam_split_n,REF_AUXILLARY)

        // Mutect2(bam_split_n,
        Mutect2(recalibrated,
                REF_AUXILLARY.out.fai,
                REF_AUXILLARY.out.dict,
                ref,
                groups)

        // Collect for each sample ID (i.e paired end read set) the (per-chromosome) scattered
        // vcfs & f1r2, to be merged. So each of these outputs will give 1vcf,1tar,1stats for each sample
        vcfs = Mutect2.out.vcfs.groupTuple(size: num_lists) | MergeVcfs
        f1r2 = Mutect2.out.f1r2.groupTuple(size: num_lists) | MergeOrientationModel
        stats = Mutect2.out.stats.groupTuple(size: num_lists) | MergeMutectStats

        // FILTERING
        // Group all needed files together by sample_id and send to Filtering process
        sample_grouped = vcfs.concat(f1r2,stats).groupTuple(size: 3)
        FilterMutect(sample_grouped,
                    REF_AUXILLARY.out.fai,
                    REF_AUXILLARY.out.dict,
                    ref)


    // haplotype_vcf
    // emit:
    //     FilterMutect.out
}


process REF_AUXILLARY {
    // Creation of .fai/.dict files for GATK
    input:
        path ref_file

    output:
        path "${ref_file}.fai", emit: fai
        path "${ref_file.baseName}.dict", emit: dict


    script:
    """
    samtools faidx ${ref_file}
    gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file.baseName}.dict
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
    label 'forking_heavy'
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
    // println "Processing SplitNCigarReads for ${interval_args} for ${sample_id}"

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
    label 'forking_heavy'
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

    // println "Processing mutect2 in interval ${interval_args}"

    """
    echo "Working on ${split_bam}"
    gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}' Mutect2 \
        --pair-hmm-implementation FASTEST_AVAILABLE \
        --smith-waterman FASTEST_AVAILABLE \
        --native-pair-hmm-threads ${task.cpus} \
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

// TODO: Avoid sending tuple from sorted_index_bam to all the others (since some of them require )
// use named outputs?
