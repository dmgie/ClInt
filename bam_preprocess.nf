include {SAMTOOLS_SORT, SAMTOOLS_INDEX} from './mapping'

// MergeBams, MarkDuplicates, SplitNCigar


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


workflow BAM_PREPROCESS {
    take:
    bam_and_idx
    ref_files
    grouping

    main:
    def (groups, listSize) = grouping
        bam_split_n = SplitNCigarReads(bam_and_idx,
                                    ref_file.out.fai,
                                    ref_file.out.dict,
                                    ref_file.out.ref,
                                    grouping)
            .groupTuple(size: listSize) | MergeBams | SAMTOOLS_SORT | MarkDuplicates | SAMTOOLS_INDEX

    emit:



}
