include { SAMTOOLS_INDEX } from './mapping';

workflow VARIANT_PREPROCESSING {
    take:
        bam_bai
        ref_files
    main:
    // TODO: Check if for each of the VCF in params.known_sites theres a tbi file
    // otherwise create it
    // def known_sites = Channel.fromFilePairs("${params.known_sites}")
        vcfs = Channel.fromPath(params.known_sites)
        list = IndexVCF(vcfs).collect()

        CalculateRecalibration(bam_bai,
                            list,
                            ref_files.out.fai,
                            ref_files.out.dict,
                            ref_files.out.ref)

        recalibrated = ApplyRecalibration(CalculateRecalibration.out,
                        ref_files.out.fai,
                        ref_files.out.dict,
                        ref_files.out.ref) | SAMTOOLS_INDEX
    emit:
        recalibrated
}


process IndexVCF {
    input:
    path vcf

    output:
    tuple path(vcf), path("*.tbi")

    script:
    println "Indexing ${vcf}"
    """
    gatk IndexFeatureFile -I ${vcf};
    """
}

process CalculateRecalibration {
    input:
        tuple val(sample_id), path(bam), path(bai)
        path vcfs
        path ref_fai
        path ref_dict
        path ref

    output:
        tuple val(sample_id), path(bam), path(bai), path("*.table")

    script:
    println "Using ${vcfs} for ${sample_id}";
    vcfgz = vcfs.collect().findAll { !it.name.endsWith(".tbi") };
    def sites_arg = "" // Give all --known-sites files as parameters
    for (vcf in vcfgz) {
        println "Adding ${vcf} to site arg"
        sites_arg += " --known-sites ${vcf}" 
    }
    """
    gatk BaseRecalibrator -I ${bam} \
        -R ${ref} \
        ${sites_arg} \
        -O ${sample_id}_recal_table.table
    """
}

process ApplyRecalibration {
    input:
        tuple val(sample_id), path(bam), path(bai), path(table)
        path ref_fai
        path ref_dict
        path ref

    output:
        tuple val(sample_id), path("*_recal.bam")

    script:
    """
    gatk ApplyBQSR \
        -R ${ref} \
        -I ${bam} \
        --bqsr-recal-file ${table} \
        -O ${sample_id}_recal.bam
    """
}
