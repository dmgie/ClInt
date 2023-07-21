include { MAPPING } from '../reuse'

workflow ASSEMBLY {
    take:
        reads
        bam_files
        method
    main:
        // Check each element of the assembly array, and run the appropriate assembly processes
        // (e.g. Trinity, RNASpades, etc.)

        // NOTE: For some reason placing Mapping outside oeach of those things, 

        // TODO: Do we need to do BAM sorting here? Or can we just use BAM files straight?
        if ("$method" == "trinity") {
            trinity_mode = params.trinity_type.toLowerCase()
            if ("$trinity_mode" == "denovo") {
                transcripts_fasta = TRINITY_DENOVO(reads) 
            } else if ("$trinity_mode" == "guided") {
                transcripts_fasta = TRINITY_GUIDED(bam_files) 
            } else {
                println "ERROR: Unknown trinity_mode"
            }
            MAPPING(REFERENCE, transcripts_fasta) // Will output bams
            // Realign transcripts (fasta) to reference again
            bams = MAPPING.out
        } else if ("$method" == "rnaspades") {
            // transcripts_fasta = RNASpades(reads)
            // RNASpades.out.view()
            // Realign transcripts (fasta) to reference again
        } else {
            println "ERROR: Assembly method \"${method}\" not recognised"
        }
    emit:
        bams
}

workflow VARIANT_CALLING {
    // Run SplitNCigarReads and HaplotypeCaller
    take:
        bam
        ref_fai
        ref_dict
        ref
    main:
        split_bam = SplitNCigarReads(bam, ref_fai, ref_dict, ref) | MarkDuplicates
        haplotype_vcf = HaplotypeCaller(split_bam, ref_fai, ref_dict, ref)
    emit:
        haplotype_vcf
}

