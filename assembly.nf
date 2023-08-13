nextflow.enable.dsl=2

include { MAPPING as REMAPPING } from './mapping.nf'

process TRINITY_DENOVO {
    maxForks 6
    input:
        path reads
    output:
        path "Trinity*.fasta"

    script:
    """
    echo "Working on ${reads}"
    Trinity \
    --seqType fq \
    --max_memory 10G \
    --single ${reads} \
    --CPU ${task.cpus} \
    --output \$PWD/trinity/

    #--left ${reads} --right ${reads} # If paired-ends
    mv trinity/Trinity.fasta TrinityDeNovo_${sorted_aligned_bam}.fasta
    """

    stub:
    """
    touch TrinityDeNovo_${sorted_aligned_bam}.fasta
    """
}

process TRINITY_GUIDED {
    maxForks 6
    input:
        path sorted_aligned_bam

    output:
        path "Trinity-GG_${sorted_aligned_bam}.fasta"
        // path "*Trinity-GG*.fasta"



    script:
    def max_intron = 10000
    """ 
    ls -lah
    Trinity \
    --genome_guided_bam \$PWD/$sorted_aligned_bam \
    --genome_guided_max_intron ${max_intron} \
    --max_memory 40G \
    --CPU ${task.cpus} \
    --output \$PWD/trinity/

    mv trinity/Trinity-GG.fasta Trinity-GG_${sorted_aligned_bam}.fasta
    """
 
    stub:
    """
    touch Trinity-GG_${sorted_aligned_bam}.fasta
    """
}

process RNASpades {
    maxForks 5
    input:
        path reads

    output:
        path "rnaspades_${reads.baseName}.fasta"

    script:
    """
    spades.py --rna -t ${task.cpus} -s ${reads} -o spades_out
    mv spades_out/transcripts.fasta rnaspades_${reads.baseName}.fasta
    """

    stub:
    """
    touch rnaspades_${reads.baseName}.fasta
    """
}


workflow ASSEMBLY {
    take:
        reads
        bam_files
        ref
    main:
        // Check each element of the assembly array, and run the appropriate assembly processes
        // (e.g. Trinity, RNASpades, etc.)

        // TODO: Do we need to do BAM sorting here? Or can we just use BAM files straight?

        // FIXME: Find a way to avoid needing to initialise empty channels & concat/mix them later
        //        as well as maybe a way to simplify the code

        // Initialising channels to empty files
        transcripts_fasta = Channel.empty()
        denovo = Channel.empty()
        guided = Channel.empty()
        spades = Channel.empty()

        params.assembly.each { method ->
            method = method.toLowerCase()

            if ("${method}" == "trinity") {
                trinity_mode = params.trinity_type.toLowerCase()
                if ("${trinity_mode}" == "denovo") {
                    denovo = TRINITY_DENOVO(reads)
                } else if ("${trinity_mode}" == "guided") {
                    guided = TRINITY_GUIDED(bam_files)
                } else {
                    error "ERROR: Assembly method not recognised"
                }
            } else if ("${method}" == "rnaspades") {
                spades = RNASpades(reads)
            } else {
                error "ERROR: Assembly method \"${method}\" not recognised"
            }
        }
        // Only accept fasta files
        REMAPPING(ref, transcripts_fasta.mix(guided,spades,denovo)) 
        bams = REMAPPING.out
    emit:
        bams
}

