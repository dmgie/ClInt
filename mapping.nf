process SAMTOOLS_SORT {
    maxForks 5
    input:
        path sam_file

    output:
        path "*.bam"

    script:

    """
    echo "Sorting ${sam_file}"
    samtools sort -@ ${task.cpus} -o ${sam_file.simpleName}.bam ${sam_file}
    """

    stub:
    """
    touch ${sam_file.simpleName}.bam
    """
}

process HISAT_BUILD {
    maxForks 5
    input:
        path ref_file
    output:
        path "ref_idx*.ht2"

    script:
    def index_name = "ref_idx"

    """
    echo "Running hisat2-build"
    hisat2-build -p ${task.cpus}  ${ref_file} ${index_name}
    """

    stub:
    """
    touch ref_idx.1.ht2
    touch ref_idx.2.ht2
    touch ref_idx.3.ht2
    """
}

process HISAT2 {
    maxForks 5
    publishDir "${params.output_dir}/hisat2_summaries", mode: 'copy', pattern: '*.txt'
    input:
        tuple val(sample_id), path(reads)
        path ref_idx
    output:
        path "*.sam"

    script:
    def (read1,read2) = [reads[0], reads[1]]
    def index_name = "ref_idx"
    def read_args = params.paired ? "-1 ${read1} -2 ${read2}" : "-U ${read1}"
    def extension_args = read1.extension == "fasta" ? "-f" : ""

    """
    hisat2 ${extension_args} -p ${task.cpus} \
    --new-summary --summary-file ${sample_id}_summary.txt \
    --rg-id ${reads} --rg SM:None --rg LB:None --rg PL:Illumina \
    -x ${index_name} \
    ${read_args} \
    -S ${sample_id}.sam
    """

    stub:
    """
    touch ${sample_id}.sam
    touch ${sample_id}_summary.txt
    """
}

process STAR_BUILD {
    maxForks 5
    input:
        path ref_file
        path annotation
    output:
        path "*"

    script:
    def READ_LENGTH = 100
    def feature = "gene"
    def extension = annotation.extension // Remove this when changing to *_args
    def extension_args = annotation.extension == "gtf" ? "" :
        "--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon $feature"

    // TODO: Modify it to be a single command but with ext_args
    if (extension == 'gtf') {
        """
            echo "GTF file detected"
            STAR --runThreadN ${task.cpus} \
            --runMode genomeGenerate \
            --genomeDir . \
            --genomeFastaFiles $ref_file \
            --sjdbGTFfile $annotation \
            --sjdbOverhang \$(($READ_LENGTH - 1))
        """
    } else {
        """
            echo "GFF file detected"
            STAR --runThreadN ${task.cpus} \
            --runMode genomeGenerate \
            --genomeDir . \
            --genomeFastaFiles $ref_file \
            --sjdbGTFfile $annotation \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbGTFfeatureExon $feature \
            --sjdbOverhang \$(($READ_LENGTH - 1))
        """
    }

    stub:
    """
    touch annotations.txt
    """
}

process STAR {
    publishDir "${params.output_dir}/star_summaries", mode: 'copy', pattern: '*.final.out'
    publishDir "${params.output_dir}/bams", mode: 'copy', pattern: '*.bam'
    maxForks 5

    input:
        tuple val(sample_id), path(reads)
        path ref_idx

    output:
        path "*.bam"

    script:
    // No strand-specific options needed here
    def (read1,read2) = [reads[0], reads[1]]
    def SAM_HEADER = "ID:aligned_${reads}\tSM:None\tLB:None\tPL:Illumina"
    def arguments = params.paired ? "--readFilesIn ${read1} ${read2}" :
                                    "--readFilesIn ${read1}"

    """
    echo "Working on ${reads}"
    STAR --runThreadN ${task.cpus} \
    --genomeDir . \
    --readFilesIn ${reads} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${aligned_fname}_ \
    --outSAMattrRGline $SAM_HEADER \
    --limitBAMsortRAM 10000000000

    # Rename file so that downstream publishDir works fine
    mv ${aligned_fname}_*.bam ${aligned_fname}.bam
    """

    stub:
    """
    touch ${reads.simpleName}.bam
    """
}

workflow MAPPING {
    take: 
        reads
        ref_file
    main:
        mapping_method = params.mapping.toLowerCase()
        if (mapping_method == "hisat2") {
            ref_idx = HISAT_BUILD(ref_file).collect() // Collect all the output files from HISAT_BUILD to be used as input for HISAT2 
                                                        // (so they don't get individually consumed)
            sorted_bams = HISAT2(ref_idx, reads) | SAMTOOLS_SORT
        } else if (mapping_method == "star") {
            ref_idx = STAR_BUILD(ref_file, Channel.fromPath(params.gff_file)).collect()
            sorted_bams = STAR(ref_idx, reads) // The command itself aligns the bams
        } else {
            println "ERROR: Mapping method not recognised"
        }
    emit:
        sorted_bams // output STAR/Hisat2 bam files (sorted)
}

