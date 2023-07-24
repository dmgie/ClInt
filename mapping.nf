process SAMTOOLS_SORT {
    maxForks 5
    input:
        path sam_file

    output:
        path "*.bam"

    script:
    def NUM_THREADS = 4

    // NOTE: If we want certain version-specific samtools features, we can do
    // do a version check by i.e running `samtools --version` and parsing the output
    // then comparing it to the version we want. 
    // Easiest way is add them to a list, sort them lexically and then compare the first element
    // to the version we want. If it's >=, then we can use the new features, otherwise we can't.
    """
    echo "Sorting ${sam_file}"
    samtools sort -@ $NUM_THREADS -o ${sam_file.simpleName}.bam ${sam_file}
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
    def NUM_THREADS = 4
    def index_name = "ref_idx"
    """
    echo "Running hisat2-build"
    hisat2-build -p $NUM_THREADS  ${ref_file} ${index_name}
    # touch ref_idx.1.ht2
    # touch ref_idx.2.ht2
    # touch ref_idx.3.ht2
    # ls -lah
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
        path ref_idx
        path reads
    output:
        path "*.sam"

    script:
    def aligned_fname = "${reads.simpleName}"
    def NUM_THREADS = 4
    def index_name = "ref_idx"
    def extension = "${reads.extension}"

    println "[LOG] HISAT2 :: Extension is ${extension}"
    if (extension == 'fasta') {
        """
        echo "Fasta file detected"
        hisat2 -f -p $NUM_THREADS --new-summary --summary-file ${aligned_fname}_summary.txt --rg-id ${reads} --rg SM:None --rg LB:None --rg PL:Illumina -x ${index_name} -U ${reads} -S ${aligned_fname}.sam
        """
    } else {
        """
        echo "Fastq file detected"
        hisat2 -p $NUM_THREADS --new-summary --summary-file ${aligned_fname}_summary.txt --rg-id ${reads} --rg SM:None --rg LB:None --rg PL:Illumina -x ${index_name} -U ${reads} -S ${aligned_fname}.sam
        """
    }

    stub:
    def aligned_fname = "${reads.simpleName}"
    """
    touch ${aligned_fname}.sam
    touch ${aligned_fname}.txt
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
    def NUM_THREADS = 4
    def READ_LENGTH = 100
    def feature = "gene"
    def extension = annotation.extension
    if (extension == 'gtf') {
        """
            echo "GTF file detected"
            echo "STAR --runThreadN $NUM_THREADS \
            --runMode genomeGenerate \
            --genomeDir . \
            --genomeFastaFiles $ref_file \
            --sjdbGTFfile $annotation \
            --sjdbOverhang \$(($READ_LENGTH - 1))"
        """
    } else {
        """
            echo "GFF file detected"
            echo "STAR --runThreadN $NUM_THREADS \
            --runMode genomeGenerate \
            --genomeDir . \
            --genomeFastaFiles $ref_file \
            --sjdbGTFfile $annotation \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbGTFfeatureExon $feature \
            --sjdbOverhang \$(($READ_LENGTH - 1))"
        """
    }

    stub:
    """
    touch one.txt
    """
}

process STAR {
    maxForks 5
    input:
        path ref_idx
        path reads

    output:
        path "aligned_*.bam"

    script:
    def NUM_THREADS = 4
    def SAM_HEADER = "@RG\tID:aligned_${reads}\tSM:None\tLB:None\tPL:Illumina"
    """
    echo "Working on ${reads}"
    echo "STAR --runThreadN $NUM_THREADS \
    --genomeDir ${ref_idx} \
    --readFilesIn ${reads} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix aligned_ \
    --outSAMattrRGline $SAM_HEADER \
    --limitBAMsortRAM 10000000000"
    #touch aligned_${reads.simpleName}.bam
    """

    stub:
    """
    touch aligned_${reads.simpleName}.bam
    """
}

workflow MAPPING {
    take: 
        ref_file
        reads
    main:
        mapping_method = params.mapping.toLowerCase()
        if (mapping_method == "hisat2") {
            ref_idx = HISAT_BUILD(ref_file).collect() // Collect all the output files from HISAT_BUILD to be used as input for HISAT2 
                                                        // (so they don't get individually consumed)
            // HISAT_BUILD.out.view()

            sorted_bams = HISAT2(ref_idx, reads) | SAMTOOLS_SORT
            //sorted_bams = SAMTOOLS_SORT(aligned_bams)

            // HISAT2.out[1].view()
        } else if (mapping_method == "star") {
            ref_idx = STAR_BUILD(ref_file, ANNOTATION).collect()
            // STAR_BUILD.out.view()

            sorted_bams = STAR(ref_idx, reads) // The command itself aligns the bams
            // STAR.out.view()
        } else {
            println "ERROR: Mapping method not recognised"
        }
    emit:
        sorted_bams // output STAR/Hisat2 bam files (sorted)
}
