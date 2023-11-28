


workflow MAPPING {
    take:
        reads
        ref_file
        annotation
    main:
        mapping_method = params.mapping.toLowerCase()
        if (mapping_method == "hisat2") {
            ref_idx = params.genome_index != '' ? Channel.fromPath("${params.genome_index}/*").collect() : HISAT_BUILD(ref_file).collect()
            sorted_index_bams = HISAT2(reads, ref_idx) | SAMTOOLS_SORT | SAMTOOLS_INDEX
        } else if (mapping_method == "star") {
            ref_idx = params.genome_index != '' ? Channel.fromPath("${params.genome_index}/*").collect() : STAR_BUILD(ref_file, annotation).collect()
            // sorted_index_bams = STAR_CIRCRNA(reads,ref_idx)
            // sorted_index_bams = SAMTOOLS_INDEX(STAR_CIRCRNA.out.bam)
            STAR_CIRCRNA(reads,ref_idx)
            STAR_CIRCRNA_CIRCTOOLS(reads,ref_idx)
            STAR(reads,ref_idx)
            sorted_index_bams = SAMTOOLS_INDEX(STAR.out.bam)
        } else {
            println "ERROR: Mapping method not recognised"
        }

    emit:
    sorted_index_bams // tuple sample_id, bam and bai
}


process SAMTOOLS_SORT {
    label 'mapping'
    input:
        tuple val(sample_id), path(xam_file) //bam or sam

    output:
        tuple val(sample_id), path("${xam_file.simpleName}.bam")

    script:

    """
    echo "Sorting ${xam_file}"
    mv ${xam_file} temp.ext # Avoid bam->bam name collision
    samtools sort -@ ${task.cpus} -o ${xam_file.simpleName}.bam temp.ext
    #samtools sort -@ ${task.cpus} -o ${xam_file.simpleName}.bam ${xam_file}
    """

    stub:
    """
    touch ${xam_file.simpleName}.bam
    """
}


process SAMTOOLS_INDEX {
    publishDir "${params.output_dir}/bam_index/", mode: 'symlink', pattern: "*.bai"
    label 'mapping'
    input:
        tuple val(sample_id), path(bam_file)
    output:
        tuple val(sample_id), path(bam_file), path("*.bai")

    script:
    """
    echo "Sorting ${bam_file} to ${bam_file.baseName}"
    samtools index -@ ${task.cpus} ${bam_file} ${bam_file.baseName}.bai
    """

    stub:
    """
    touch ${bam_file.simpleName}.bai
    """
}

process HISAT_BUILD {
    label 'mapping'
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
    label 'mapping'
    publishDir "${params.output_dir}/hisat2_summaries", mode: 'copy', pattern: '*.txt'
    input:
        tuple val(sample_id), path(reads)
        path ref_idx
    output:
        tuple val(sample_id), path("*.sam")

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
    // NOTE: REMOVE v37 ONCE PATIENT HAS RUN THROUGH
    publishDir "${params.output_dir}/star_index", mode: 'copy', pattern: '*'
    label 'mapping'
    input:
        path ref_file
        path annotation
    output:
        path "*"

    script:
    def READ_LENGTH = 100
    def feature = "exon"
    def extension = annotation.extension // Remove this when changing to *_args
    def extension_args = annotation.extension == "gtf" ? "" :
        "--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon $feature"
    """
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir . \
        --genomeFastaFiles $ref_file \
        --sjdbGTFfile $annotation \
        --sjdbOverhang \$(($READ_LENGTH - 1)) ${extension_args}
    """

    stub:
    """
    touch annotations.txt
    """
}

process STAR {
    label 'mapping'
    publishDir "${params.output_dir}/star_mapping/", mode: 'symlink', pattern: '*.bam'
    publishDir "${params.output_dir}/star_log/", mode: 'copy', pattern: '*.{final,SJ}*'

    input:
        tuple val(sample_id), path(reads)
        path ref_idx

    output:
    tuple val(sample_id),  path("*.bam"), emit: bam
    path "*{final,SJ}*", emit: logs

    script:
    // No strand-specific options needed here
    def (read1,read2) = [reads[0], reads[1]]
    def SAM_HEADER = "ID:aligned_${sample_id}\tSM:None\tLB:None\tPL:Illumina"
    def read_arguments = params.paired ? "--readFilesIn ${read1} ${read2}" :
                                    "--readFilesIn ${read1}"
    def two_pass = params.star_two_pass ? "--twopassMode Basic" : ""
    """
    echo "Working on ${reads}"
    STAR --runThreadN ${task.cpus} \
        --genomeDir . \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMattrRGline $SAM_HEADER \
        --outSAMmapqUnique 60 \
        --limitBAMsortRAM 10000000000 \
        ${read_arguments} ${two_pass}

    # Rename file so that downstream (i.e HaplotypeCaller) publishDir works fine
    mv ${sample_id}_*.bam ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}


process STAR_CIRCRNA {
    label 'mapping'
    publishDir "${params.output_dir}/star_circ_mapping/", mode: 'symlink', pattern: '*.bam'
    publishDir "${params.output_dir}/star_circ_log/", mode: 'copy', pattern: '*.{final,SJ}*'
    publishDir "${params.output_dir}/star_circ_chimeric/", mode: 'copy', pattern: '*[Cc]himeric.out.junction*'
    input:
    tuple val(sample_id), path(reads)
    path ref_idx

    output:
    tuple val(sample_id),  path("*.bam"), emit: bam
    path "*{final,SJ}*", emit: logs
    path "*[Cc]himeric*", emit: chimera

    // Options adopted from `STAR-Fusion` documentation
    script:
    def (read1,read2) = [reads[0], reads[1]]
    def SAM_HEADER = "ID:aligned_${sample_id}\tSM:None\tLB:None\tPL:Illumina"
    def read_arguments = params.paired ? "--readFilesIn ${read1} ${read2}" :
        "--readFilesIn ${read1}"
    def two_pass = params.star_two_pass ? "--twopassMode Basic" : ""
    """
    echo "Working on ${reads}"
    STAR --runThreadN ${task.cpus} \
        --genomeDir . \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMattrRGline $SAM_HEADER \
        --outSAMmapqUnique 60 \
        --limitBAMsortRAM 10000000000 \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        ${two_pass} ${read_arguments}

    # Rename file so that downstream (i.e HaplotypeCaller) publishDir works fine
    mv ${sample_id}_*.bam ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}

process STAR_CIRCRNA_CIRCTOOLS {
    label 'mapping'
    publishDir "${params.output_dir}/star_circtools/${sample_id}", mode: 'symlink', pattern: 'main/*'
    publishDir "${params.output_dir}/star_circtools/${sample_id}/mate1", mode: 'symlink', pattern: 'mate1/*'
    publishDir "${params.output_dir}/star_circtools/${sample_id}/mate2", mode: 'symlink', pattern: 'mate2/*'
    input:
    tuple val(sample_id), path(reads)
    path annotation
    path ref_idx

    output:
    path "*"

    // This is specific for the `circtools` tool.
    // To be more sensitive, it runs STAR 3 times (for paired-end data)
    // Once with both mates, and then separately for each one
    //
    script:
    def (read1,read2) = [reads[0], reads[1]]
    def SAM_HEADER = "ID:aligned_${sample_id}\tSM:None\tLB:None\tPL:Illumina"
    def read_arguments = params.paired ? "--readFilesIn ${read1} ${read2}" :
        "--readFilesIn ${read1}"
    """
    # Requires almost deprecated SeparateSAMold function to work
    # Run once with pair, then separately
    STAR --runThreadN ${task.cpus} \
        --genomeDir . \
        --genomeLoad NoSharedMemory \
        --readFilesCommand zcat  \
        ${read_arguments}  \
        --outFileNamePrefix main/${sample_id}_ \
        --sjdbGTFfile ${annotation} \
        --outReadsUnmapped Fastx \
        --outSAMattributes NH HI AS nM NM MD jM jI XS \
        --outSJfilterOverhangMin 15   15   15   15 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNminOverLread 0.7 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.05 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 10 \
        --alignSoftClipAtReferenceEnds No \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15 \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --chimOutType Junctions SeparateSAMold

    STAR --runThreadN ${task.cpus} \
        --genomeDir . \
        --genomeLoad NoSharedMemory \
        --readFilesCommand zcat  \
        --readFilesIn ${read1}  \
        --outFileNamePrefix mate1/${read1.baseName}_ \
        --sjdbGTFfile ${annotation} \
        --outReadsUnmapped Fastx \
        --outSAMattributes NH HI AS nM NM MD jM jI XS \
        --outSJfilterOverhangMin 15   15   15   15 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNminOverLread 0.7 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.05 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 10 \
        --alignSoftClipAtReferenceEnds No \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15 \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --chimOutType Junctions SeparateSAMold

    STAR --runThreadN ${task.cpus} \
        --genomeDir . \
        --genomeLoad NoSharedMemory \
        --readFilesCommand zcat  \
        --readFilesIn ${read1}  \
        --outFileNamePrefix mate2/${read2.baseName}_ \
        --sjdbGTFfile ${annotation} \
        --outReadsUnmapped Fastx \
        --outSAMattributes NH HI AS nM NM MD jM jI XS \
        --outSJfilterOverhangMin 15   15   15   15 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNminOverLread 0.7 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.05 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 10 \
        --alignSoftClipAtReferenceEnds No \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15 \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --chimOutType Junctions SeparateSAMold

    # Rename file so that downstream (i.e HaplotypeCaller) publishDir works fine
    mv ${sample_id}_*.bam ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}

// process BamQC {
//     input:
//     output:
//     script:
//     """
//     """
// }
