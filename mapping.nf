// All mapping-related workflow

workflow MAPPING {
    take:
        reads
        ref_file
        annotation
    main:
        circ = false;
        mapping_method = params.mapping.toLowerCase()
        if (mapping_method == "hisat2") {
            ref_idx = params.hisat_index != '' ? Channel.fromPath("${params.hisat_index}/*").collect() : HISAT_BUILD(ref_file).collect()
            sorted_index_bams = HISAT2(reads, ref_idx) | SAMTOOLS_SORT | SAMTOOLS_INDEX
        } else if (mapping_method == "star") {
            ref_idx = params.star_index != '' ? Channel.fromPath("${params.star_index}/*").collect() :
                    STAR_BUILD(ref_file, annotation).collect()
            STAR(reads,ref_idx,false)

            sorted_index_bams = SAMTOOLS_INDEX(STAR.out.bam)
        } else {
            println "ERROR: Mapping method not recognised"
        }

    // emit:
    // sorted_index_bams // tuple sample_id, bam and bai
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

// NOTE: Output section adapted from nf-core's STAR module
process STAR {
    label 'mapping'
    publishDir "${params.output_dir}/star_mapping/", mode: 'symlink', pattern: '*.bam'
    publishDir "${params.output_dir}/star_log/", mode: 'copy', pattern: '*.{final,SJ}*'

    input:
        tuple val(sample_id), path(reads)
        path  ref_idx
        val   opt                               // Optional for circRNA

    output:
    tuple val(sample_id), path('*.bam')      , emit: bam
    tuple val(sample_id), path('*.tab')                   , optional:true, emit: tab
    tuple val(sample_id), path('*.SJ.out.tab')            , optional:true, emit: spl_junc_tab
    tuple val(sample_id), path('*.out.junction')          , optional:true, emit: junction
    tuple val(sample_id), path('*.out.sam')               , optional:true, emit: sam


    script:
    def sam_header = "ID:aligned_${sample_id}\tSM:None\tLB:None\tPL:Illumina"
    def read_args = "--readFilesIn ${reads}"
    def two_pass = params.star_two_pass ? "--twopassMode Basic" : ""
    def circ_opts = opt ? "--outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15" : ""


    println "STAR on $reads"
    
    """
    echo "Working on ${reads}"
    STAR --runThreadN ${task.cpus} \
        --genomeDir . \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMattrRGline ${sam_header} \
        --outSAMmapqUnique 60 \
        --limitBAMsortRAM 10000000000 \
        ${circ_opts} ${read_args} ${two_pass}

    # Rename file so that downstream (i.e HaplotypeCaller) publishDir works fine
    mv ${sample_id}_*.bam ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}
