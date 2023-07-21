process QUALITYCONTROL {
    maxForks 5
    
    input:
        path read

    output:
        path "trimmed_*"
        // stdout emit: verb


    script:
    """
    echo "Working on ${read}"
    fastp -i ${read} -o trimmed_${read} -j /dev/null -h /dev/null
    """
}

process REFERENCE_HELP_FILES {
    // This process mainly relates to the necessary side-files
    // such as the .fai and .dict files for the reference genome
    // as these are required in some of the processes

    input: 
        path ref_file

    output:
        path "${ref_file}.fai" 
        path "${ref_file.baseName}.dict" 
        // stdout emit: verbo


    def docker = docker_current_dir()
    script:
    """
    echo "Running samtools faidx and docker"
    samtools faidx ${ref_file}
    ${docker} broadinstitute/gatk gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file.baseName}.dict
    #${docker} broadinstitute/gatk gatk CreateSequenceDictionary -R \$PWD/${ref_file} -O \$PWD/${ref_file}.dict
    ls -lah
    """
}


process SAMTOOLS_SORT {
    maxForks 5
    input:
        path sam_file

    output:
        path "*.bam"

    script:
    def NUM_THREADS = 4
    """
    echo "Sorting ${sam_file}"
    samtools sort -@ $NUM_THREADS -o ${sam_file.simpleName}.bam ${sam_file}
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
}

process HISAT2 {
    maxForks 5
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
        hisat2 -f -p $NUM_THREADS --rg-id ${reads} --rg SM:None --rg LB:None --rg PL:Illumina -x ${index_name} -U ${reads} -S ${aligned_fname}.sam
        """
    } else {
        """
        echo "Fastq file detected"
        hisat2 -p $NUM_THREADS --rg-id ${reads} --rg SM:None --rg LB:None --rg PL:Illumina -x ${index_name} -U ${reads} -S ${aligned_fname}.sam
        """
    }
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
            touch one.txt
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
            touch one.txt
        """
    }
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
    touch aligned_${reads.simpleName}.bam
    """
}


process TRINITY_DENOVO {
    maxForks 6
    input:
        path reads
    output:
        path "Trinity*.fasta"

    script:
    def NUM_THREADS = 4
    """
    echo "Working on ${reads}"
    $docker trinityrnaseq/trinityrnaseq Trinity \
    --seqType fq \
    --max_memory 10G \
    --single ${reads} \
    --CPU $NUM_THREADS \
    --output \$PWD/trinity/

    #--left ${reads} --right ${reads} # If paired-ends
    mv trinity/Trinity.fasta TrinityDeNovo_${sorted_aligned_bam}.fasta
    """
}

process TRINITY_GUIDED {
    maxForks 6
    input:
        path sorted_aligned_bam

    output:
        path "*Trinity-GG*.fasta"
        // path "Trinity-GG_${sorted_aligned_bam}.fasta"


    def docker = docker_current_dir()

    script:
    def NUM_THREADS = 4
    def max_intron = 10000
    """ 
    ls -lah
    $docker trinityrnaseq/trinityrnaseq Trinity \
    --genome_guided_bam \$PWD/$sorted_aligned_bam \
    --genome_guided_max_intron ${max_intron} \
    --max_memory 40G \
    --CPU $NUM_THREADS \
    --output \$PWD/trinity/

    mv trinity/Trinity-GG.fasta Trinity-GG_${sorted_aligned_bam}.fasta
    """
        
}

process RNASpades {
    maxForks 5
    input:
        path reads

    output:
        path "rnaspades_${reads.baseName}.fasta"

    """
    spades.py --rna -ss-fr ${reads} -o spades_out
    mv spades_out/transcripts.fasta rnaspades_${reads.baseName}.fasta
    """
}

process MarkDuplicates {
    // NOTE: We an add --REMOVE_DUPLICATES=true to remove duplicates from the final BAM file
    //       intead of just switching the flag for that read
    maxForks 5
    input:
        path aligned_bam

    output:
        path "dedup_${aligned_bam}"

    def docker = docker_current_dir()

    script:
    """
    echo "Working on ${aligned_bam}"
    ${docker} broadinstitute/gatk gatk MarkDuplicates -I \$PWD/${aligned_bam} -O \$PWD/dedup_${aligned_bam} -M \$PWD/dedup_${aligned_bam}.metrics
    """
}

process SplitNCigarReads {
    // TODO: Change aligned_bam to dedup_bam, and change the output to dedup_${aligned_bam}, because we wanna remove duplicates first
    maxForks 8
    input:
        path aligned_bam
        path ref_fai
        path ref_dict
        path ref

    output:
        path "snc_${aligned_bam}"
        // stdout emit: temp
    
    def docker = docker_current_dir()

    script:
    def NUM_THREADS = 4
    """
    echo "Working on ${aligned_bam}"

    ${docker} broadinstitute/gatk gatk SplitNCigarReads -R \$PWD/${ref} -I \$PWD/${aligned_bam} -O \$PWD/snc_${aligned_bam}

    # ${docker} broadinstitute/gatk bash -c "cd \$PWD/; ls -lah ../../"
    # ${docker} broadinstitute/gatk gatk SplitNCigarReads -R \$PWD/${ref} -I \$PWD/${aligned_bam} -O \$PWD/snc_${aligned_bam}
    # echo "gatk SplitNCigarReads -R ${ref} -I ${aligned_bam} -O split_${aligned_bam}"
    # touch snc_${aligned_bam}
    # ls -lah
    """
}

process HaplotypeCaller {
    maxForks 8
    publishDir "${params.output_dir}/haplotype_vcf/", mode: 'copy', overwrite: true
    input:
        path split_bam
        path ref_fai
        path ref_dict
        path ref

    output:
        path "haplotype_*.vcf"

    def docker = docker_current_dir()

    script:
    def NUM_THREADS = 4
    """
    echo "Working on ${split_bam}"
    samtools index ${split_bam}
    $docker broadinstitute/gatk gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
    --pair-hmm-implementation FASTEST_AVAILABLE \
    --smith-waterman FASTEST_AVAILABLE \
    -R \$PWD/${ref} -I \$PWD/${split_bam} -O \$PWD/haplotype_${split_bam}.vcf
    #touch haplotype_${split_bam.simpleName}.vcf
    ls -lah
    """
}
