#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Set up some params
params.input_reads = file("${params.input_reads_path}/*.f*q.gz");
params.output_path = file("./workflow_output"); // Default to subdir in current dir
REFERENCE = file(params.reference_file); // Require as file() can't be empty
ANNOTATION = file(params.gff_file); // Require as file() can't be empty

// FIXME: Replace "ref_genome" with variables instead
// FIXME: GATK interacts weirdly with .fna files, as in it looks for
//        fname.dict instead of fname.fna.dict, but it works with .fasta; but for ".fai" it looks for fname.fna.fai
//        https://gatk.broadinstitute.org/hc/en-us/community/posts/9688967985691-Does-GATK4-BaseRecalibrator-not-understand-relative-paths-
//        Possible workaround:
//        1. Rename genome files to have .fasta extension when creating/reading (since this initial symlink is propogated)
//        2. Check for the extensions during the Help_file creation, and then adjust the file names accordingly ({ref.baseName}.dict or {ref}.dict)
//        3. TO CHECK: Maybe this happens with fasta as well --> look into it
//        4. Use samtools dict to create the dict file
// FIXME: Simplify ref_fai, ref_dict, and ref_genome to just be ref_files
// FIXME: Tinity-GG output file should have some sort of identifier in the name,
//        since it will be the same for all samples. Use sample name. This causes problems in VCF creation / publishDir

def docker_current_dir() {
    // Simple function to return the docker command for the current process (as the PWD is different for each process)

    // This is a bit of a hack, since the docker container would need to get the reference file from the host,
    // but inputs into processes are symlinked to their original location, so the docker container can't see them

    // That is the reason why we are mounting REF_FULLPATH into the container, so that it can see the reference file
    // from the symlink that is given as input
    REF_FULLPATH = "realpath ${params.reference_file}".execute().text.trim()
    // "docker run -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
    "docker run -v $PWD:$PWD -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
}

def CHECKPARAMS() {
    // Simple function to check if all parameters are set

    println "Checking parameters..."
    if (params.input_reads_path == '') {
        error "ERROR: input_reads_path is not set"
    } else if (params.output_path == '') {
        error "ERROR: output_path is not set"
    } else if (params.reference_file == '') {
        error "ERROR: reference_file is not set"
    } else if (params.gff_file == '') {
        error "ERROR: gff_file is not set"
    } else {
        println "All parameters are set"
        println "   User ID:\t ${params.USER_ID}"
        println "   Group ID:\t ${params.GROUP_ID}"
        // println "   Reference file:\t ${REF_FULLPATH}"
        // println "   Docker command:\t ${DOCKER_COMMAND}"
    }
}

// --------------------Program-specific Processes-------------------------

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
        stdout emit: temp

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

    // """
    // echo "Using hisat2 to align reads ${reads}"
    // hisat2 -p $NUM_THREADS --rg-id ${reads} --rg SM:None --rg LB:None --rg PL:Illumina -x ${index_name} -U ${reads} -S ${aligned_fname}.sam
    // # samtools sort -@ $NUM_THREADS -o ${aligned_fname}.bam ${aligned_fname}.sam
    // # touch aligned_${reads}.bam
    // # ls -lah
    // """
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
        path "Trinity.fasta"

    script:
    def NUM_THREADS = 4
    """
    echo "Working on ${reads}"
    echo "Trinity --seqType fq --max_memory 10G --left ${reads} --right ${reads} --CPU $NUM_THREADS --output Trinity.fasta"
    touch Trinity.fasta
    """
}

process TRINITY_GUIDED {
    maxForks 6
    input:
        path sorted_aligned_bam

    output:
        path "trinity/Trinity-GG*.fasta"


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
    #samtools sort -@ $NUM_THREADS -o sorted_${sorted_aligned_bam} ${sorted_aligned_bam}
    #touch Trinity_GG_${sorted_aligned_bam}.fasta
    """
        
}

process RNASpades {
    maxForks 5
    input:
        path reads

    output:
        path "transcripts.fasta"

    """
    touch transcripts.fasta
    """
}

process SplitNCigarReads {
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
    $docker broadinstitute/gatk gatk --java-options '-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
    --pair-hmm-implementation FASTEST_AVAILABLE \
    --smith-waterman FASTEST_AVAILABLE \
    -R \$PWD/${ref} -I \$PWD/${split_bam} -O \$PWD/haplotype_${split_bam}.vcf
    #touch haplotype_${split_bam.simpleName}.vcf
    ls -lah
    """
}

// -----------------------------------------------------------------------

process QUALITYCONTROL {
    maxForks 5
    // if (${PUBLISH_DIRECTORIES} == 1) {
    //     publishDir "${params.output_path}/QC" // mode: 'copy',
    // }
    
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
        path "${ref_file}.fai" // val used to allow for multiple consumption
        path "${ref_file.baseName}.dict" // val used to allow for multiple consumption
        // path "${ref_file.baseName}.dict" // val used to allow for multiple consumption
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

            aligned_bams = HISAT2(ref_idx, reads)[0]
            sorted_bams = SAMTOOLS_SORT(aligned_bams)

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


workflow ASSEMBLY {
    take:
        reads
        bam_files
    main:
        // Check each element of the assembly array, and run the appropriate assembly processes
        // (e.g. Trinity, RNASpades, etc.)

        // NOTE: For some reason placing Mapping outside oeach of those things, 
        params.assembly.each { method ->
            
            method = method.toLowerCase()
            if ("${method}" == "trinity") {
                trinity_mode = params.trinity_type.toLowerCase()
                if ("${trinity_mode}" == "denovo") {
                    transcripts_fasta = TRINITY_DENOVO(reads) // Returns the Trinity.fasta files
                    // TRINITY_DENOVO.out.view()
                } else if ("${trinity_mode}" == "guided") {
                    transcripts_fasta = TRINITY_GUIDED(bam_files) // Returns the Trinity.fasta files
                    // TRINITY_GUIDED.out.view()
                } else {
                    println "ERROR: Assembly method not recognised"
                }

                // Realign transcripts (fasta) to reference again
                MAPPING(REFERENCE, transcripts_fasta) // Will output bams

            } else if ("${method}" == "rnaspades") {
                // transcripts_fasta = RNASpades(reads)
                // RNASpades.out.view()
                MAPPING(REFERENCE, transcripts_fasta) // Will output bams
            } else {
                println "ERROR: Assembly method \"${method}\" not recognised"
            }
        }
        bams = MAPPING.out
    emit:
        bams
}

workflow GATK {
    take:
        bam
        ref_fai
        ref_dict
        ref
    main:
        split_bam = SplitNCigarReads(bam, ref_fai, ref_dict, ref)
        SplitNCigarReads.out.view()
        haplotype_vcf = HaplotypeCaller(split_bam, ref_fai, ref_dict, ref)
    emit:
        // split_bam
        haplotype_vcf
}


workflow {
    CHECKPARAMS()
    READS = QUALITYCONTROL(Channel.fromPath(params.input_reads))
    


    REFERENCE_HELP_FILES(REFERENCE)
    (ref_fai, ref_dict) = REFERENCE_HELP_FILES.out

    MAPPING(REFERENCE, READS) // Will output bams
    ASSEMBLY(READS, MAPPING.out) // Will output fasta

    // FIXME: Give GATK both the MAPPING.out and ASSEMBLY.out (both bam files)
    GATK(MAPPING.out.concat(ASSEMBLY.out), ref_fai, ref_dict, REFERENCE)
    GATK.out.view()
}

/*
    * Notes on nextflow:
        * If you wanna see stdout, its "PROCESS.out.view()"
        * If you want to define a variable within a process, use "def" after the input/output section
        * the 'val' keyword can be used an ulimited number of times, in comparison to channles which can only be used once 
            (and are consumed)
            Actually, it seems that DSL2 allows for re-use of channels
        * For passing i.e ".ht2" files all at once, we can do a collect() on the output of the process which produces these index files

*/
