// After circRNA do some analyses on miRNA prediction on those sites

workflow MIRNA_PREDICTION{

    take:
    circrna_fasta
    circrna_bed12
    
    main:
    //
    // TARGETSCAN WORKFLOW:
    //
    
    def mature = Channel.fromPath("${params.mirna_mature}").first()
    println "THIS IS THE MATURE MIRNA: $params.mirna_mature"
    TARGETSCAN_DATABASE(mature)
    TARGETSCAN(circrna_fasta, TARGETSCAN_DATABASE.out.mature_txt)

    //
    // MIRANDA WORKFLOW:
    // Split fasta, run miRanda and merge results back together
    split = circrna_fasta.splitFasta(by: 10, file: true, elem: 1)
    MIRANDA(split, mature)
    MIRANDA_MERGE(MIRANDA.out.groupTuple())


    //
    // CONSOLIDATE PREDICTIONS WORKFLOW:
    //

    // consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt).join(circrna_bed12)
    // MIRNA_TARGETS( consolidate_targets )
}

process TARGETSCAN_DATABASE {
    tag "$mature"

    container 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'

    input:
    path(mature)

    output:
    path("mature.txt")  , emit: mature_txt
    
    script:
    """
    targetscan_format.sh $mature
    """
}

process TARGETSCAN {
    tag "${sample}"
    queue "long" // It takes a very long time
    container 'https://depot.galaxyproject.org/singularity/targetscan:7.0--pl5321hdfd78af_0'

    input:
    tuple val(sample), path(fasta)
    path(mature_txt)

    output:
    tuple val(sample), path("${sample}.txt"), emit: txt

    script:
    """
    ##format for targetscan
    cat $fasta | grep ">" | sed 's/>//g' > id
    cat $fasta | grep -v ">" > seq
    paste id seq | awk -v OFS="\t" '{print \$1, "0000", \$2}' > ${sample}_ts.txt
    # run targetscan
    targetscan_70.pl mature.txt ${sample}_ts.txt ${sample}.txt
    """
}

process MIRANDA {
    tag "${sample}"
    
    container 'https://depot.galaxyproject.org/singularity/miranda:3.3a--h779adbc_3'

    input:
    tuple val(sample), path(query)
    path(mirbase)

    output:
    tuple val(sample), path("*.txt"), emit: txt

    script:
    """
    miranda \\
        $mirbase \\
        $query \\
        $args \\
        -out tmp_${sample}.out

    echo "miRNA\tTarget\tScore\tEnergy_KcalMol\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tAln_len\tSubject_Identity\tQuery_Identity" > ${sample}.txt
    grep -A 1 "Scores for this hit:" ${sample}.out | sort | grep ">"  | cut -c 2- | tr ' ' '\t' >> ${sample}.txt
    """
}

process MIRANDA_MERGE {
    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("${sample}.txt")

    script:
    """
    echo tmp_* > ${sample}.txt
    """
}

process MIRNA_TARGETS {
    tag "$sample"

    // container 'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2'

    input:
    tuple val(sample), path(targetscan), path(miranda), path(bed12)

    output:
    tuple val(sample), path("${sample}.mirna_targets.txt"), emit: results


    script:
    """
    ## reformat and sort miRanda, TargetScan outputs, convert to BED for overlaps.
    tail -n +2 $targetscan | sort -k1,1 -k4n | awk -v OFS="\t" '{print \$1, \$2, \$4, \$5, \$9}' | awk -v OFS="\t" '{print \$2, \$3, \$4, \$1, "0", \$5}' > targetscan.bed
    tail -n +2 $miranda | sort -k2,2 -k7n | awk -v OFS="\t" '{print \$2, \$1, \$3, \$4, \$7, \$8}' | awk -v OFS="\t" '{print \$2, \$5, \$6, \$1, \$3, \$4}' | sed 's/^[^-]*-//g' > miranda.bed

    ## intersect, consolidate miRanda, TargetScan information about miRs.
    ## -wa to output miRanda hits - targetscan makes it difficult to resolve duplicate miRNAs at MRE sites.
    bedtools intersect -a miranda.bed -b targetscan.bed -wa > ${sample}.mirnas.tmp
    bedtools intersect -a targetscan.bed -b miranda.bed | awk '{print \$6}' > mirna_type

    ## remove duplicate miRNA entries at MRE sites.
    ## strategy: sory by circs, sort by start position, sort by site type - the goal is to take the best site type (i.e rank site type found at MRE site).
    paste ${sample}.mirnas.tmp mirna_type | sort -k3,3 -k2n -k7r | awk -v OFS="\t" '{print \$4,\$1,\$2,\$3,\$5,\$6,\$7}' | awk -F "\t" '{if (!seen[\$1,\$2,\$3,\$4,\$5,\$6]++)print}' | sort -k1,1 -k3n > ${sample}.mirna_targets.tmp
    echo -e "circRNA\tmiRNA\tStart\tEnd\tScore\tEnergy_KcalMol\tSite_type" | cat - ${sample}.mirna_targets.tmp > ${sample}.mirna_targets.txt

    """
}