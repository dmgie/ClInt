// This would contain all the needed things for circRNA related workflow

include { STAR as STAR_PAIRED } from "./mapping.nf"
include { STAR as STAR_MATE1 } from "./mapping.nf"
include { STAR as STAR_MATE2 } from "./mapping.nf"
include { MIRNA_PREDICTION } from "./mirna.nf"

circ = true;

workflow CIRCRNA {
         take:
            reads      // Sample IDs + Paired end reads
            ref_file
            annotation
         // Preperation (Paired, and mate-independent mapping)

         main:
                ref_idx = Channel.fromPath("${params.star_index}/*").collect()

         read1 = reads.map { name, reads -> tuple(name, [reads[0]]) }
         read2 = reads.map { name, reads -> tuple(name, [reads[1]]) }

         STAR_PAIRED(reads, ref_idx, circ)
         STAR_MATE1(read1, ref_idx, circ)
         STAR_MATE2(read2, ref_idx, circ)

         prep_dcc = STAR_PAIRED.out.junction.map { name, junction -> tuple(name, junction) }
                                      .join(STAR_MATE1.out.junction, remainder: true)
                                      .join(STAR_MATE2.out.junction, remainder: true)
         prep_dcc.view()
         

         DCC(prep_dcc, ref_file, annotation)
         DCC.out.view()
         DCC_FILTER(DCC.out.txt, 2) // 2 = min number of bsj_reads to filter by

         // Annotation
         ch_biotypes = Channel.fromPath("${projectDir}/bin/unwanted_biotypes.txt")
         ANNOTATION(DCC_FILTER.out.results, annotation, ch_biotypes.collect(), 200) // 200 = exon boundary

         // Sequence Extraction
         FASTA(ANNOTATION.out.bed, ref_file)


         // miRNA Prediction
         MIRNA_PREDICTION(FASTA.out.analysis_fasta, ANNOTATION.out.bed)
         
         // Identification (DCC)
         // Quantification (CIRIquant?)
         // Characterisation (sequence extraction e.g FUCHS)
}


process DCC {
    tag "${sample}"

    input:
    tuple val(sample), path(paired), path(mate1, stageAs: "mate1.Chimeric.out.junction"), path(mate2, stageAs: "mate2.Chimeric.out.junction")
    path ref_fasta
    path annotation
    // path repeats

    output:
    tuple val(sample), path("${sample}.txt"), emit: txt

    script:
    // NOTE: The "-R ${repeats}" doesn't seem to work well, causing some errors with HTSeq which DCC uses for GTF file parsing
    """

    mkdir ${sample} && mv ${sample}_Chimeric.out.junction ${sample} && printf "${sample}/${sample}_Chimeric.out.junction" > samplesheet
    mkdir ${sample}_mate1 && mv mate1.Chimeric.out.junction ${sample}_mate1 && printf "${sample}_mate1/mate1.Chimeric.out.junction" > mate1file
    mkdir ${sample}_mate2 && mv mate2.Chimeric.out.junction ${sample}_mate2 && printf "${sample}_mate2/mate2.Chimeric.out.junction" > mate2file

    DCC @samplesheet -mt1 @mate1file -mt2 @mate2file -D -an ${annotation} -Pi -ss -F -M -Nr 1 1 -fg -A ${ref_fasta} -T ${task.cpus}

    awk '{print \$6}' CircCoordinates >> strand
    paste CircRNACount strand | tail -n +2 | awk -v OFS="\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${sample}.txt

    """

}


process DCC_FILTER {
    tag "${sample}"
    label 'process_single'

    input:
    tuple val(sample), path(txt)
    val bsj_reads

    output:
    tuple val(sample), path("${sample}_dcc_circs.bed"), emit: results
    tuple val(sample), path("${sample}_dcc.bed")      , emit: matrix

    script:
    """
    awk '{if(\$5 >= ${bsj_reads}) print \$0}' ${sample}.txt > ${sample}_dcc.filtered
    awk -v OFS="\t" '{\$2-=1;print}' ${sample}_dcc.filtered > ${sample}_dcc.bed
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${sample}_dcc.bed > ${sample}_dcc_circs.bed
    """
}

process ANNOTATION {
    tag "${sample}"
    container "" // Doesn't need one since it uses gtfTogenePred

    input:
    tuple val(sample), path(bed)
    path gtf
    path biotypes
    val exon_boundary


    output:
    tuple val(sample), path("${sample}.bed"), emit: bed
    path("*.log")                         , emit: log
    
    script:
    """
    grep -vf $biotypes $gtf > filt.gtf
    mv $bed circs.bed

    annotate_outputs.sh $exon_boundary &> ${sample}.log
    mv master_bed12.bed ${sample}.bed.tmp

    awk -v FS="\t" '{print \$11}' ${sample}.bed.tmp > mature_len.tmp
    awk -v FS="," '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' mature_len.tmp > mature_length

    paste ${sample}.bed.tmp mature_length > ${sample}.bed
    """
}

process FASTA {
    tag "${sample}"

    // own container was missing "sed" command
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2'


    input:
    tuple val(sample), path(bed)
    path fasta

    output:
    tuple val(sample), path("${sample}.fa"), emit: analysis_fasta
    path("${sample}.fasta")              , emit: publish_fasta

    script:
    """
    ## FASTA sequences (bedtools does not like the extra annotation info - split will not work properly)
    cut -d\$'\t' -f1-12 ${sample}.bed > bed12.tmp
    bedtools getfasta -fi $fasta -bed bed12.tmp -s -split -name > circ_seq.tmp

    ## clean fasta header
    grep -A 1 '>' circ_seq.tmp | cut -d: -f1,2,3 > ${sample}.fa && rm circ_seq.tmp

    ## add backsplice sequence for miRanda Targetscan, publish canonical FASTA to results.
    rm $fasta
    bash ${workflow.projectDir}/bin/backsplice_gen.sh ${sample}.fa
    """
}

