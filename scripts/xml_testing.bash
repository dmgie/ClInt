#!/usr/bin/env bash

#xmllint --xpath '/*' arguments.xml

#xmllint  --xpath '//programs/program[*]/arg[*]/text()' arguments.xml

#for program in "${!rna_programs[@]}" 
#do
#    xmllint  --xpath '//programs/program/[@type=$program]/arg[*]/text()' arguments.xml
#done


#test=$(xmllint --xpath "string(//body/value/@name)" test.xml)

declare -A rna_programs;         
declare -a rna_programs_order;
declare -A rna_categories;                                      
declare -a rna_categories_order;

declare -A rna_program_arguments;


 rna_programs["FastQC"]=0;                   rna_programs_order+=("FastQC")             rna_program_arguments["FastQC"]="$INPUT_GENOME_PATH $OUTPUT_DIR"
    rna_programs["MultiQC"]=0;                  rna_programs_order+=("MultiQC")
    rna_programs["fastp"]=0;                    rna_programs_order+=("fastp")              rna_program_arguments["fastp"]="$INPUT_GENOME_PATH $OUTPUT_DIR" 
    rna_programs["TrimGalore"]=0;               rna_programs_order+=("TrimGalore")
    rna_programs["SPAdes"]=0;                   rna_programs_order+=("SPAdes")
    rna_programs["Trinity"]=0;                  rna_programs_order+=("Trinity")
    rna_programs["STRING"]=0;                   rna_programs_order+=("STRING")
    rna_programs["minimap2"]=0;                 rna_programs_order+=("minimap2")           rna_program_arguments["minimap2"]="$HUMAN_REFERENCE_PATH $INPUT_GENOME_PATH $OUTPUT_DIR"
    rna_programs["HISAT"]=0;                    rna_programs_order+=("HISAT")
    rna_programs["STAR"]=0;                     rna_programs_order+=("STAR")
    rna_programs["BWA"]=0;                      rna_programs_order+=("BWA")
    rna_programs["FreeBayes"]=0;                rna_programs_order+=("FreeBayes")
    rna_programs["BCFTools"]=0;                 rna_programs_order+=("BCFTools")
    rna_programs["SnpEff"]=0;   


for program in "${!rna_programs[@]}" 
do
    echo $program
    xmlstarlet sel -t \
        --var prog_name="'$program'" \
        -v '//programs/program[@type = $prog_name]' \
        -nl arguments.xml
    echo $output
done

 