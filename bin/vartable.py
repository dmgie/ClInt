from utils import *
import os
import csv
import sys

def run_vartable(dna_vcf, rna_vcf, bam, output_dir, gff, patient, gff_filtering, snpEff, agreement, dna_startswith, rna_startswith):

    dir_dict = {
        "dna":dna_vcf,
        "rna":rna_vcf,
        "bam":bam,
        "out":output_dir,
        "gff":gff
    }
    
    print(f"\nStarting comparative .vcf analysis on files from patient: {patient}\n")
    print(f"RNA files from {dir_dict['rna']}:", rna_startswith)
    print(f"DNA files from {dir_dict['dna']}:", dna_startswith)

    #### Initialization of variables

    variant_positions = search_vcf_position_matches(dir_dict, dna_startswith, rna_startswith, snpEff)

    if variant_positions == None:
        print("\n>>>>>No matching dna/rna files found - quitting vartable analysis<<<<<")
        sys.exit()

    patient = patient.replace(" ", "")
    dir_dict["out"] += f"/{patient}"
    os.makedirs(f"{dir_dict['out']}", exist_ok=True)
    os.makedirs(f"{dir_dict['out']}/featureCounts", exist_ok=True)
   
    gff_file_new = f"{dir_dict['out']}/gff_new.gff"
    feature_counts_output = f"{dir_dict['out']}/featureCounts/fc_output.txt"
    feature_counts_tpm = f"{dir_dict['out']}/featureCounts/fc_output_tpm.txt"

    bam_filenames = read_bam_filenames(dir_dict["bam"], rna_startswith)
    gff_lines = []
    extracted_gff_lines = []
    bam_file_string = []
    output_lines = []
    agreement = {}
    variants_dict = {}

    matching_count = 0 
    total_count = 0
    dna_count = 0
    rna_count = 0

    ## Starting comparative .vcf analysis
    ## Create writeable output .tsv file, write header line
    with open(f'./{dir_dict["out"]}/{patient}_vartable_output.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Gene_ID', "Gene_Name", "Chromosome", 'Start', 'End', 'Variant_Position', 'Extension', 'REF', "DNA_ALT", "RNA_ALT", "Counts_Alt", "Counts_TPM", "Annotation"])
        
        ## Iterate over all found variant positions
        for position in variant_positions:
            
            ## No variant in RNA but variant in DNA at position#
            if not "rna" in variant_positions[position] and "dna" in variant_positions[position]:
                dna_count += 1

            ## Variant in RNA but no variant in DNA at position#
            if "rna" in variant_positions[position] and not "dna" in variant_positions[position]:
                rna_count += 1

            ## Matching variants in both RNA and DNA at position#
            if "rna" in variant_positions[position] and "dna" in variant_positions[position]:
                ## Console output
                print(f"Variant {matching_count+1}")
                print(position, "REF:", variant_positions[position]["ref"], "ALT DNA:", variant_positions[position]["dna"], "ALT RNA:", variant_positions[position]["rna"])

                ## DNA / RNA match counting
                dna_count += 1
                rna_count += 1
                matching_count += 1
                for hit in variant_positions[position]["hits"]:
                    if hit in agreement:
                        agreement[hit] += 1
                    else:
                        agreement[hit] = 1

                extracted_gff_lines, extend = get_gff_lines(dir_dict["gff"], variant_positions[position]["chromosome"], position, 0, "gene", [])
                gff_lines.append(extracted_gff_lines)
 
                for line in extracted_gff_lines:
                    values = line.split('\t')
                    start, end = values[3], values[4]

                    variant_line = {
                        "out": {
                            "gene_id":   extract_attribute(values[8], 'ID='),
                            "gene_name": extract_attribute(values[8], 'Name='),
                            "chromosome":variant_positions[position]["chromosome"],
                            "start":     start,
                            "end":       end,
                            "var_pos":   position,
                            "extension": get_extend_character(start, end, position),
                            "ref":       variant_positions[position]["ref"],
                            "dna_alt":   calculate_allele_percentages(variant_positions[position]["dna"]),
                            "rna_alt":   calculate_allele_percentages(variant_positions[position]["rna"]),
                            "counts_alt":[],
                            "counts_tpm":"",
                            "annotation":variant_positions[position]["annotation"][0]
                        },
                        "files": variant_positions[position]["hits"]
                    }
                    
                    for file in variant_positions[position]["hits"]:
                        if file not in bam_file_string:
                            bam_file_string.append(file)
                    
                    output_lines.append(variant_line)  
            total_count += 1


        ## Create DNA / RNA agreement output file
        if agreement: create_agreement(agreement, dna_count, rna_count, matching_count, dir_dict["out"], patient)

        ## Run featureCounts to get read counts per feature (gene)
        #### -> only on features that contain at least one variant position
        #### -> OPTIONAL: gff filtering - create new .gff containing only features with variants
        ####              runtime optimization / might imply normalization issues

        all_bam_filename_string = remove_prefix_and_suffix(" ".join(dir_dict["bam"] + bam for bam in bam_filenames))

        if gff_filtering == True:
            print("Run featureCounts using gff_filtering")
            input_gff = create_gff(gff_file_new, get_gff_header(dir_dict["gff"]), gff_lines)
                
        else:
            print("Run featureCounts without using gff_filtering")
            input_gff = dir_dict["gff"]

        run_featurecounts(all_bam_filename_string, input_gff, feature_counts_output, "gene")  

        ## Calculate TPM-normalized values (transcripts-per-million)
        #### -> Make samples directly comparable
        #### -> Print to featureCoults-like file

        calculate_tpm(feature_counts_output, dir_dict["out"]+"/featureCounts") 

        ## Generate TSV output file
        #### -> File contains one position per line, matching in position
        #### -> feature holding position at both dna and rna sequences
        #### -> Only possible if dna and rna were aligned against same genomic refernce (same coordinate system)
        
        ## Iterate over all output lines
        for line in output_lines:

            ## Get read information from featureCounts output file
            if line["out"]["gene_id"] != 'No_Annotation_Found':
            
                # print("######", line, output_lines)
                alternative = get_alt_files(list(line["files"].keys()), dir_dict["bam"])
                
                ## Counts for variant matching files
                for file in line["files"]:
                    line["out"]["counts_tpm"] += line["files"][file] + ":" + \
                            str(get_expression_count(feature_counts_tpm, \
                            line["out"]["gene_id"].replace("_gene", ""), \
                            dir_dict["bam"] + remove_prefix_and_suffix(file))) + ","  
                            
                ## Counts for all other files, calculate average
                for file in alternative:
                    line["out"]["counts_alt"].append( int(get_expression_count(feature_counts_tpm, \
                            line["out"]["gene_id"].replace("_gene", ""), \
                            dir_dict["bam"] + remove_prefix_and_suffix(file))))
                    
                pos = line["out"]["var_pos"]
                chrom = line["out"]["chromosome"]
                
                variant = line["out"]["var_pos"]
                    
                variants_dict[chrom+":"+str(pos)] = {"variant": line["out"]["rna_alt"]}
                variants_dict[chrom+":"+str(pos)].update({"counts":  line["out"]["counts_tpm"].rsplit(',', maxsplit=1)[0]})
                
            ## Add / character if no matching annotation feature is found 
            elif line["out"]["gene_id"] == 'No_Annotation_Found':
                line["out"]["counts_tpm"] = "/"
                line["out"]["counts_alt"] = "/"

            ## Add reads counts to each output line (if available)  
            ## Concatenate output strings (variant information & counts)
            line["out"]["counts_tpm"] = line["out"]["counts_tpm"].rsplit(',', maxsplit=1)[0]
            out_line = list(line["out"].values())
            out_file.write('\t'.join(map(str, out_line)) + '\n')

    ## Console output, summary statistics
    if matching_count != 0:
        print(f"\nCounted {matching_count} dna/rna variant position matches, {total_count} variant_positions in total ({round(matching_count/total_count,2)}).")
        print(f"Generated meta variant comparison output file at './{dir_dict['out']}/output.tsv'")

    if agreement: print(f"Generated dna/rna agreement output file at './{dir_dict['out']}/agreement.tsv'")
    
    return variants_dict