from itertools import count
from pickle import TRUE
from utils import *
import os
import csv
import argparse
import pandas as pd

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/RawData/",
        "rna":"../../../../local_scratch/ClINT/working_files/deduped_vcfs/",
        "bam":"../../../../local_scratch/ClINT/working_files/deduped_bams/",
        "out":"../../output",
        "gff":"../../../../local_scratch/ClINT/RawData/ref_genome.gff"
    }

    parser = argparse.ArgumentParser(description='Tool for comparative RNA/DNA variant analysis')

    parser.add_argument('--dna', required=False, help='DNA input folder -- required')
    parser.add_argument('--rna', required=False, help='RNA input folder -- required')
    parser.add_argument('--bam', required=False, help='BAM input folder -- required')
    parser.add_argument('--gff', required=False, help='GFF annotation file -- required')
    parser.add_argument('--out', required=False, help='Output folder -- required')
    parser.add_argument('--gff_filter', required=False, help='GFF filtering -- optional')

    args = parser.parse_args()

    for arg in vars(args):
        if getattr(args, arg) != None and arg in dir_dict.keys():
            dir_dict[arg] = getattr(args, arg)

    os.makedirs(dir_dict["out"], exist_ok=True)
   
    gff_file_new = f"{dir_dict['out']}/gff_new.gff"
    feature_counts_output = f"{dir_dict['out']}/fc_output.txt"
    feature_counts_tpm = f"{dir_dict['out']}/fc_output_tpm.txt"

    variant_positions = search_vcf_position_matches(dir_dict)
    bam_filenames = read_bam_filenames(dir_dict["bam"])
    gff_lines = []
    extracted_lines = []
    bam_file_string = []
    output_lines = []
    matches_count = 0
    total_count = 0

    gff_filtering = False

    ## Starting comparative .vcf analysis
    print("\nStarting comparative .vcf analysis\n")
    
    ## Create writeable output .tsv file, write header line
    with open(f'./{dir_dict["out"]}/output.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Gene_ID', "Gene_Name", 'Start', 'End', 'Variant_Position', 'Extension', 'REF', "DNA_ALT", "RNA_ALT", "Counts_Alt", "Counts"])
        
        for position in variant_positions:
            if "rna" in variant_positions[position] and "dna" in variant_positions[position]:
                
                print(f"Variant {matches_count+1}")
                print(position, "REF:", variant_positions[position]["ref"], "ALT DNA:", variant_positions[position]["dna"], "ALT RNA:", variant_positions[position]["rna"])

                extracted_lines, extend = get_gff_lines(dir_dict["gff"], variant_positions[position]["chromosome"], position, 0, "gene", [])
                gff_lines.append(extracted_lines)
 
                for line in extracted_lines:
                    values = line.split('\t')
                    start, end = values[3], values[4]

                    variant_line = {
                        "out": {
                            "gene_id":   extract_attribute(values[8], 'ID='),
                            "gene_name": extract_attribute(values[8], 'Name='),
                            "start":     start,
                            "end":       end,
                            "var_pos":   position,
                            "extension": get_extend_character(start, end, position),
                            "ref":       variant_positions[position]["ref"],
                            "dna_alt":   calculate_allele_percentages(variant_positions[position]["dna"]),
                            "rna_alt":   calculate_allele_percentages(variant_positions[position]["rna"])
                        },
                        "files": variant_positions[position]["hits"]
                    }
                    
                    for file in variant_positions[position]["hits"]:
                        if file not in bam_file_string:
                            bam_file_string.append(file)
                    
                    output_lines.append(variant_line) 
                matches_count += 1    
            total_count += 1


        ## Run featureCounts to get read counts per feature (gene)
        #### -> only on features that contain at least one variant position
        #### -> OPTIONAL: gff filtering - create new .gff containing only features with variants
        ####              runtime optimization / might imply normalization issues

        all_bam_filename_string = remove_prefix_and_suffix(" ".join(dir_dict["bam"] + bam for bam in bam_filenames))

        if gff_filtering is True:
            input_gff = create_gff(gff_file_new, get_gff_header(dir_dict["gff"]), gff_lines)
                
        else:
            input_gff = dir_dict["gff"]

        # run_featurecounts(all_bam_filename_string, input_gff, feature_counts_output, "gene")  

        ## Generate TSV output file
        #### -> File contains one position per line, matching in position
        #### -> feature holding position at both dna and rna sequences
        #### -> Only possible if dna and rna were aligned against same genomic refernce (same coordinate system)
        
        # Berechne die TPM-normierten Werte
        calculate_tpm(feature_counts_output, dir_dict["out"]) 

        for line in output_lines:
            counts = ""
            counts_alt = []

            if line["out"]["gene_id"] != 'No_Annotation_Found':
            
                alternative = get_alt_files(list(line["files"].keys()), dir_dict["bam"])
                
                ## Counts for variant matching files
                for file in line["files"]:
                    counts += line["files"][file] + ":" + \
                            str(get_expression_count(feature_counts_tpm, \
                            line["out"]["gene_id"].replace("_gene", ""), \
                            dir_dict["bam"] + remove_prefix_and_suffix(file))) + ","  
                            
                ## Counts for all other files, calculate average
                for file in alternative:
                    counts_alt.append( int(get_expression_count(feature_counts_tpm, \
                            line["out"]["gene_id"].replace("_gene", ""), \
                            dir_dict["bam"] + remove_prefix_and_suffix(file))))
                
                counts_alt_average = int(round(sum(counts_alt)/len(counts_alt),0))

            elif line["out"]["gene_id"] == 'No_Annotation_Found':
                counts = "/"
                counts_alt_average = "/"
                
            counts = counts.rsplit(',', maxsplit=1)[0]
            out_line = list(line["out"].values())+[counts_alt_average]+[counts]
            out_file.write('\t'.join(map(str, out_line)) + '\n')

    ## Console output, summary statistics
    print(f"\nCounted {matches_count} dna/rna variant position matches, {total_count} variant_positions in total ({round(matches_count/total_count,2)}).")
    print(f"Generated output file at './{dir_dict['out']}/output.tsv'")
