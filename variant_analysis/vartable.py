import os
import csv
import numpy as np
from utils import *

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/RawData/",
        "rna":"../../../../local_scratch/ClINT/working_files/deduped_vcfs/",
        "bam":"../../../../local_scratch/ClINT/working_files/deduped_bams/"
    }

    output_dir = "../../output"
    os.makedirs(output_dir, exist_ok=True)
    gff_file = "../../../../local_scratch/ClINT/RawData/ref_genome.gff"
    gff_file_new = f"{output_dir}/gff_new.gff"
    feature = "gene"
    feature_counts_output = f"{output_dir}/fc_output.txt"

    header_lines=[]
    gff_lines=[]
    variants = read_vcf_file(dir_dict)
    extracted_lines=[]
    matches_count = 0
    total_count = 0
    bam_file_string = []

    output_line=[]
    
    with open(f'./{output_dir}/output.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Gene_ID', "Gene_Name", 'Start', 'End', 'Variant_Position', 'REF', "DNA_ALT", "RNA_ALT", "Counts"])
        print("Variants appearing at the same position in DNA and RNA\n")
        
        for idx, variant in enumerate(variants, start=1):
            if "rna" in variants[variant] and "dna" in variants[variant]:
                print("#####", variant, variants[variant])
                print(f"\nVariant {idx}")
                print(variant, "REF:", variants[variant]["ref"], "ALT DNA:", variants[variant]["dna"], "ALT RNA:", variants[variant]["rna"])
                extracted_lines = get_gff_lines(gff_file, variants[variant]["chromosome"], variant, 0, feature, [])
                gff_lines.append(extracted_lines)
                  
                for line in extracted_lines:
                    values = line.split('\t')

                    variant_line = {
                        "out": {
                            "gene_id":   extract_attribute(values[8], 'ID='),
                            "gene_name": extract_attribute(values[8], 'Name='),
                            "start":     values[3],
                            "end":       values[4],
                            "var_pos":   variant,
                            "ref":       variants[variant]["ref"],
                            "dna_alt":   calculate_allele_percentages(variants[variant]["dna"]),
                            "rna_alt":   calculate_allele_percentages(variants[variant]["rna"])
                        },
                        "files": variants[variant]["hits"]
                    }
                    
                    for file in variants[variant]["hits"]:
                        if file not in bam_file_string:
                            bam_file_string.append(file)
                    
                    output_line.append(variant_line) 
                matches_count += 1     
            total_count += 1
        
        print(f"\n\nCounted {matches_count} dna/rna variant matches, {total_count} variants in total ({round(matches_count/total_count,2)}%).")

        create_gff(gff_file_new, get_gff_header(gff_file), gff_lines)
        #run_featurecounts(remove_prefix_and_suffix(" ".join(dir_dict["bam"]+bam for bam in bam_file_string)), gff_file_new, feature_counts_output, "gene")


        for row in output_line:
            # print(row)
            counts = ""
            for file in row["files"]:
                counts += str(get_expression_count(feature_counts_output, row["out"]["gene_id"].replace("_gene", ""), dir_dict["bam"] + remove_prefix_and_suffix(file)))+","  
            counts=counts.rsplit(',', maxsplit=1)[0]
            out_line=list(row["out"].values())+[counts]
            out_file.write('\t'.join(map(str, out_line)) + '\n')
