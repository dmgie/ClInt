import os
import csv
from utils import *

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/RawData/",
        "rna":"../../../../local_scratch/ClINT/working_files/deduped_vcfs/",
        "bam":"../../../../local_scratch/ClINT/working_files/deduped_bams/"
    }

    output_dir = "./output"
    os.makedirs(output_dir, exist_ok=True)
    gff_file = "../../../../local_scratch/ClINT/RawData/ref_genome.gff"
    gff_file_new = f"{output_dir}gff_new.gff"
    feature = "gene"
    feature_counts_output = f"{output_dir}/fc_output.txt"

    header_lines=[]
    gff_lines=[]
    variants = read_vcf_file(dir_dict)
    extracted_lines=[]
    matches=[]
    matches_count = 0
    total_count = 0
    
    with open(f'./{output_dir}/output.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Gene_ID', "Gene_Name", 'Start', 'End', 'Variant_Position', 'REF', "DNA_ALT", "RNA_ALT", "Counts"])
        print("Variants appearing at the same position in DNA and RNA\n")
        for idx, variant in enumerate(variants, start=1):
           
            if "rna" in variants[variant] and "dna" in variants[variant]:
                print(f"\nVariant {idx}")
                print(variant, "REF:", variants[variant]["ref"], "ALT DNA:", variants[variant]["dna"], "ALT RNA:", variants[variant]["rna"])
                extracted_lines=get_gff_lines(gff_file, variants[variant]["chromosome"], variant, 0, feature, [])
                gff_lines.append(extracted_lines)
                for line in extracted_lines:
                    values = line.split('\t')
                    gene_name=extract_attribute(values[8], 'Name=')
                    gene_id=extract_attribute(values[8], 'ID=')
                    start=values[3]
                    end=values[4]

                    tsv_writer.writerow([gene_id, gene_name, start, end, variant, variants[variant]["ref"], variants[variant]["dna"], variants[variant]["rna"]])

                    print(f"Found Variant at gene: {gene_name}, ID:{gene_id}, Start:{start}, End:{end}")
            
                print("Variants appearing in files:")
                for hit in variants[variant]["hits"]:
                    print(hit)
                    if hit not in matches:
                        matches.append(hit) 
                
                matches_count += 1     
            total_count += 1
        
        print(f"\n\nCounted {matches_count} dna/rna variant matches, {total_count} variants in total ({round(matches_count/total_count,2)}%).")

    create_gff(gff_file_new, get_gff_header(gff_file), gff_lines)


    bam_file_string = " ".join(str(dir_dict["bam"]+match) for match in matches)
    bam_file_string = remove_prefix_and_suffix(bam_file_string)

    #run_featurecounts(bam_file_string, gff_file_new, feature_counts_output, "gene")
    count = get_expression_count(feature_counts_output, "ACCFDFCE_00647", "../../../../local_scratch/ClINT/working_files/deduped_bams/dedup_snc_trimmed_Q12SA022AP_20211014141226__2109SHr1022_10_S22_R1_001.bam")

    

    






    # for variant in variants:
    #     if "rna" in variants[variant] and "dna" in variants[variant]:
    #         extracted_lines=get_gff_lines(gff_file, variants[variant]["chromosome"], variant, 0, feature, [])
    #         for line in extracted_lines:
    #             print("#########Count: ", line)
    #             for hit in variants[variant]["hits"]:
    #                 print(hit, "Expression Counts:", get_expression_count(feature_counts_output, {extract_attribute(values[8], 'ID=').replace("_gene", "")}, dir_dict["bam"]+remove_prefix_and_suffix(hit)))
    #                 if hit not in matches:
    #                     matches.append(hit) 
    
