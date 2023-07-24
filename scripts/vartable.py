import os
from utils import *

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/RawData-RNA/",
        "rna":"../../../../local_scratch/ClINT/working_files/deduped_vcfs/"
    }
    
    gff_file = "../../../../local_scratch/ClINT/RawData-RNA/ref_genome.gff"
    gff_file_new = "./gff_new.gff"
    feature = "gene"

    header_lines=[]
    gff_lines=[]
    variants = read_vcf_file(dir_dict)
    matches_count = 0
    total_count = 0
    
    print("Variants appearing at the same position in DNA and RNA\n")
    for idx, variant in enumerate(variants, start=1):
           
        if "rna" in variants[variant] and "dna" in variants[variant]:
            print(f"\nVariant {idx}")
            print(variant, "REF:", variants[variant]["ref"], "ALT DNA:", variants[variant]["dna"], "ALT RNA:", variants[variant]["rna"])
            extracted_lines=get_gff_lines(gff_file, variants[variant]["chromosome"], variant, 0, feature, [])
            gff_lines.append(extracted_lines)
            for line in extracted_lines:
                values = line.split('\t')
                print(f"Found Variant at gene: {extract_attribute(values[8], 'Name=')}, ID:{extract_attribute(values[8], 'ID=')}, Start:{values[3]}, End:{values[4]}")
            
            print("Variants appearing in files:")
            for hit in variants[variant]["hits"]:
                print(hit)
            matches_count += 1     
        total_count += 1
        
    print(f"\n\nCounted {matches_count} dna/rna variant matches, {total_count} variants in total ({round(matches_count/total_count,2)}%).")
    create_gff(gff_file_new, get_gff_header(gff_file), gff_lines)
    
