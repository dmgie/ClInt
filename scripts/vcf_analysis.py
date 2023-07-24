import pysam
import os

def read_vcf_file(types):
    variants_dict = {}

    for type in types:
        directory = os.fsencode(types[type])
        for file in os.listdir(directory):
            filename = os.fsencode(file).decode()
            if filename.endswith(".vcf"):
                with pysam.VariantFile(directory+file) as vcf:
                    for record in vcf:
                        chrom = record.chrom
                        pos = record.pos
                        ref = record.ref
                        alt = ",".join(map(str, record.alts))
                        qual = record.qual
                        filter_status = record.filter
                        info = record.info
                        sample_info = record.samples

                        # print(f"Chromosom: {chrom}")
                        # print(f"Position: {pos}")
                        # print(f"Referenz: {ref}")
                        # print(f"Alternativen: {alt}")
                        # print(f"Qualität: {qual}")
                        # print(f"Filter-Status: {filter_status}")
                        # print(f"Info: {info}")
                        # print(f"Sample-Informationen: {sample_info}")
                        # print("-" * 30)
                        
                        if pos in variants_dict:
                            if type in variants_dict[pos]:
                                variants_dict[pos][type] += ","+alt
                                variants_dict[pos]["hits"] += [filename]
                            else:
                                variants_dict[pos][type] = alt
                                variants_dict[pos]["hits"] = [filename]
                        else:
                            variants_dict[pos] = {
                             "ref": ref,
                              type: alt,
                              "hits": [filename],
                              "chromosome":chrom
                            }         
    return variants_dict
    

def print_results(variants, gff_file, feature):
    matches_count = 0
    total_count = 0
    
    print("Variants appearing at the same position in DNA and RNA\n")
    for idx, variant in enumerate(variants, start=1):
           
        if "rna" in variants[variant] and "dna" in variants[variant]:
            print(f"\nVariant {idx}")
            print(variant, "REF:", variants[variant]["ref"], "ALT DNA:", variants[variant]["dna"], "ALT RNA:", variants[variant]["rna"])
            extracted_lines=get_gff_lines(gff_file, variants[variant]["chromosome"], variant, 0, feature, [])
            for line in extracted_lines:
                values = line.split('\t')
                print(f"Found Variant at gene: {extract_attribute(values[8], 'Name=')}, ID:{extract_attribute(values[8], 'ID=')}, Start:{values[3]}, End:{values[4]}")
            
            print("Variants appearing in files:")
            for hit in variants[variant]["hits"]:
                print(hit)
            matches_count += 1

            
            
        total_count += 1
        
    print(f"\n\nCounted {matches_count} dna/rna variant matches, {total_count} variants in total ({round(matches_count/total_count,2)}%).")
    

def get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines):
    
    #header_lines=[]
    with open(gff_file, 'r+') as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5: 
                    chr_name, source, feature_type, start, end, *_ = fields
                    if feature_type == feature:
                        if chr_name == str(chromosome): 
                            start_pos, end_pos = int(start), int(end)
                            if (start_pos - extend) <= position <= (end_pos + extend):
                                extracted_lines.append(line.strip())
            #else:
             #   header_lines.append(line.strip())
                                
        if not extracted_lines:
            extend += 10
            return get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines)
            
            
    # try:
    #     mode = 'a' if os.path.exists("./gff_new.gff") else 'w'
    #     with open("gff_new.gff", mode) as outfile:
    #         for line in header_lines:
    #             outfile.write(f"{line}\n")
    #         for line in extracted_lines:
    #             outfile.write(f"{line}\n")
    # except Exception as e:
    #     print("Error at writing .gff file")   
        
        
    print(f"Result found at extension {extend}")                                                
    return extracted_lines
    

def extract_attribute(input_string, query):
    start = input_string.find(query)
    if start == -1:
        return None

    start += len("Name=")
    end = input_string.find(";", start)
    if end == -1:
        return None

    return input_string[start:end]
    
    
           

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/RawData-RNA/",
        "rna":"../../../../local_scratch/ClINT/working_files/deduped_vcfs/"
    }
    
    gff_file = "../../../../local_scratch/ClINT/RawData-RNA/ref_genome.gff"
    gff_file_new = "./gff_new.gff"
    feature = "gene"
    
    header_lines=[]
    
    # mode = 'a' if os.path.exists(gff_file_new) else 'w'
    # with open(gff_file_new, mode) as infile: 
    #     for line in infile:   
    #         header_lines.append(line.strip() + "\n")
    

    variants = read_vcf_file(dir_dict)
    print_results(variants, gff_file, feature)
    
