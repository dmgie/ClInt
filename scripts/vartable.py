import pysam
import os

def read_vcf_file(types):
    variants_dict = {}

    for type in types:
        directory = os.fsencode(types[type])
        for file in os.listdir(directory):
            filename = os.fsencode(file).decode()
            if filename.endswith(".vcf"):
                print("############", filename,"############")
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

                        print(f"Chromosom: {chrom}")
                        print(f"Position: {pos}")
                        print(f"Referenz: {ref}")
                        print(f"Alternativen: {alt}")
                        print(f"Qualität: {qual}")
                        print(f"Filter-Status: {filter_status}")
                        print(f"Info: {info}")
                        print(f"Sample-Informationen: {sample_info}")
                        print("-" * 30)
                        
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
                              "hits": [filename]
                            }         
    return variants_dict


#def get_matching_variants():
    

def print_results(variants):
    matches_count = 0
    total_count = 0
    
    print("Variants appearing at the same position in DNA and RNA\n")
    for idx, variant in enumerate(variants, start=1):
           
        if "rna" in variants[variant] and "dna" in variants[variant]:
            print(f"\nVariant {idx}")
            print(variant, "REF:", variants[variant]["ref"], "ALT DNA:", variants[variant]["dna"], "ALT RNA:", variants[variant]["rna"])
            print("Variants appearing in files:")
            for hit in variants[variant]["hits"]:
                print(hit)
            matches_count += 1
            
        total_count += 1
        
    print(f"\n\nCounted {matches_count} dna/rna variant matches, {total_count} variants in total ({round(matches_count/total_count,2)}%).")
    

def get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines):
    with open(gff_file, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5: 
                    chr_name, source, feature_type, start, end, *_ = fields
                    if feature_type == feature:
                        if chr_name == str(chromosome): 
                            start_pos, end_pos = int(start), int(end)
                            if (start_pos - extend) <= position <= (end_pos + extend):
                                extracted_lines.append(line)
                                
        if not extracted_lines:
            extend += 10
            return get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines)
            
        
        print(f"Result found at extension {extend}")                                                
        return extracted_lines
    
    
    
           

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/vcfs/dna_vcf/",
        "rna":"../../../../local_scratch/ClINT/deduped_vcfs/"
    }

    variants = read_vcf_file(dir_dict)
    print_results(variants)
    
    # extend=0
    # extracted_lines=[]
    
    # extracted_lines=get_gff_lines("../../../../local_scratch/ClINT/vcfs/ref_genome.gff", 1, 1666529, extend, "gene", extracted_lines)
    # print("Results")
    
    # for line in extracted_lines:
    #         print(line)
