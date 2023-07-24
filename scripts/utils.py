import os
import pysam

def extract_attribute(input_string, query):
    start = input_string.find(query)
    if start == -1:
        return None

    start += len("Name=")
    end = input_string.find(";", start)
    if end == -1:
        return None

    return input_string[start:end]

def create_gff(gff_file_new, header_lines, extracted_lines):
    with open(gff_file_new, 'w') as f:
        # Schreibe die Header-Zeilen in die Datei
        for header_line in header_lines:
            f.write(f"{header_line[0]}\n")
        
        # Schreibe die Inhaltszeilen in die Datei
        for line in extracted_lines:
            gff_line = "\t".join(line)
            f.write(f"{gff_line}\n")

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

                        # print_variants(chrom, pos, ref, alt, qual, filter_status, info, sample_info)
                        
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

def print_variants():
    print(f"Chromosom: {chrom}")
    print(f"Position: {pos}")
    print(f"Referenz: {ref}")
    print(f"Alternativen: {alt}")
    print(f"Qualität: {qual}")
    print(f"Filter-Status: {filter_status}")
    print(f"Info: {info}")
    print(f"Sample-Informationen: {sample_info}")
    print("-" * 30)


def get_gff_header(gff_file):
    header_lines=[]
    with open(gff_file, 'r+') as infile:
        for line in infile:
            if line.startswith('#'):
               header_lines.append(line.strip().split("\t"))
    return header_lines

def get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines):
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
                                
        if not extracted_lines:
            extend += 10
            return get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines)
            
    print(f"Result found at extension {extend}")                                                
    return extracted_lines