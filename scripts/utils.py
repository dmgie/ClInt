import os

def extract_attribute(input_string, query):
    start = input_string.find(query)
    if start == -1:
        return None

    start += len("Name=")
    end = input_string.find(";", start)
    if end == -1:
        return None

    return input_string[start:end]

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