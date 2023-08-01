from itertools import count
from utils import *
import os
import csv

if __name__ == "__main__":

    dir_dict = {
        "dna":"../../../../local_scratch/ClINT/RawData/",
        "rna":"../../../../local_scratch/ClINT/working_files/deduped_vcfs/",
        "bam":"../../../../local_scratch/ClINT/working_files/deduped_bams/",
        "out":"../../output",
        "gff":"../../../../local_scratch/ClINT/RawData/ref_genome.gff"
    }

    os.makedirs(dir_dict["out"], exist_ok=True)
   
    gff_file_new = f"{dir_dict['out']}/gff_new.gff"
    feature_counts_output = f"{dir_dict['out']}/fc_output.txt"
    feature_counts_output_alt = f"{dir_dict['out']}/fc_output_alt.txt"

    variant_positions = read_vcf_file(dir_dict)
    gff_lines = []
    extracted_lines = []
    bam_file_string = []
    output_lines = []
    matches_count = 0
    total_count = 0
    
    ## Create writeable output .tsv file, write header line
    with open(f'./{dir_dict["out"]}/output.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Gene_ID', "Gene_Name", 'Start', 'End', 'Variant_Position', 'Extension', 'REF', "DNA_ALT", "RNA_ALT", "Counts"])
        
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

        
        bam_files = os.listdir(dir_dict["bam"])
        
        bam_files_filitered = []
        bam_files_total = []
        
        for file in bam_file_string:
            file = remove_prefix_and_suffix(file)
            bam_files_filitered.append(file)
            
        for file in bam_files:
            if file.startswith("dedup_snc_trimmed"):
                file = remove_prefix_and_suffix(file)
                bam_files_total.append(file)
        
        alt_files = set(bam_files_total) - set(bam_files_filitered)
        
            
            
        
        # for file in os.listdir(dir_dict["bam"]):
        #     found = False
        #     for bam_file in bam_file_string:
        #         filename = os.fsencode(file).decode()
        #         print(filename)
        #         if filename.endswith(".bam") and filename in bam_file:
        #             found = True
        #             break
        #         if not found:
        #             alt_files.append(filename)


        # counter = 0
        # filenames = []

        # for file in os.listdir(dir_dict["bam"]):
        #     filename = os.fsencode(file).decode()
        #     filenames.append(filename)
                

        # list1 = set(remove_prefix_and_suffix(str(bam_file_string)))
        # list2 = set(remove_prefix_and_suffix(str(filenames)))
        # #print("######", main_list)
                    # with pysam.VariantFile(directory+file) as vcf:
                    #     for record in vcf:
                    #         chrom = record.chrom
                    #         pos = record.pos
                    #         ref = record.ref
                    #         alt = ",".join(map(str, record.alts))
                        
                    #         if pos in variants_dict:
                    #             if type in variants_dict[pos]:
                    #                 variants_dict[pos][type] += ","+alt
                    #                 variants_dict[pos]["hits"].update({filename:alt})
                    #             else:
                    #                 variants_dict[pos][type] = alt
                    #                 variants_dict[pos]["hits"] = {filename:alt}
                    #         else:
                    #             variants_dict[pos] = {
                    #             "ref": ref,
                    #             type: alt,
                    #             "hits": {filename:alt},
                    #             "chromosome":chrom
                    #             }       
        # print("#####", found, len(found))  
        # print("#####", alt_files, len(alt_files))  
        # print("1", list1)
        # print("2", list2)
        
        ## Console output, summary statistics
        print(f"\nCounted {matches_count} dna/rna variant position matches, {total_count} variant_positions in total ({round(matches_count/total_count,2)}).")
        print(f"Generated output file at './{dir_dict['out']}/output.tsv'")

        ## Create GFF file of features containing variant_positions
        create_gff(gff_file_new, get_gff_header(dir_dict["gff"]), gff_lines)

        ## Run featureCounts only on features that contain at least one position
        run_featurecounts(remove_prefix_and_suffix(" ".join(dir_dict["bam"]+bam for bam in bam_file_string)), gff_file_new, feature_counts_output, "gene")
        
        ## Run featureCounts for files thtat do not contain variant
        run_featurecounts(remove_prefix_and_suffix(" ".join(dir_dict["bam"]+bam for bam in alt_files)), gff_file_new, feature_counts_output_alt, "gene")

        ## Generate TSV output file
        #### -> File contains one position per line, matching in position
        #### -> feature holding position at both dna and rna sequences
        #### -> Only possible if dna and rna were aligned against same genomic refernce (same coordinate system)

        for line in output_lines:
            counts = ""
            for file in line["files"]:
                counts += line["files"][file] + ":" + \
                          str(get_expression_count(feature_counts_output, \
                          line["out"]["gene_id"].replace("_gene", ""), \
                          dir_dict["bam"] + remove_prefix_and_suffix(file))) + ","  
                
            counts = counts.rsplit(',', maxsplit=1)[0]
            out_line = list(line["out"].values())+[counts]
            out_file.write('\t'.join(map(str, out_line)) + '\n')
