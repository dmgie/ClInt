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

    variant_positions = search_vcf_position_matches(dir_dict)
    bam_files = read_bam_filenames(dir_dict["bam"])
    gff_lines = []
    extracted_lines = []
    bam_file_string = []
    output_lines = []
    matches_count = 0
    total_count = 0

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

        

        ## Create GFF file of features containing variant_positions
        create_gff(gff_file_new, get_gff_header(dir_dict["gff"]), gff_lines)

        ## Run featureCounts only on features that contain at least one position
        # run_featurecounts(remove_prefix_and_suffix(" ".join(dir_dict["bam"]+bam for bam in bam_files)), gff_file_new, feature_counts_output, "gene")

        ## Generate TSV output file
        #### -> File contains one position per line, matching in position
        #### -> feature holding position at both dna and rna sequences
        #### -> Only possible if dna and rna were aligned against same genomic refernce (same coordinate system)

        for line in output_lines:
            counts = ""
            counts_alt = []
            
            alternative = get_alt_files(list(line["files"].keys()), dir_dict["bam"])
            
            ## Counts for variant matching files
            for file in line["files"]:
                counts += line["files"][file] + ":" + \
                          str(get_expression_count(feature_counts_output, \
                          line["out"]["gene_id"].replace("_gene", ""), \
                          dir_dict["bam"] + remove_prefix_and_suffix(file))) + ","  
                          
            ## Counts for all other files, calculate average
            for file in alternative:
                counts_alt.append( int(get_expression_count(feature_counts_output, \
                         line["out"]["gene_id"].replace("_gene", ""), \
                         dir_dict["bam"] + remove_prefix_and_suffix(file))))
            
            counts_alt_average = int(round(sum(counts_alt)/len(counts_alt),0))
                
            counts = counts.rsplit(',', maxsplit=1)[0]
            out_line = list(line["out"].values())+[counts_alt_average]+[counts]
            out_file.write('\t'.join(map(str, out_line)) + '\n')

    ## Console output, summary statistics
    print(f"\nCounted {matches_count} dna/rna variant position matches, {total_count} variant_positions in total ({round(matches_count/total_count,2)}).")
    print(f"Generated output file at './{dir_dict['out']}/output.tsv'")
