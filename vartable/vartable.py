from utils import *
import os
import re
import csv
import multiprocessing
from tqdm import tqdm

#### Initialization of global variables

manager = multiprocessing.Manager()

variant_positions = manager.dict()
gff_lines = []
extracted_gff_lines = []
bam_file_string = []
output_lines = manager.list()
agreement_dict = manager.dict()
variants_dict = manager.dict()

matching_count = 0
total_count = 0
dna_count = 0
rna_count = 0
both_dna_and_rna_count = 0
rna_prefixes = []

gff_parsed = []
feature_counts_tpm = ""
dir_dict = {}

feature = ""

lock = manager.Lock()

def run_vartable(dna_vcf, rna_vcf, bam, output_dir, gff, patient, gff_filtering, snpEff, agreement, dna_startswith, rna_startswith, num_processes, _feature, condition):

    global variant_positions
    global gff_lines
    global extracted_gff_lines
    global bam_file_string
    global output_lines
    global agreement_dict
    global variants_dict
    global feature_counts_tpm
    global rna_prefixes

    global matching_count
    global total_count
    global dna_count
    global rna_count
    global both_dna_and_rna_count
    global dir_dict
    global gff_parsed
    global feature
    
    feature = _feature
    condition = condition.replace(" ", "")

    dir_dict = {
        "dna":dna_vcf,
        "rna":rna_vcf,
        "bam":bam,
        "out":output_dir,
        "gff":gff
    }
    
    rna_prefixes = rna_startswith
    
    print(f"\nStarting comparative .vcf analysis on files from patient: {patient}\n")
    print(f"RNA files from {dir_dict['rna']}:", rna_startswith)
    print(f"DNA files from {dir_dict['dna']}:", dna_startswith)

    print(f"\n################# PARSE VARIANT POSITIONS FOR FEATURE {feature} FOR CONDITIONS {condition}\n")
    variant_positions = search_vcf_position_matches(dir_dict, dna_startswith, rna_startswith, snpEff, num_processes)

    ## Quit program if no matching variants were found
    if variant_positions == None:
        print("\n>>>>>Some of the input files were not found - quitting vartable analysis<<<<<")
        return {}
    
    ## Quit program if no matching variants were found
    if total_positions == 0:
        print("\n>>>>>No matching dna/rna files found - quitting vartable analysis<<<<<")
        return {}
    
    total_positions = len(variant_positions)

    ## Pre-calculate dna/rna variant matching rate
    #### -> Only consider matching dna/rna variant positions
    #### -> Discard unmatching variant positions
     
    for position, data in list(variant_positions.items()):
        
        ## DNA / RNA match
        if data.get("dna") is not None and data.get("rna") is not None:
            matching_count += 1
            dna_count += 1
            rna_count += 1

        ## DNA variant only
        elif data.get("dna") is not None and data.get("rna") is None:
            dna_count += 1
            del variant_positions[position]

        ## RNA variant only
        elif data.get("dna") is None and data.get("rna") is not None:
            rna_count += 1
            del variant_positions[position]
        
        total_count += 1

    percentage_both = (matching_count / total_positions) * 100

    print(f"matching percentage: {percentage_both:.2f}%\nmatching count: {matching_count}\ntotal positions: {total_positions}\n")


    print("\n##################### PARSE GFF FILE ######################")
    gff_parsed = parse_gff(dir_dict["gff"])
    print("\n################## FINISHED GFF PARSING ###################\n\n")

    patient = patient.replace(" ", "")
    dir_dict["out"] += f"/{patient}"
    os.makedirs(f"{dir_dict['out']}", exist_ok=True)
    os.makedirs(f"{dir_dict['out']}/featureCounts", exist_ok=True)
   
    gff_file_new = f"{dir_dict['out']}/gff_new.gff"
    feature_counts_output = f"{dir_dict['out']}/featureCounts/fc_output.txt"
    feature_counts_tpm = f"{dir_dict['out']}/featureCounts/fc_output_tpm.txt"


    ## Starting comparative .vcf analysis
    #### -> Processing of variant positions
    #### -> Create match entry if variant appears in both dna and rna
    #### -> Fully parallelized with multiprocessing
    
    
    print("\n################## PROCESSING POSITIONS ###################")
    with multiprocessing.Pool(processes=num_processes) as pool:
        list(tqdm(pool.imap(process_position, variant_positions), total=len(variant_positions)))

    
    ## Create DNA / RNA agreement output file
    if agreement: create_agreement(agreement_dict, dna_count, rna_count, matching_count, dir_dict["out"], patient, condition)

    ## Run featureCounts to get read counts per feature (gene)
    #### -> only on features that contain at least one variant position
    #### -> OPTIONAL: gff filtering - create new .gff containing only features with variants
    ####              runtime optimization / might imply normalization issues


    print("#####Finished variant matching, proceeding with featureCounts#####")
    bam_filenames = read_bam_filenames(dir_dict["bam"], rna_startswith)
    all_bam_filename_string = remove_prefix_and_suffix(" ".join(dir_dict["bam"] + bam for bam in bam_filenames))

    if gff_filtering == True:
        print("Run featureCounts using gff_filtering")
        input_gff = create_gff(gff_file_new, get_gff_header(dir_dict["gff"]), gff_lines)
            
    else:
        print("Run featureCounts without using gff_filtering")
        input_gff = dir_dict["gff"]

    run_featurecounts(all_bam_filename_string, input_gff, feature_counts_output, "gene", num_processes)  

    ## Calculate TPM-normalized values (transcripts-per-million)
    #### -> Make samples directly comparable
    #### -> Print to featureCoults-like file

    calculate_tpm(feature_counts_output, dir_dict["out"]+"/featureCounts") 

    ## Generate TSV output file
    #### -> File contains one position per line, matching in position
    #### -> feature holding position at both dna and rna sequences
    #### -> Only possible if dna and rna were aligned against same genomic refernce (same coordinate system)
    #### -> Fully parallelized with multiprocessing
    
    
    print("Starting count extraction from featureCounts output file")

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = list(tqdm(pool.imap(extract_counts, output_lines), total=len(output_lines)))

        # Create writeable output .tsv file 
        ### -> write header line
        ### -> write each result line below
        
        with open(f'./{dir_dict["out"]}/{condition}_{patient}_vartable_output.tsv', 'wt') as out_file:
            writer = csv.writer(out_file, delimiter="\t")
            writer.writerow(['Gene_ID', "Gene_Name", "Chromosome", 'Start', 'End', 'Variant_Position', 'Extension', 'REF', "DNA_ALT", "RNA_ALT", "Counts_Alt", "Counts_TPM", "Annotation"])
            for data in results:
                out_file.write(data)

    ## Console output, summary statistics
    if matching_count != 0:
        print(f"\nCounted {matching_count} dna/rna variant position matches, {total_count} variant_positions in total ({round(matching_count/total_count,2)}).")
        print(f"Generated meta variant comparison output file at './{dir_dict['out']}/output.tsv'")

    if agreement: print(f"Generated dna/rna agreement output file at './{dir_dict['out']}/agreement.tsv'")
    
    ## Return results (dict of matching variants) to vartable file handler
    return variants_dict

def extract_counts(line):

    ## Get read information from featureCounts output file
    
    global feature_counts_tpm
    global dir_dict
    global variants_dict
    global rna_prefixes
    

    if line["out"]["gene_id"] != 'No_Annotation_Found':
        
        alternative = filter_filenames_by_prefixes(line["files"],rna_prefixes)
        
        ## Counts for matching files
        for file in line["files"]:
            line["out"]["counts_tpm"] += int(get_expression_count(feature_counts_tpm, line["out"]["gene_id"].replace("_gene", ""), dir_dict["bam"] + remove_prefix_and_suffix(file)))

        line["out"]["counts_tpm"] = round(line["out"]["counts_tpm"] / len(line["files"]),0)

        ## Counts for all other files, calculate average
        for file in alternative:
            line["out"]["counts_alt"].append(get_expression_count(feature_counts_tpm, \
                    line["out"]["gene_id"].replace("_gene", ""), \
                    dir_dict["bam"] + remove_prefix_and_suffix(file)))
        
        if line["out"]["counts_alt"] != []:
            line["out"]["counts_alt"] = ";".join(str(x) for x in line["out"]["counts_alt"])
        else:
            line["out"]["counts_alt"] = "N/A"
            
        pos = line["out"]["var_pos"]
        chrom = line["out"]["chromosome"]
        rna_variant = line["out"]["rna_alt"]
        dna_variant = line["out"]["dna_alt"]
        is_snp = line["is_snp"]
        feature_type = ""
        
        if is_snp:
            feature_type = "SNP"
        else:
            feature_type = "NON-SNP"
        
        with lock:
            variants_dict[chrom+":"+str(pos)] = {feature_type: {"rna":string_to_dict(rna_variant, is_snp), "dna":string_to_dict(dna_variant, is_snp)},
                                                    "counts"    :  line["out"]["counts_tpm"]}#.rsplit(';', maxsplit=1)[0]}
        
    ## Add / character if no matching annotation feature is found 
    elif line["out"]["gene_id"] == 'No_Annotation_Found':
        line["out"]["counts_tpm"] = "/"
        line["out"]["counts_alt"] = "/"

    ## Add reads counts to each output line (if available)  
    ## Concatenate output strings (variant information & counts)
    out_line = list(line["out"].values())
    return '\t'.join(map(str, out_line)) + '\n'
    
    

def process_position(position):
            
    global variant_positions
    global gff_lines
    global extracted_gff_lines
    global bam_file_string
    global output_lines
    global agreement_dict
    global gff_parsed
    global feature

    variant = variant_positions[position]
    
    ## Matching variants in both RNA and DNA at position
        
    for hit in variant["hits"]:
        if hit in agreement_dict:
            agreement_dict[hit] += 1
        else:
            agreement_dict[hit] = 1

    extracted_gff_lines, extend = get_gff_lines(gff_parsed, variant["chromosome"], variant["position"], 0, "gene", [])
    gff_lines.append(extracted_gff_lines)

    for line in extracted_gff_lines:
        values = line.split('\t')
        
        extracted_feature = extract_attribute(values[8], 'gene_name=')
        
        ## Filter variants
        #### --> Only consider variants that contain feature for gene_name
        if extracted_feature in feature:
        
            start, end = values[3], values[4]

            variant_line = {
                "out": {
                    "gene_id":   extract_attribute(values[8], 'ID='),
                    "gene_name": extracted_feature, # gene_name for human data, Name for bacteria
                    "chromosome":variant["chromosome"],
                    "start":     start,
                    "end":       end,
                    "var_pos":   variant["position"],
                    "extension": get_extend_character(start, end, variant["position"]),
                    "ref":       variant["ref"],
                    "dna_alt":   calculate_allele_percentages(variant["dna"], 1),
                    "rna_alt":   calculate_allele_percentages(variant["rna"], len(variant["hits"])),
                    "counts_alt":[],
                    "counts_tpm":0,
                    "annotation":variant["annotation"][0]
                },
                "files": variant["hits"],
                "is_snp": check_for_snp(variant["dna"], variant["rna"])
            }
            
            for file in variant["hits"]:
                if file not in bam_file_string:
                    bam_file_string.append(file)
            
            with lock:
                output_lines.append(variant_line)
                    
def reset_global_variables():
    ## Reset gloabl variables
    global manager
    global variant_positions
    global gff_lines
    global extracted_gff_lines
    global bam_file_string
    global output_lines
    global agreement_dict
    global variants_dict
    global feature_counts_tpm
    global rna_prefixes

    global matching_count
    global total_count
    global dna_count
    global rna_count
    global both_dna_and_rna_count
    global dir_dict
    global gff_parsed
    global feature
    
    
    variant_positions = manager.dict()
    gff_lines = []
    extracted_gff_lines = []
    bam_file_string = []
    output_lines = manager.list()
    agreement_dict = manager.dict()
    variants_dict = manager.dict()

    matching_count = 0
    total_count = 0
    dna_count = 0
    rna_count = 0
    both_dna_and_rna_count = 0
    rna_prefixes = []

    gff_parsed = []
    feature_counts_tpm = ""
    dir_dict = {}

    feature = ""