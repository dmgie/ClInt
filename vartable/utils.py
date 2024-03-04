## File containing all of the helper functions needed in vartable.py

import os
import pysam
import subprocess
import pandas as pd
from tqdm import tqdm

def search_vcf_position_matches(directories, dna_startswith, rna_startswith, use_snpeff, num_processes):
    """Get DNA/RNA position matches
    """

    variants_dict = {}
    prefix_to_files = {
        "rna":{prefix: [] for prefix in rna_startswith},
        "dna":{prefix: [] for prefix in dna_startswith}
    }
    prefix_dict = {
        "rna": rna_startswith,
        "dna": dna_startswith
    }
    
    for type in directories:
        if type == "dna" or type == "rna":
            if not os.path.exists(directories[type]):
                return None
            directory = os.fsencode(directories[type])
            for root, _, files in os.walk(directory):
                for file in files:
                    file = file.decode()
                    for prefix in prefix_dict[type]:
                        if file.startswith(prefix):
                            prefix_to_files[type][prefix].append(os.path.join(root.decode(), file))
                            if use_snpeff: 
                                print("START SNPEFF")
                                snpeff(directories[type]+file, directories["out"])
                                directory = os.fsencode(directories["out"]+"/vcf_annotated/")
                                file+=".ann.vcf"
                            with pysam.VariantFile(directory+file.encode(), threads=num_processes) as vcf:
                                for variant in vcf:
                                    
                                    genotype = str(variant.samples.values()[0]["GT"])
                                    pos = int(variant.pos)
                                    
                                    chrom = variant.chrom
                                    if "chr" not in chrom:
                                        chrom = "chr" + chrom
                                    
                                    chrom_and_pos = chrom+":"+str(variant.pos)
                                    
                                    ref = variant.ref
                                    alt = ",".join(map(str, variant.alts))
                                    
                                    alleles = ""
                                    
                                    ## Heterozygous variants
                                    if genotype == "(0, 1)":
                                        alleles += ref+"," + alt
                                    
                                    # ## Homozygous variants 
                                    elif genotype == "(1, 1)":
                                        alleles += alt+"," + alt
                                        
                                    alt = alleles
                                    
                                    if use_snpeff: info = variant.info["ANN"]
                                    else: info = "/"
                                
                                    if chrom_and_pos in variants_dict:
                                        if type in variants_dict[chrom_and_pos]:
                                            variants_dict[chrom_and_pos][type] += ","+alt
                                            variants_dict[chrom_and_pos]["hits"].update({file:alt})
                                        else:
                                            variants_dict[chrom_and_pos][type] = alt
                                            variants_dict[chrom_and_pos]["hits"] = {file:alt}
                                    else:
                                        variants_dict[chrom_and_pos] = {
                                        "ref": ref,
                                        type: alt,
                                        "hits": {file:alt},
                                        "chromosome":chrom,
                                        "position":pos,
                                        "annotation":info,
                                        "is_snp":False
                                        } 

    for field in ['rna', 'dna']:
        if any(prefix_to_files.get(field, {}).values()):
            return variants_dict
        else:
            return None

def read_bam_filenames(bam_path, prefix):
    """Get filenames of BAM files corresponding to patient

    Keyword arguments:
    @bam_path  -- BAM directory
    @prefix    -- BAM prefix from cohort metafile
    """
    bam_files = []
    for bam_file in os.listdir(bam_path):
        if bam_file.endswith(".bam") and bam_file.startswith(tuple(prefix)):
            bam_files.append(bam_file)

    return bam_files

def get_gff_header(gff_file):
    """Get GFF header lines

    Keyword arguments:
    @gff_file  -- GFF filename
    """
    header_lines=[]
    with open(gff_file, 'r+') as infile:
        for line in infile:
            if line.startswith('#'):
               header_lines.append(line.strip().split("\t"))
    return header_lines

def parse_gff(gff_file):
    """Parse GFF to dict

    Keyword arguments:
    @gff_file  -- GFF filename
    """

    gff_data = {}
    with open(gff_file, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chromosome, _, feature_type, start, end, *_ = fields
                    start, end = int(start), int(end)
                    gff_data.setdefault(chromosome, []).append((start, end, feature_type, line.strip()))
    return gff_data

def get_gff_lines(gff_data, chromosome, position, extend, feature, extracted_lines, max_extend=0):
    """Extract gene annotation for variant position
    """
    
    for (start, end, feature_type, line) in gff_data.get(chromosome, []):
        if extend <= max_extend:
            if feature_type == feature and start - extend <= position <= end + extend:
                extracted_lines.append(line)
                
    if extracted_lines == []:
        return ["/\t/\t/\t/\t/\t.\t-\t.\tID=No_Annotation_Found;locus_tag=No_Annotation_Found"], max_extend

    ## OPTIONAL: Recursively search for variant matches in EXTEND +/- mode
    ## IMPORTANT: Drastically reduces runtime
           
    # if not extracted_lines:
    #     extend += 10
    #     return get_gff_lines(gff_data, chromosome, position, extend, feature, extracted_lines)
    
    # print(f"Result found at extension {extend}")
    return extracted_lines, extend

def create_gff(gff_file, header_lines, extracted_lines):
    """Create new GFF file of extracted lines
    
    Keyword arguments:
    @gff_file        -- new .gff file name
    @header_lines    -- all header lines from original .gff
    @extracted_lines -- extracted .gff lines
    """
    with open(gff_file, 'w') as f:
        for header_line in header_lines:
            f.write(f"{header_line[0]}\n")
        
        for line in extracted_lines:
            gff_line = "\t".join(line)
            f.write(f"{gff_line}\n")

    return gff_file

def extract_attribute(input_string, query):
    """Extracts strings from GFF attribute column (9)

    Keyword arguments:
    @input_string -- attribute string
    @query        -- query term, i.e. "Name="
    """
    start = input_string.find(query)
    if start == -1:
        return None

    start += len(query)
    end = input_string.find(";", start)
    if end == -1:
        return None

    return input_string[start:end]

def run_featurecounts(input_bam, annotation_gff, output_counts_file, feature, num_processes):
    """Run featureCounts as shell subprocess
    """
    cmd = f"featureCounts -p -a {annotation_gff} -o {output_counts_file} -t {feature} -g {'gene_id'} {input_bam} -T {num_processes}" #-p gene_id for human data, locus_tag for bacteria

    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(e)

def get_expression_count(fc_file, gene_id, filename):
    """Extract expression counts for certain genes from featureCounts file
    
    Keyword arguments:
    @fc_file      -- featureCounts filename
    @gene_id      -- ID of gene to be extracted
    @filename     -- BAM filename
    """

    expression_count = None
    with open(fc_file, "r") as file:
        lines = file.readlines()
        header = lines[1].strip().split("\t")
        index_geneid = header.index("Geneid")
        index_filename = None
        
        for index, field in enumerate(header):
            if filename in field:
                index_filename = index
                
        if index_filename == None:
            return "/"
        
        for line in lines[1:]:
            data = line.strip().split("\t")
            if data[index_geneid] == gene_id:
                expression_count = int(data[index_filename])
                break
    
    if expression_count is not None:
        return expression_count
    else:
        print("Gene ID not found", gene_id)
        return "/"

def remove_prefix_and_suffix(filename, prefix = "haplotype_", suffix = ".ann.vcf"):
    """Removes substrings to trim fileendinding from .vcf files

    Keyword arguments:
    @filename  -- input string
    @prefix    -- string prefix, default: "haplotype_"
    @sufffix   -- string sufffix, default: ".vcf"
    """
    
    filename = filename.replace(prefix, "")
    filename = filename.replace(suffix, "")
    filename = filename.replace(".vcf", "")
    
    return filename

def remove_after_bam(input_string):
    """Removes substrings to trim fileendinding from .bam files

    Keyword arguments:
    @input_string  -- BAM input file
    """

    index = input_string.find(".bam")
    if index != -1:
        return input_string[:index + 4]
    return input_string

def string_to_bool(input):
    """Turn string to boolean variable
    """
    if input == "True": return True
    else: return False

def calculate_allele_percentages(allele_sequence, len_hits):
    """Calculate allele percentages in format [A,C,G,T]

    Keyword arguments:
    @allele_sequence  -- base sequence of variant
    @len_hits         -- variant character length
    """
    allele_sequence = allele_sequence.split(",")
    allele_count = {}
    total_alleles = len(allele_sequence)

    for allele in allele_sequence:
        allele_count[allele] = allele_count.get(allele, 0) + 1

    allele_percentages = {}
    allele_percentages["TOT"] = len_hits
    
    for allele, count in reversed(dict(sorted(allele_count.items(), key=lambda item: item[1])).items()):
        percentage = round(count / total_alleles, 2)
        allele_percentages[allele] = percentage

    counts_string = []
    for count in allele_percentages:
        counts_string.append(f"{count}:{allele_percentages[count]}")

    return " ".join(counts_string)


def get_extend_character(start, end, position):
    """Get annotation character for EXTEND +/- mode

    Keyword arguments:
    @start     -- feature start position
    @end       -- feature end position
    @position  -- variant position
    """

    extend_character = "" 
    
    if start == "/":
        return "/"
    
    if int(start) <= position <= int(end):
        extend_character = f"."
    elif int(start) > position:
        extend_character = f"-"
    else:
        extend_character = f"+"
        
    return extend_character

def check_for_snp(split_nucleotides_dna, split_nucleotides_rna):
    """Classify variant as SNP or non-SNP based on variant length

    Keyword arguments:
    @split_nucleotides_dna -- DNA variant nucleotides
    @split_nucleotides_rna -- RNA variant nucleotides
    """
    for nucleotide in split_nucleotides_dna.split(",") + split_nucleotides_rna.split(","):
        if len(nucleotide) > 1:
            return False
    return True

def filter_filenames_by_prefixes(filenames_dict, prefixes):
    """Filter filenames for patient file prefixes

    Keyword arguments:
    @filenames_dict -- patient filenames
    @prefixes       -- file prefixes
    """
    
    filtered_filenames = []
    
    for prefix in prefixes:
        if not any(key.startswith(prefix) for key in filenames_dict.keys()):
            filtered_filenames.append(prefix)
        
    return filtered_filenames

def calculate_tpm(feature_counts_output, output_dir):
    """TPM normalization of featureCounts file

    Keyword arguments:
    @feature_counts_output -- featureCounts .txt file
    @output_dir            -- directory for TPM file
    """

    counts_df = pd.read_csv(feature_counts_output, sep='\t', skiprows=1, index_col=0)

    for col in counts_df.columns[5:]:
        counts_df[col] = counts_df[col] / counts_df[counts_df.columns[4]]
    
    for col in counts_df.columns[5:]:
        counts_df[col] = (counts_df[col] * 1e6 / counts_df[col].sum()).round(0).astype(int)

    with open(f'{output_dir}/fc_output_tpm.txt', 'w') as f:
        with open(feature_counts_output, 'r') as input_file:
            first_line = input_file.readline().strip()
            f.write(first_line + '\n')

        counts_df.to_csv(f, sep='\t', index=True)

def create_agreement(agreement, dna_count, rna_count, matching_count, output, sample_name, condition):
    """Create agreement file of DNA/RNA and RNA/DNA agreement rates

    Keyword arguments:
    @agreement      -- agreement rate dict
    @dna_count      -- total DNA variant count
    @rna_count      -- total RNA variant count
    @matching_count -- Concording DNA/RNA count
    @output         -- output dir
    @sample_name    -- Patient Name
    @condition      -- Severity condition
    """
    with open(f"{output}/{condition}_{sample_name}_agreement.tsv", 'w', encoding='utf-8') as file:
        header = "####Agreement rate: proportion of dna variants also found in rna per rna .vcf file "
        header += f"RNA_COUNT: {rna_count} DNA_COUNT: {dna_count} MATCHING: {matching_count}\n"
        header += "####vcf_filename\tagreement_rate_dna\tagreement_rate_rna\n"
        file.write(header)

        for filename in agreement:
            agreement_dna = round(agreement[filename]/dna_count, 2)
            agreement_rna = round(agreement[filename]/rna_count, 2)
            file.write(filename + "\t" + str(agreement_dna) + "\t" + str(agreement_rna) + "\n")

def snpeff(path, out):
    """snpEff annotation of VCF files

    Keyword arguments:
    @path   -- input VCF file
    @out    -- output directory
    """
    print(f"Start: SNPEff analysis for {path}.")
    file_name = os.path.basename(path)
    os.makedirs(out+"/vcf_annotated", exist_ok=True)
    snpeff_command = f"java -Xmx8g -jar ../../snpEff/snpEff.jar -noLog -noStats GRCh37.75 {path} > {out}/vcf_annotated/{file_name}.ann.vcf"
    subprocess.run(snpeff_command, shell=True, capture_output=True, text=True)


def string_to_dict(input_string, is_snp):
    """Turn SNP variant nucleotide string into allele dict

    Keyword arguments:
    @input_string   -- input nucleotide string
    @is_snp         -- snp_classification
    """
    
    if is_snp == False:
        return input_string
    
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    parts = input_string.split()
    
    for part in parts:
        key, value = part.split(':')
        if key.upper() in nucleotide_counts:
            nucleotide_counts[key.upper()] = float(value)
            
    return nucleotide_counts

def find_gene_by_variant_position(gff_data, chromosome, variant_position):
    """Return gene name by variant position and chromosome name from GFF

    Keyword arguments:
    @gff_data         -- parsed GFF file
    @chromosome       -- chromosome name
    @variant_position -- variant position on chromosome
    """
    if chromosome in gff_data:
        entries = gff_data[chromosome]
        for start, end, feature_type, line in entries:
            if start <= variant_position <= end and feature_type == 'gene':
                gene_name_start = line.find('gene_name=') + len('gene_name=')
                gene_name_end = line.find(';', gene_name_start)
                gene_name = line[gene_name_start:gene_name_end]
                return gene_name
    return None
    