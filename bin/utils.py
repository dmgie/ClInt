import os
import pysam
import subprocess
import pandas as pd

def search_vcf_position_matches(directories, dna_startswith, rna_startswith, use_snpeff):
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
            directory = os.fsencode(directories[type])
            for root, _, files in os.walk(directory):
                for file in files:
                    file = file.decode()
                    print("###### FILE", file)
                    for prefix in prefix_dict[type]:
                        if file.startswith(prefix):
                            prefix_to_files[type][prefix].append(os.path.join(root.decode(), file))
                            if use_snpeff: 
                                snpeff(directories[type]+file, directories["out"])
                                directory = os.fsencode(directories["out"]+"/vcf_annotated/")
                                file+=".ann.vcf"
                            with pysam.VariantFile(directory+file.encode()) as vcf:
                                for variant in vcf:
                                    chrom = variant.chrom
                                    pos = variant.pos
                                    ref = variant.ref
                                    alt = ",".join(map(str, variant.alts))
                                    if use_snpeff: info = variant.info["ANN"]
                                    else: info = "/"
                                
                                    if pos in variants_dict:
                                        if type in variants_dict[pos]:
                                            variants_dict[pos][type] += ","+alt
                                            variants_dict[pos]["hits"].update({file:alt})
                                        else:
                                            variants_dict[pos][type] = alt
                                            variants_dict[pos]["hits"] = {file:alt}
                                    else:
                                        variants_dict[pos] = {
                                        "ref": ref,
                                        type: alt,
                                        "hits": {file:alt},
                                        "chromosome":chrom,
                                        "annotation":info
                                        }   
    
    return variants_dict

def read_bam_filenames(bam_path, prefix):


    bam_files = []
    for bam_file in os.listdir(bam_path):
        # bam_file = os.fsencode(bam_file).decode()
        if bam_file.endswith(".bam") and bam_file.startswith(tuple(prefix)):
            bam_files.append(bam_file)

    return bam_files

def get_gff_header(gff_file):
    header_lines=[]
    with open(gff_file, 'r+') as infile:
        for line in infile:
            if line.startswith('#'):
               header_lines.append(line.strip().split("\t"))
    return header_lines

def get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines, max_extend=100):
    """Search for occurences of variant position in genomic annotation
       -> Extract .gff lines from original
       -> 1st iteration: look for direct matches
       -> n'th iteration. recursively search for matches in extended search
    """
    with open(gff_file, 'r+') as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 5: 
                    chr_name, source, feature_type, start, end, *_ = fields
                    if extend < max_extend:
                        if feature_type == feature:
                            if chr_name == str(chromosome): 
                                start_pos, end_pos = int(start), int(end)
                                if (start_pos - extend) <= position <= (end_pos + extend):
                                    extracted_lines.append(line.strip())
                    else:
                        return ["/\t/\t/\t/\t/\t.\t-\t.\tID=No_Annotation_Found;locus_tag=No_Annotation_Found"], max_extend
                                           
        if not extracted_lines:
            extend += 10
            return get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines)
            
    print(f"Result found at extension {extend}")                                                
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

def run_featurecounts(input_bam, annotation_gtf, output_counts_file, feature):
    """Run featureCounts as shell subprocess
    """
    cmd = f"featureCounts -a {annotation_gtf} -o {output_counts_file} -t {feature} -g {'locus_tag'} {input_bam} -T 10"

    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(e)

def get_expression_count(fc_file, gene_id, filename):
    print("Expression Count File:", fc_file, filename)
    expression_count = None
    with open(fc_file, "r") as file:
        lines = file.readlines()
        header = lines[1].strip().split("\t")
        index_geneid = header.index("Geneid")
        index_filename = header.index(filename)
        
        for line in lines[1:]:
            data = line.strip().split("\t")
            if data[index_geneid] == gene_id:
                expression_count = int(data[index_filename])
                break
    
    if expression_count is not None:
        return expression_count
    else:
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
    index = input_string.find(".bam")
    if index != -1:
        return input_string[:index + 4]
    return input_string

def string_to_bool(input):
    if input == "True": return True
    else: return False

def calculate_allele_percentages(allele_sequence):
    allele_sequence = allele_sequence.replace(",", "")
    allele_count = {}
    total_alleles = len(allele_sequence)

    for allele in allele_sequence:
        allele_count[allele] = allele_count.get(allele, 0) + 1

    allele_percentages = {}
    allele_percentages["TOT"] = total_alleles
    
    for allele, count in reversed(dict(sorted(allele_count.items(), key=lambda item: item[1])).items()):
        percentage = (count / total_alleles)
        allele_percentages[allele] = percentage

    counts_string = []
    for count in allele_percentages:
        counts_string.append(f"{count}:{allele_percentages[count]}")

    return " ".join(counts_string)


def get_extend_character(start, end, position):
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

def get_alt_files(bam_file_string, bam_path):
    bam_files = os.listdir(bam_path)
        
    bam_files_filitered = []
    bam_files_total = []
        
    for file in list(bam_file_string):
        file = remove_prefix_and_suffix(file)
        bam_files_filitered.append(file)
            
    for file in bam_files:
        if file.startswith("dedup_snc_trimmed"):
            file = remove_prefix_and_suffix(file)
            bam_files_total.append(file)
        
    alt_files = set(bam_files_total) - set(bam_files_filitered)
    return list(alt_files)

def calculate_tpm(feature_counts_output, output_dir):

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

def create_agreement(agreement, dna_count, rna_count, matching_count, output, sample_name):
    with open(f"{output}/{sample_name}_agreement.tsv", 'w', encoding='utf-8') as file:
        header = "####Agreement rate: proportion of dna variants also found in rna per rna .vcf file "
        header += f"RNA_COUNT: {rna_count} DNA_COUNT: {dna_count} MATCHING: {matching_count}\n"
        header += "####vcf_filename\tagreement_rate\n"
        file.write(header)

        for filename in agreement:
            agreement[filename] = round(agreement[filename]/dna_count, 2)
            file.write(filename + "\t" + str(agreement[filename]) + "\n")

def snpeff(path, out):
    print(f"Start: SNPEff analysis for {path}.")
    file_name = os.path.basename(path)
    os.makedirs(out+"/vcf_annotated", exist_ok=True)
    snpeff_command = f"java -Xmx8g -jar ../../snpEff/snpEff.jar -noLog -noStats GRCh37.75 {path} > {out}/vcf_annotated/{file_name}.ann.vcf"
    subprocess.run(snpeff_command, shell=True, capture_output=True, text=True)
