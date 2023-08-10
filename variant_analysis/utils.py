import os
import pysam
import subprocess

def search_vcf_position_matches(directories):
    variants_dict = {}

    for type in directories:
        if type == "dna" or type == "rna":
            directory = os.fsencode(directories[type])
            for file in os.listdir(directory):
                filename = os.fsencode(file).decode()
                if filename.endswith(".vcf") and filename.startswith("haplotype_dedup_snc_trimmed") or filename.startswith("S_aureus"):
                    # print("found", filename)
                    with pysam.VariantFile(directory+file) as vcf:
                        for record in vcf:
                            chrom = record.chrom
                            pos = record.pos
                            ref = record.ref
                            alt = ",".join(map(str, record.alts))
                        
                            if pos in variants_dict:
                                if type in variants_dict[pos]:
                                    variants_dict[pos][type] += ","+alt
                                    variants_dict[pos]["hits"].update({filename:alt})
                                else:
                                    variants_dict[pos][type] = alt
                                    variants_dict[pos]["hits"] = {filename:alt}
                            else:
                                variants_dict[pos] = {
                                "ref": ref,
                                type: alt,
                                "hits": {filename:alt},
                                "chromosome":chrom
                                }         
    return variants_dict

def read_bam_filenames(bam_path):


    bam_files = []
    for bam_file in os.listdir(bam_path):
        if bam_file.endswith(".bam") and bam_file.startswith("dedup_snc_trimmed"):
            bam_files.append(bam_file)

    return bam_files

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

    print("Works #########")

    cmd = f"featureCounts -a {annotation_gtf} -o {output_counts_file} -t {feature} -g {'locus_tag'} {input_bam}"

    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(e)

def get_expression_count(fc_file, gene_id, filename):
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

def remove_prefix_and_suffix(filename, prefix = "haplotype_", suffix = ".vcf"):
    """Removes substrings to trim fileendinding from .vcf files

    Keyword arguments:
    @filename  -- input string
    @prefix    -- string prefix, default: "haplotype_"
    @sufffix   -- string sufffix, default: ".vcf"
    """
    
    filename = filename.replace(prefix, "")
    filename = filename.replace(suffix, "")
    return filename

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
    
    if start is "/":
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


def normalize_read_count(featurecounts_file):

    gene_counts = {}  # Dictionary zur Speicherung der Genzählungen

    with open(featurecounts_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                gene_id = parts[0]
                counts = [int(x) for x in parts[6:]]  # Zählungen für jede Probe
                gene_counts[gene_id] = counts

# Berechne die Summe der Counts pro Gen
    gene_sums = {gene_id: sum(counts) for gene_id, counts in gene_counts.items()}

# Lese die Gesamtzahl der gelesenen Reads aus der FeatureCounts-Datei
    total_reads = None
    with open(featurecounts_file, 'r') as f:
        for line in f:
            if line.startswith("Status"):
                total_reads = int(line.split(":")[1].strip())
                break

    if total_reads is None:
        raise ValueError("Total reads count not found in the FeatureCounts file.")

# Berechne die TPM-Werte für jedes Gen und speichere sie in einer Ausgabedatei
    tpm_output_file = "tpm_normalized_output.txt"
    with open(tpm_output_file, 'w') as output:
        output.write("GeneID\tTPM\n")
        for gene_id, counts in gene_counts.items():
            tpm = [(count / gene_sums[gene_id]) * 1_000_000 / (total_reads / 1_000_000) for count in counts]
            output.write(f"{gene_id}\t{' '.join(map(str, tpm))}\n")

    print(f"TPM-normalisierte Werte wurden in '{tpm_output_file}' gespeichert.")


## TODO: check the expression in files with differente allele