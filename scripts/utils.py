import os
import pysam
import subprocess

def read_vcf_file(directories):
    variants_dict = {}

    for type in directories:
        directory = os.fsencode(directories[type])
        for file in os.listdir(directory):
            filename = os.fsencode(file).decode()
            if filename.endswith(".vcf"):
                with pysam.VariantFile(directory+file) as vcf:
                    for record in vcf:
                        chrom = record.chrom
                        pos = record.pos
                        ref = record.ref
                        alt = ",".join(map(str, record.alts))
                        
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
                    
                                
        if not extracted_lines:
            extend += 10
            return get_gff_lines(gff_file, chromosome, position, extend, feature, extracted_lines)
            
    print(f"Result found at extension {extend}")                                                
    return extracted_lines

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
    
    return expression_count

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



## check the expression in files with differente allele
