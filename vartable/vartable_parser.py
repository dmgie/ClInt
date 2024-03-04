import os
from collections import defaultdict
from utils import *

patient_vector = defaultdict(lambda: defaultdict(dict))
expression_dict = defaultdict(lambda: defaultdict(dict))

print("VARTABLE PARSER - PARSE GFF")
# gff = parse_gff("./gencode.v44lift37.annotation.gff3")
gff = parse_gff("./ClInt/vartable/gencode.v44.chr_patch_hapl_scaff.annotation.gff3")
print("DONE")

# Passe den Pfad zu deinem Verzeichnis an, das die Patientenordner enthält
# directory = '../../vartable_output_grch37'
directory = './vartable_output'

def get_genotype(input_string):
    result_list = [0, 0, 0, 0]
    split_string = input_string.split()

    for part in split_string:
        key_value = part.split(':')
        if len(key_value) == 2:
            key, value = key_value
            if key == 'A':
                result_list[0] = float(value)
            elif key == 'C':
                result_list[1] = float(value)
            elif key == 'G':
                result_list[2] = float(value)
            elif key == 'T':
                result_list[3] = float(value)

    return result_list

def parse_vartable(file_path):
    condition = os.path.basename(file_path).split('_')[0]  # Extrahiere den Condition-Namen aus dem Dateinamen
    patient = os.path.basename(file_path).split('_')[1]  # Extrahiere den Patientennamen aus dem Dateinamen
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('Gene_ID'):
                continue
            data = line.strip().split('\t')
            gene = data[1]
            chromosome = data[2]
            start = data[5]
            variant_key = f'{chromosome}:{start}'

            counts = float(data[-2])
            
            dna_alt = data[8].split(':')[1]
            rna_alt = data[9].split(':')[1]

            dna_nucleotides = get_genotype(data[8])
            rna_nucleotides = get_genotype(data[9])

            
            dna_counts = {x.split(':')[0]: float(x.split(':')[1]) for x in dna_alt.split() if len(x.split(':')) == 2}
            rna_counts = {x.split(':')[0]: float(x.split(':')[1]) for x in rna_alt.split() if len(x.split(':')) == 2}
            
            variant_data = {'SNP': {'rna_nuc': rna_nucleotides,
                                    'dna_nuc': dna_nucleotides
                                    }, 
                            'counts': counts
                            }
            
            patient_vector[condition][patient][variant_key] = variant_data
            expression_dict[condition][patient][gene] = counts

def get_all_variants():

    variant_positions = set()

    def parse_vartable(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('Gene_ID'):
                    continue
                data = line.strip().split('\t')
                chromosome = data[2]
                start = data[5]
                variant_positions.add(f'{chromosome}:{start}')

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.tsv') and 'vartable' in file:
                file_path = os.path.join(root, file)
                parse_vartable(file_path)

    positions_list = list(variant_positions)
    return positions_list
            

def return_variants_list():
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.tsv') and 'vartable' in file:
                file_path = os.path.join(root, file)
                parse_vartable(file_path)

    variants_list = get_all_variants()

    return patient_vector, variants_list











