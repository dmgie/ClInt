from cmath import exp
import json
import os
from collections import defaultdict
from utils import *
from scipy import stats
import numpy as np
from tqdm import tqdm

patient_vector = defaultdict(lambda: defaultdict(dict))
expression_dict = defaultdict(lambda: defaultdict(dict))

print("PARSE GFF")
gff = parse_gff("./gencode.v44lift37.annotation.gff3")
print("DONE")
directory = './vartable_output_grch37'

def get_genotype(input_string):
    ## Get genotype from nucleotide characters

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
    ## Parse all VarTable output files into dict
    #### -> keys are patient name and severity condition

    condition = os.path.basename(file_path).split('_')[0]
    patient = os.path.basename(file_path).split('_')[1]
    
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
            
            variant_data = {'SNP': {'rna': rna_counts, 
                                    'dna': dna_counts,
                                    'rna_nuc': rna_nucleotides,
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
            

print("PARSE VARTABLE")
for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith('.tsv') and 'vartable' in file:
            file_path = os.path.join(root, file)
            parse_vartable(file_path)
print("DONE")


variants_list = get_all_variants()


def calculate_log_fold_change(list1, list2):
    mean1 = np.mean(list1) # has variant
    mean2 = np.mean(list2) # hasnt variant
    
    # add small value if 0 to be able to calculate fold change
    epsilon = 1e-8 
    if mean1 < epsilon:
        mean1 += epsilon
    if mean2 < epsilon:
        mean2 += epsilon
    
    fold_change = mean1 / mean2 
    log_fold_change = np.log2(fold_change)

    return log_fold_change

significant = {}

def ks(variants_liste):
    ## Calculate Kolmogorov-Smirnov test for each variant position
    #### -> Expression values of patients that have variant vs. have no variant

    for variant in tqdm(variants_liste):

        has_variant = []
        hasnt_variant = []

        has_variant_condition = []
        hasnt_variant_condition = []

        has_variant_rna_nucleos = []
        has_variant_dna_nucleos = []

        hasnt_variant_rna_nucleos = []
        hasnt_variant_dna_nucleos = []

        for condition in patient_vector:
            for patient in patient_vector[condition]:                
                    
                ## Patients that have variant
                if variant in patient_vector[condition][patient]:

                    has_variant.append(patient_vector[condition][patient][variant]["counts"])
                    has_variant_condition.append(condition)

                    rna_nucleo = patient_vector[condition][patient][variant]["SNP"]["rna_nuc"]
                    dna_nucleo = patient_vector[condition][patient][variant]["SNP"]["dna_nuc"]

                    has_variant_rna_nucleos.append(rna_nucleo)
                    has_variant_dna_nucleos.append(dna_nucleo)

                ## Patients that do not have variant
                else:
                    
                    gene_name = find_gene_by_variant_position(gff, variant.split(":")[0], int(variant.split(":")[1]))

                    hasnt_variant_condition.append(condition)

                    hasnt_variant_rna_nucleos.append("REF")
                    hasnt_variant_dna_nucleos.append("REF")

                    if gene_name in expression_dict[condition][patient]:

                        hasnt_variant.append(expression_dict[condition][patient][gene_name])
                    else:
                        hasnt_variant.append(0)

        ## Only report variants present in at least 2 patients
        if len(has_variant) >= 2 and len(hasnt_variant) >= 2:
            ks_statistic, p_value = stats.ks_2samp(has_variant, hasnt_variant)

            ## Report all variants, p-value filtering later
            if p_value < 1:

                entry = {
                    "pos": variant,
                    "p_value": p_value,
                    "log_fc": calculate_log_fold_change(has_variant, hasnt_variant),
                    "has_var": has_variant,
                    "hasnt_var": hasnt_variant,
                    "has_var_cond": has_variant_condition,
                    "hasnt_var_cond": hasnt_variant_condition,
                    "has_var_rna_nuc": has_variant_rna_nucleos,
                    "has_var_dna_nuc": has_variant_dna_nucleos,
                    "hasnt_var_rna_nuc": hasnt_variant_rna_nucleos,
                    "hasnt_var_dna_nuc": hasnt_variant_dna_nucleos
                }

                significant[variant] = entry

print("CALCULATE KS VALUES")
ks(variants_list)

with open('output_all_variants_new.json', 'w') as out_file:
    json.dump(significant, out_file)
