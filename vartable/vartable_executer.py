import argparse
import time
from vartable import *
from vectorization import *
from collections import OrderedDict

## Executer for vartable.py - COVID patient cohort dataset
#### -> Enables parallel processing of variant samples using Popen
#### -> Takes the top level VCF directory as input, runs vartable.py on all folders below
#### -> Define vartable input arguments in this file

def main():

    ## Argument parsing
    parser = argparse.ArgumentParser(description="File handler for VarTable Analysis")
    parser.add_argument('--vcf_dna_path', required=False, help="Top level VCF folder containing dna samples")
    parser.add_argument('--vcf_rna_path', required=False, help="Top level VCF folder containing rna samples")
    parser.add_argument('--bam_path', required=False, help="Top level VCF folder containing bams")
    parser.add_argument('--gff_path', required=False, help="Top level VCF folder containing bams")
    parser.add_argument('--out_dir', required=False, default=".", help="output directory, DEFAULT='.'")
    parser.add_argument('--feature_file_path', required=False, default="./features.txt", help="Path to feature file: containing features to search variants for")
    
    args = parser.parse_args()
    
    vcf_dna_path = getattr(args, "vcf_dna_path")
    vcf_rna_path = getattr(args, "vcf_rna_path")
    bam_path = getattr(args, "bam_path")
    out_dir = getattr(args, "out_dir")
    gff_path = getattr(args, "gff_path")
    feature_file_path = getattr(args, "feature_file_path")


    ## Preset file paths
    vcf_rna_path="../../OUTPUT/vcf/"
    vcf_dna_path="../combined_folder_data/vcf_dna/"
    bam_path = "../../OUTPUT/bam/"
    gff_path = "./gencode.v44lift37.annotation.gff3"
    
    ## Read sample preparation file
    #### -> Read in each line of meta file 
    #### -> Create dictionary patient:{rna/dna:{qbic_id, alt_id}}
    #### -> Read dict into filename_prefixes: Per patient dna/rna ids (prefixes of actual filenames)

    metafile_path = "./Q001H_sample_preparations_20230803115337.tsv"
    meta_dict_sorted = get_meta_dict(metafile_path)

    ## Read in feature file
    #### -> List of (differentially expressed) genes
    #### -> VarTable only considers variants within these genes

    features = parse_feature_file(feature_file_path)
    
    ## Execute vartable
    #### -> Execute per patient
    #### -> Only consider matched patients, i.e. DNA and RNA samples are availbale
    #### -> Give file prefixes of RNA/DNA files to vartable executer as arguments
    #### -> Use information from sample preparation file
    
    variants_dict_per_condition = {}
    variants_dict_per_patient = {}
    variants_list = []
    patient_vectors = {}
    runtimes = {}

    for idx, patient in enumerate(meta_dict_sorted):
        for condition in meta_dict_sorted[patient]:

            if "DNA" in meta_dict_sorted[patient][condition] and "RNA" in meta_dict_sorted[patient][condition]:
                if meta_dict_sorted[patient][condition]["DNA"] != {} and meta_dict_sorted[patient][condition]["RNA"] != {}:  
                
                    start_time = time.time()
                    print("\n############################", patient, condition, "##########################")  
                    
                    filename_prefixes = get_filname_prefixes(meta_dict_sorted, condition)
                    
                    reset_global_variables()
                    vartable_results = execute_vartable(vcf_dna_path, vcf_rna_path, bam_path, gff_path, filename_prefixes[patient], patient, out_dir, features, condition)
                    
                    ## OPTIONAL: Save vartable output in single, binary output file
                    ## Was not used since parsing of individual VarTable files showed to be more effective
                    ## Still might be useful when integrated into a pipeline

                    # if condition not in variants_dict_per_patient:
                    #     variants_dict_per_patient[condition] = {}

                    # if patient not in variants_dict_per_patient[condition]:
                    #     variants_dict_per_patient[condition][patient] = vartable_results
                    # else:
                    #     variants_dict_per_patient[condition][patient].update(vartable_results)
                    
                    
                    # for condition2 in variants_dict_per_patient:
                    #     for patient2 in variants_dict_per_patient[condition2]:

                    #         if condition2 not in variants_dict_per_condition:
                    #             variants_dict_per_condition[condition2] = []
                            
                    #         for variant_position2 in variants_dict_per_patient[condition2][patient2].keys():
                    #             if variant_position2 not in variants_dict_per_condition[condition2]:
                    #                 variants_dict_per_condition[condition2].append(variant_position2)
            
                    end_time = time.time()
                    runtimes[patient+"_"+condition] = end_time-start_time
                    print("RUNTIME:", end_time-start_time)     


    ## OPTIONAL: Save vartable output in single, binary output file
    ## Was not used since parsing of individual VarTable files showed to be more effective
    ## Still might be useful when integrated into a pipeline

    # for variants in variants_dict_per_condition.values():
    #     for variant in variants:
    #         if variant not in variants_list:
    #             variants_list.append(variant)

    # patient_vectors = create_patient_vectors(variants_dict_per_condition, variants_dict_per_patient)
    
    # save_patient_vectors(patient_vectors, out_dir)
    # save_variant_vectors(variants_list, out_dir)
    # save_runtime(runtimes, out_dir)


def execute_vartable(vcf_dna_path, vcf_rna_path, bam_path, gff_path, filename_prefixes, patient, output_dir, features, condition):
    
    ## Command to run vartable.py
    variants_dict = run_vartable(vcf_dna_path, \
                                 vcf_rna_path, \
                                 bam_path, \
                                 f'{output_dir}/vartable_output', \
                                 gff_path, \
                                 patient, \
                                 False, \
                                 False, \
                                 True, \
                                 filename_prefixes["DNA"], \
                                 filename_prefixes["RNA"], \
                                 64, \
                                 features, \
                                 condition)
    
    return variants_dict  
        
def get_filname_prefixes(meta_dict_sorted, condition_type):
    filename_prefixes = {}
    
    for patient in meta_dict_sorted:
        for condition in meta_dict_sorted[patient]:
            if condition == condition_type:
                for type in meta_dict_sorted[patient][condition]:
                    for sample in meta_dict_sorted[patient][condition][type]:
                        if patient in filename_prefixes:
                            if type in filename_prefixes[patient]:
                                filename_prefixes[patient][type].append(meta_dict_sorted[patient][condition][type][sample]["alt_identifier"])
                            else:
                                filename_prefixes[patient][type] = [meta_dict_sorted[patient][condition][type][sample]["alt_identifier"]]
                        else:
                            filename_prefixes[patient] = {type: [meta_dict_sorted[patient][condition][type][sample]["alt_identifier"]]}
        
    return filename_prefixes

def get_meta_dict(metafile_path):
    meta_dict = {}

    with open(metafile_path, "r") as metafile:
        lines = metafile.readlines()

        for line in lines[1:]:
            line_content = line.strip().split('\t')

            patient_id = line_content[4]
            qbic_identifier = line_content[0]
            alt_identifier = line_content[2]
            type = line_content[3]
            condition = line_content[15]
            
            if "TCR" not in alt_identifier and "BCR" not in alt_identifier:
                if patient_id in meta_dict:
                    if condition in meta_dict[patient_id]:
                        if type in meta_dict[patient_id][condition]:
                            sample_number = f"Sample_{len(meta_dict[patient_id][condition][type]) + 1}_{type}"
                            dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier, 'condition': condition}
                            meta_dict[patient_id][condition][type][sample_number] = dict_entry
                        else:
                            sample_number = f"Sample_1_{type}"
                            dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier, 'condition': condition}
                            meta_dict[patient_id][condition][type] = {sample_number: dict_entry}
                    else:
                        sample_number = f"Sample_1_{type}"
                        dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier, 'condition': condition}
                        meta_dict[patient_id][condition] = {type: {sample_number: dict_entry}}
                else:
                    sample_number = f"Sample_1_{type}"
                    dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier, 'condition': condition}
                    meta_dict[patient_id] = {condition: {type: {sample_number: dict_entry}}} 

                dna_entry = {}
                for condition in meta_dict[patient_id]:
                    if "DNA" in meta_dict[patient_id][condition]:
                        dna_entry = meta_dict[patient_id][condition]["DNA"]   

                for condition in meta_dict[patient_id]:
                    if "DNA" not in meta_dict[patient_id][condition]:
                        meta_dict[patient_id][condition]["DNA"] = dna_entry
                
            
    return sort_dict_by_keys(meta_dict)

def natural_sort_key(s):
    def atoi(text):
        return int(text) if text.isdigit() else text
    return [atoi(c) for c in s.split()]

def sort_dict_by_keys(dictionary):
    sorted_items = sorted(dictionary.items(), key=lambda x: natural_sort_key(x[0]))
    sorted_dict = OrderedDict(sorted_items)
    return sorted_dict

def parse_feature_file(path):
    line_list = []

    with open(path, 'r') as file:
        for line in file:
            line = line.strip()
            line_list.append(line)

    return line_list

main()
