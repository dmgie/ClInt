import argparse
from vartable import *

## File handler for vartable.py - COVID dataset
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
    parser.add_argument('--out_dir', required=False, default="../..", help="output directory, DEFAULT='.'")
    
    args = parser.parse_args()
    
    vcf_dna_path = getattr(args, "vcf_dna_path")
    vcf_rna_path = getattr(args, "vcf_rna_path")
    bam_path = getattr(args, "bam_path")
    out_dir = getattr(args, "out_dir")
    gff_path = getattr(args, "gff_path")
    
    ## Read sample preparation file
    #### -> Read in each line of meta file 
    #### -> Create dictionary patient:{rna/dna:{qbic_id, alt_id}}
    #### -> Read dict into filename_prefixes: Per patient dna/rna ids (prefixes of actual filenames)

    metafile_path = "./Q001H_sample_preparations_20230803115337.tsv"
    meta_dict_sorted = get_meta_dict(metafile_path)
    filename_prefixes = get_filname_prefixes(meta_dict_sorted)
    
    # vcf_rna_path="../../test_dir/vcf_annotated_rna/"
    # vcf_dna_path="../../test_dir/vcf_annotated_dna/"
    # bam_path = "../../test_dir/bams/"
    gff_path = "./gencode.v44.chr_patch_hapl_scaff.annotation.gff"
    
    ## Execute vartable
    #### -> Execute per patient
    #### -> Only consider matched patients, i.e. DNA and RNA samples are availbale
    #### -> Give file prefixes of RNA/DNA files to vartable executer as arguments
    #### -> Use information from sample preparation file
    
    variants_dict_complete = {}
    variants_dict_per_patient = {}
    
    
    
    for patient in filename_prefixes:
        print("#####", patient)
        if "DNA" in filename_prefixes[patient] and "RNA" in filename_prefixes[patient]:
            print("\n##############################", patient, "##############################")
            variants_dict_per_patient[patient] = execute_vartable(vcf_dna_path, vcf_rna_path, bam_path, gff_path, filename_prefixes[patient], patient, out_dir)
        # break
    
    # print("vardict from 2", variants_dict_per_patient)

def execute_vartable(vcf_dna_path, vcf_rna_path, bam_path, gff_path, filename_prefixes, patient, output_dir):
    
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
                                 filename_prefixes["RNA"])
    
    return variants_dict  
        
def get_filname_prefixes(meta_dict_sorted):
    filename_prefixes = {}
    
    for patient in meta_dict_sorted:
        for type in meta_dict_sorted[patient]:
            for sample in meta_dict_sorted[patient][type]:
                if patient in filename_prefixes:
                    if type in filename_prefixes[patient]:
                        filename_prefixes[patient][type].append(meta_dict_sorted[patient][type][sample]["alt_identifier"])
                    else:
                        filename_prefixes[patient][type] = [meta_dict_sorted[patient][type][sample]["alt_identifier"]]
                else:
                    filename_prefixes[patient] = {type: [meta_dict_sorted[patient][type][sample]["alt_identifier"]]}
        
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
            
            if "TCR" not in alt_identifier and "BCR" not in alt_identifier:
                if patient_id in meta_dict:
                    if type in meta_dict[patient_id]:
                        sample_number = f"Sample_{len(meta_dict[patient_id][type]) + 1}_{type}"
                        dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier}
                        meta_dict[patient_id][type][sample_number] = dict_entry
                    else:
                        sample_number = f"Sample_1_{type}"
                        dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier}
                        meta_dict[patient_id][type] = {sample_number: dict_entry}
                else:
                    sample_number = f"Sample_1_{type}"
                    dict_entry = {'qbic_identifier': qbic_identifier, 'alt_identifier': alt_identifier}
                    meta_dict[patient_id] = {type: {sample_number: dict_entry}}        
            
    return sort_dict_by_keys(meta_dict)

def natural_sort_key(s):
        def atoi(text):
            return int(text) if text.isdigit() else text
        return [atoi(c) for c in s.split()]

def sort_dict_by_keys(dictionary):
    sorted_items = sorted(dictionary.items(), key=lambda x: natural_sort_key(x[0]))
    sorted_dict = dict(sorted_items)
    return sorted_dict

main()
