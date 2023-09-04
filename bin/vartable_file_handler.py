import argparse
from re import T
from subprocess import Popen, PIPE
import glob
import subprocess

## File handler for vartable.py - COVID dataset
#### -> Enables parallel processing of variant samples using Popen
#### -> Takes the top level VCF directory as input, runs vartable.py on all folders below
#### -> Define vartable input arguments in this file

def main():
    ## Argument parsing

    parser = argparse.ArgumentParser(description="File handler for VarTable Analysis")
    parser.add_argument('--vcf_path', required=True, help="Top level VCF folder containing sample folders")
    args = parser.parse_args()
    vcf_path = getattr(args, "vcf_path")
    
    ## Read sample preparation file
    #### -> Read in each line of meta file 
    #### -> Create dictionary patient:{rna/dna:{qbic_id, alt_id}}
    #### -> Read dict into filename_prefixes: Per patient dna/rna ids (prefixes of actual filenames)

    metafile_path = "./Q001H_sample_preparations_20230803115337.tsv"
    meta_dict_sorted = get_meta_dict(metafile_path)
    filename_prefixes = get_filname_prefixes(meta_dict_sorted)
    
    # patient_match_counter = 0
    
    ## Execute vartable
    #### -> Execute per patient
    #### -> Only consider matched patients, i.e. DNA and RNA samples are availbale
    #### -> Give file prefixes of RNA/DNA files to vartable executer as arguments
    #### -> Use information from sample preparation file
    
    for patient in filename_prefixes:
            if "DNA" in filename_prefixes[patient] and "RNA" in filename_prefixes[patient]:
                print("\n##############################", patient, "##############################")
                execute_vartable(vcf_path, filename_prefixes[patient])
                # patient_match_counter += 1
            #break
              
    # print("Number of patients:", len(meta_dict_sorted))
    # print("Number of DNA/RNA matches:", patient_match_counter)
  

def natural_sort_key(s):
        def atoi(text):
            return int(text) if text.isdigit() else text
        return [atoi(c) for c in s.split()]

def sort_dict_by_keys(dictionary):
    sorted_items = sorted(dictionary.items(), key=lambda x: natural_sort_key(x[0]))
    sorted_dict = dict(sorted_items)
    return sorted_dict

def execute_vartable(vcf_path, filename_prefixes):
    # Get all directories below vcf_path
    # directories = glob.glob(f'{vcf_path}/*/')

    # for d in directories: print(d)
    directory = "."

    ## Command to run vartable.py
    #### -> Define vartable arguments
    cmd_list = [["python", "./vartable.py", "--dna", "../../test_dir/vcf_annotated_dna/", \
                                            "--rna", "../../test_dir/vcf_annotated_rna/", \
                                            "--bam", "../../test_dir/bams/", \
                                            "--out", f'{directory}/vartable_output', \
                                            "--gff", "../../../../local_scratch/ClINT/RawData/ref_genome.gff", \
                                            "--gff_filter", "False", \
                                            "--snpEff", "True", \
                                            "--agreement", "True", \
                                            "--dna_startswith", *filename_prefixes["DNA"], \
                                            "--rna_startswith", *filename_prefixes["RNA"]] \
                                                
                                            ]# for directory in directories]

    ## Collect commands, create Popen objects
    # prc_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmd_list]
    prc = Popen(cmd_list[0], universal_newlines=True)
    return_code = prc.wait()
    
    # # Execute each of the processes, enable console output
    # for prc in prc_list:
    #     stdout, stderr = prc.communicate()
    #     print(stdout.decode('ascii'), stderr.decode())
    #     prc.wait()
    # subprocess.run(cmd_list[0], shell=True, capture_output=True, text=True)
        
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

main()
