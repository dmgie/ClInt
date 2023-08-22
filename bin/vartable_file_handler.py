import argparse
from subprocess import Popen, PIPE
import glob

## File handler for vartable.py
#### -> Enables parallel processing of variant samples using Popen
#### -> Takes the top level VCF directory as input, runs vartable.py on all folders below
#### -> Define vartable input arguments in this file

def main():
    ## Argument parsing

    parser = argparse.ArgumentParser(description="File handler for VarTable Analysis")
    parser.add_argument('--vcf_path', required=True, help="Top level VCF folder containing sample folders")
    args = parser.parse_args()
    vcf_path = getattr(args, "vcf_path")

    # metafile_path = "./Q001H_sample_preparations_20230803115337.tsv"
    # meta_file_sorted = get_meta_file(metafile_path)

    # for patient in meta_file_sorted:
    #     print("Patient:", patient)
    #     for sample in meta_file_sorted[patient]:
    #         print(sample, meta_file_sorted[patient][sample])
    #     print("\n")

    # print("Number of patients:", len(meta_file_sorted))

    execute_vartable(vcf_path)


def natural_sort_key(s):
        def atoi(text):
            return int(text) if text.isdigit() else text
        return [atoi(c) for c in s.split()]

def sort_dict_by_keys(dictionary):
    sorted_items = sorted(dictionary.items(), key=lambda x: natural_sort_key(x[0]))
    sorted_dict = dict(sorted_items)
    return sorted_dict

def get_meta_file(metafile_path):
    meta_file = {}

    with open(metafile_path, "r") as metafile:
        lines = metafile.readlines()

        for line in lines[1:]:
            line_content = line.strip().split('\t')

            patient_id = line_content[4]
            qbic_identifier = line_content[0]
            alt_identifier = line_content[1]

            if "BCR" not in alt_identifier and "TCR"  not in alt_identifier:
                if patient_id in meta_file:
                    ## Patient # already in dict, extend
                    sample_number = f"Sample_{len(meta_file[patient_id])+1}"
                    dict_entry = {sample_number:{'qbic_identifier': qbic_identifier,'alt_identifier': alt_identifier}}
                    meta_file[patient_id].update(dict_entry)

                else:
                    ## New dict entry
                    meta_file[patient_id] = {"Sample_1":{
                        'qbic_identifier': qbic_identifier,
                        'alt_identifier': alt_identifier
                    }}
            
    return sort_dict_by_keys(meta_file)

def execute_vartable(vcf_path):
    # Get all directories below vcf_path
    directories = glob.glob(f'{vcf_path}/*/')

    for d in directories: print(d)

    ## Command to run vartable.py
    #### -> Define vartable arguments
    cmd_list = [["python", "./vartable.py", "--vcf", directory, \
                                            "--bam", "../../../../local_scratch/ClINT/working_files/deduped_bams/", \
                                            "--out", f'{directory}/vartable_output', \
                                            "--gff", "../../../../local_scratch/ClINT/RawData/ref_genome.gff", \
                                            "--gff_filter", "False", \
                                            "--snpEff", "True", \
                                            "--agreement", "True"] \
                                            for directory in directories]

    ## Collect commands, create Popen objects
    prc_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmd_list]

    ## Execute each of the processes, enable console output
    for prc in prc_list:
        stdout, stderr = prc.communicate()
        print(stdout.decode('ascii'), stderr.decode())
        prc.wait()

main()
