import argparse
from subprocess import Popen, PIPE
import glob

## File handler for vartable.py
#### -> Enables parallel processing of variant samples using Popen
#### -> Takes the top level VCF directory as input, runs vartable.py on all folders below
#### -> Define vartable input arguments in this file

if __name__ == "__main__":

    ## Argument parsing

    parser = argparse.ArgumentParser(description="File handler for VarTable Analysis")
    parser.add_argument('--vcf_path', required=True, help="Top level VCF folder containing sample folders")
    args = parser.parse_args()
    vcf_path = getattr(args, "vcf_path")

    ## Get all directories below vcf_path
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
