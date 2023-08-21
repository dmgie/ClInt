import cmd
import os
import argparse
import subprocess
from subprocess import Popen, PIPE
import glob

if __name__ == "__main__":

    print("Starting VarTable file handler")

    parser = argparse.ArgumentParser(description="File handler fpr VarTable Analysis")
    parser.add_argument('--vcf_path', required=True, help="Top level VCF folder containing sample folders")
    args = parser.parse_args()
    vcf_path = getattr(args, "vcf_path")

    print("Top level .vcf directory path:", vcf_path)

    directories = glob.glob(f'{vcf_path}/*/')

    for d in directories: print(d)

    cmd_list = [["python", "./vartable.py", "--vcf", directory, "--out", f'{directory}/vartable_output', "--snpEff", "True"] for directory in directories]
    prc_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmd_list]


    for prc in prc_list:
        stdout, stderr = prc.communicate()
        print(stdout.decode('ascii'), stderr.decode())
        prc.wait()
