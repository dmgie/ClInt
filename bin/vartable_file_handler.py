import os
import argparse

if __name__ == "__main__":

    print("Starting VarTable file handler")

    parser = argparse.ArgumentParser(description="File handler fpr VarTable Analysis")
    parser.add_argument('--vcf_path', required=True, help="Top level VCF folder containing sample folders")

    args = parser.parse_args()

    vcf_path = getattr(args, "vcf_path")

    print("Top level .vcf directory path:", vcf_path)

