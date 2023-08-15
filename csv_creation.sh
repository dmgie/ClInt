#!/usr/bin/env bash
# set -x
# The input to the script will be a list of *.fastq files 
# They will be in the format of <sample_name>_<r1_pattern>.fastq.gz and <sample_name>_<r2_pattern>.fastq.gz

# The output will be a single csv file with the following columns:
# 1. Sample name 
# 2. read1_path
# 3. read2_path


# Get the input arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            shift
            # If given a folder, recursively find all fastq files
            if [[ -d $1 ]]; then
                input_path=$(realpath $1)
                input_files=($(find $input_path -name "*.f*q.*"))
                shift
                # Otherwise add all the files given as input (usually a regex i.e *.fastq.gz)
            else
                while [[ $# -gt 0 && ! $1 =~ ^-.* ]]; do
                    input_files+=("$1")
                    shift
                done
            fi
            ;;
        -1|--r1)
            r1_pattern="$2"
            echo "R1 Pattern: $r1_pattern"
            shift 2
            ;;
        -2|--r2)
            r2_pattern="$2"
            shift 2
            ;;
        *)
            if [[ ! $1 =~ ^-.* ]]; then
                echo "Unknown input (check again): $1"
                exit 1
            fi
            shift
            ;;
    esac
done

# Remove non-paired files i.e files which do not have both *r1* and *r2* patterns
for file in ${input_files[@]}; do
    if [[ $file =~ $r1_pattern ]]; then
        # check if the corresponding r2 file exists
        r2_file=${file/$r1_pattern/$r2_pattern}
        if [[ -f $r2_file ]]; then
            sample_name=${file/$r1_pattern/}
            sample_name=${sample_name%_}
            sample_list+=($sample_name)
            r1_files+=($file)
            r2_files+=($r2_file)
        fi
    fi
done

# Remove paths and extension from the sample names
for i in ${!sample_list[@]}; do
    sample_list[$i]=${sample_list[$i]##*/}
    sample_list[$i]=${sample_list[$i]%%.*}
done



# Create the csv file with the sample names and paths
echo "sample_name,r1_path,r2_path" > clint_metadata.csv # Header
for i in ${!sample_list[@]}; do
    echo "Adding sample ${sample_list[$i]}"
    echo "${sample_list[$i]},${r1_files[$i]},${r2_files[$i]}" >> clint_metadata.csv
done



# echo "  R1 pattern: $r1_pattern"
# echo "  R2 pattern: $r2_pattern"
# echo "  Input files: ${input_files[@]}"
# echo "  Sample list: ${sample_list[@]}"
# echo "  R1 files: ${r1_files[@]}"
# echo "  R2 files: ${r2_files[@]}"

