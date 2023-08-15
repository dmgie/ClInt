process CSV_CREATION {
    publishDir "${baseDir}", mode: 'copy', overwrite: true
    // Create a CSV file detailing the sample_name and the file paths to the read(s)
    // This is depending on whether the "paired" flag was enabled or disabled when launching the command

    input: 
        path read_dir
    output:
        path "*.csv"
        stdout emit: std

    shell:
    r1 = "_R1" 
    r2 = "_R2"
    depth = 1

    shell:
    if (params.paired) {
        '''#!/usr/bin/env bash
        # Get the input arguments
        input_path=$(realpath !{read_dir})
        input_files=($(find $input_path -name "*.f*q.*"))
        r1_pattern=!{r1}
        r2_pattern=!{r2}

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
            echo "${sample_list[$i]},${r1_files[$i]},${r2_files[$i]}" >> clint_metadata.csv
        done

        cat clint_metadata.csv
        '''
    } else 
        '''
        input_path=$(realpath !{read_dir})
        input_files=($(find $input_path -maxdepth !{depth} -type f -regex '.*f*q\\(.gz\\|.bz2\\)?'))
        files_path=($(find $input_path -maxdepth !{depth} -type f -regex '.*f*q\\(.gz\\|.bz2\\)?'))

        echo ${input_files[@]}

        for i in ${!input_files[@]}; do
            sample_name[$i]=${input_files[$i]##*/}
            sample_name[$i]=${sample_name[$i]%%.*}
        done
        echo "sample_name,r1_path,r2_path" > clint_metadata.csv # Header
        for i in ${!sample_name[@]}; do
            echo "${sample_name[$i]}, ${files_path[$i]}" >> clint_metadata.csv
        done

        cat clint_metadata.csv
        '''
}
