#!/usr/bin/env bash

## This is a file that contains shared functions throughout the scripts
display_menu() {
    local prompt=$1
    local options_string=$2

    # Set the Internal Field Separator (IFS) to comma
    IFS=',' read -ra options <<< "$options_string"

    local num_options=${#options[@]}

    echo "$prompt" >&2
    for ((i=0; i<num_options; i++)); do
        echo "$((i+1)). ${options[i]}" >&2
    done

    exec < /dev/tty

    local selected_options=()
    local choice

    while true; do
        read -rp $'[1mEnter your choice [1-'$num_options'] (Enter to finish)[0m: ' choice
        if [[ -z $choice ]]; then
            break
        elif [[ $choice =~ ^[1-$num_options]$ ]]; then
            selected_options+=("${options[choice-1]}")
        else
            echo "Invalid choice. Please try again." >&2
        fi
    done

    echo "${selected_options[@]}"
}


# Used as thread/core variable in the workflow
get_core_count() {
    local prompt="Enter the number of cores to use"
    # Check how many cores the system has
    local num_cores=$(grep -c ^processor /proc/cpuinfo)

    # 1/8th of the cores would seem like a good amount?
    local default=$((num_cores/8))

    local core_count

    while true; do
        read -rp "$prompt [Default: $default / Available: $num_cores]: " core_count
        if [[ -z $core_count ]]; then
            core_count=$default
            break
        elif [[ $core_count =~ ^[1-9][0-9]*$ ]]; then
            break
        else
            echo "Invalid core count. Please try again." >&2
        fi
    done

    echo $core_count
}

