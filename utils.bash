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
