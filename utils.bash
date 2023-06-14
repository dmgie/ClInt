#!/usr/bin/env bash

## This is a file that contains shared functions throughout the scripts
#
# display_menu() {
#     # Create a display menu for the user to select from, and return the selected option (echo'd)
#     local prompt=$1
#     shift
#     local options=("$@")
#     local num_options=${#options[@]}
#
#     echo "$prompt" >&2
#     for ((i=0; i<num_options; i++)); do
#         echo "$((i+1)). ${options[i]}" >&2
#     done
#
#     exec < /dev/tty
#
#     read -rp $'Enter your choice [1-'$num_options']: ' choice
#
#     if [[ $choice =~ ^[1-$num_options]$ ]]; then
#         selected_option=${options[choice-1]}
#         echo "$selected_option"
#     else
#         echo "Invalid choice. Please try again." >&2
#         display_menu "$prompt" "${options[@]}"
#     fi
# }
#


display_menu() {
    local prompt=$1
    shift
    local options=("$@")
    local num_options=${#options[@]}

    echo "$prompt" >&2
    for ((i=0; i<num_options; i++)); do
        echo "$((i+1)). ${options[i]}" >&2
    done

    exec < /dev/tty

    local selected_options=()
    local choice

    while true; do
        read -rp $'Enter your choice [1-'$num_options'] (Enter to finish): ' choice
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
