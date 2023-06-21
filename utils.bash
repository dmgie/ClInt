#!/usr/bin/env bash

## This is a file that contains shared functions throughout the scripts
display_menu() {
    local prompt=$1
    local options_string=$2

    # Set the Internal Field Separator (IFS) to comma
    IFS=',' read -ra options <<< "$options_string"

    # trim newlines
    prompt=${prompt//$'\n'/}

    # trim whitespace
    prompt=${prompt//$'\t'/}

    local num_options=${#options[@]}

    echo "$prompt" >&2
    for ((i=0; i<num_options; i++)); do
        echo "    $((i+1)). ${options[i]}" >&2
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
        # if comma separated list is entered
        elif [[ $choice =~ ^[1-$num_options](,[1-$num_options])*$ ]]; then
            IFS=',' read -ra choices <<< "$choice"
            for choice in "${choices[@]}"; do
                selected_options+=("${options[choice-1]}")
            done
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

get_available_ram() {
    local prompt="Enter the amount of RAM to use (in GB)"
    # Check how much RAM the system has
    local total_ram=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    local total_ram_gb=$((total_ram/1024/1024))

    # 1/8th of the RAM would seem like a good amount?
    local default=$((total_ram_gb/8))

    local ram_gb

    while true; do
        read -rp "$prompt [Default: $default / Available: $total_ram_gb]: " ram_gb
        if [[ -z $ram_gb ]]; then
            ram_gb=$default
            break
        elif [[ $ram_gb =~ ^[1-9][0-9]*$ ]]; then
            break
        else
            echo "Invalid RAM amount. Please try again." >&2
        fi
    done

    echo $ram_gb

}

_reset_programs(){
    local -n programs=$1
    local -n order_programs=$2
    for program in "${order_programs[@]}"; do
        programs[$program]=0
    done
}
_print_selected(){
    local -n categories=$1
    local -n categories_order=$2
    local -n programs=$3

    echo -e "\e[4m\e[1mSelected Programs: \e[0m\e[0m"
    for category in "${categories_order[@]}"; do
        echo -e "\e[1m$category: \e[0m"
        ## For each program in the category, check if it is selected, and if so, print it
        category_programs=${categories[$category]}

        # Ignore existing spaces, and split on commas, used for the loop below
        IFS=',' read -r -a category_programs <<< "$category_programs"

        # program_full being the full description, e.g fastp (All-in-one)
        for program_full in "${category_programs[@]}"; do
            # Take only the first word of the program, as that is the actual program name e.g fastp
            program=$(echo $program_full | cut -d' ' -f1)
            if [ "${programs[$program]}" -eq 1 ]; then
                echo -e "       $program_full"
            fi
        done
    done
}

## Display a menu for each category, and allow the user to select which programs they want to run
## If multiple are selected, they are returned as space delimited
## Switch the names in programs to 1 if they are selected
category_chooser() {
    # TODO: Maybe make sure at least one option is selected?
    local -n categories=$1
    local -n categories_order=$2
    local -n programs=$3
    local -n order_programs=$4

    for category in "${categories_order[@]}"; do
        echo #newline
        echo -e "\e[1m$category: \e[0m"
        selected_options=$(display_menu "    Select the programs you want to run:" "${categories[$category]}")
        for program in "${order_programs[@]}"; do
            if [[ " ${selected_options[@]} " =~ " $program " ]]; then
                programs[$program]=1
            fi
        done
    done
    # Ask if the choices were correct
    echo -e "\e[1mAre these choices correct? \e[0m"
    _print_selected $1 $2 $3
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) break;;
            No ) _reset_programs $3 $4; category_chooser $1 $2 $3 $4; break;;
        esac
    done
}
