#!/usr/bin/env bash

## This is a file that contains shared functions throughout the scripts
display_menu() {
    local prompt=$1
    local options_string=$2

    # Options_string is given as a comma delimited list
    # echo $2 >&2

    # Set the Internal Field Separator (IFS) to comma
    IFS='§' read -ra options <<<"$options_string"

    # trim newlines
    prompt=${prompt//$'\n'/}

    # trim whitespace
    prompt=${prompt//$'\t'/}

    local num_options=${#options[@]}

    echo "$prompt" >&2
    for ((i = 0; i < num_options; i++)); do
        descript=$(get_description options[i])
        if [[ -z $descript ]]; then
            echo -e "    $((i + 1)). ${options[i]}" >&2
        else
            echo -e "    $((i + 1)). ${options[i]} - $descript" >&2
        fi
    done

    exec </dev/tty

    local selected_options=()
    local current_selection=""
    local choice

    # NOTE: We can replace the "Selected: " part of the prompt with the actual name (using `options[choice-1]`)
    # but it might fill up the horizontal space a lot if the names are long
    # Though it might be a bit more advnatageous so you don't need to keep looking at the list
    # of programs the entire time
    # NOTE: Maybe remove the deselection if its repeated? Don't know
    # FIXME: For some reason we can for one of them, enter it twice. Fix this somehow
    # FIXME: Comma is leftover after deselecting, fix this somehow.
    while true; do
        read -rp $'[1mEnter your choice [1-'$num_options'][Selected: '"$current_selection"']: [0m' choice
        if [[ -z $choice ]]; then
            break
        # if choice is already selected, deselect
        elif [[ $choice =~ ^[1-$num_options]$ ]] && [[ $current_selection =~ $choice ]]; then
            selected_options=("${selected_options[@]/${options[choice - 1]}/}")
            current_selection=${current_selection//$choice/}
            # if selected_options=(); then
            #     current_selection=""
            # fi
            echo "    deselected: ${options[choice - 1]}" >&2
        elif [[ $choice =~ ^[1-$num_options]$ ]]; then
            selected_options+=("${options[choice - 1]}")

            # If its the first, don't add a comma
            if [[ -z $current_selection ]]; then
                current_selection="$choice"
            else
                current_selection="${current_selection},$choice"
            fi
            echo "    selected: ${options[choice - 1]}" >&2

        # if comma separated list is entered
        elif [[ $choice =~ ^[1-$num_options](,[1-$num_options])*$ ]]; then
            IFS=',' read -ra choices <<<"$choice"
            for choice in "${choices[@]}"; do
                selected_options+=("${options[choice - 1]}")
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
    local default=$((num_cores / 8))

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
    local total_ram_gb=$((total_ram / 1024 / 1024))

    # 1/8th of the RAM would seem like a good amount?
    local default=$((total_ram_gb / 8))

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

_reset_programs() {
    local -n programs=$1
    local -n order_programs=$2
    for program in "${order_programs[@]}"; do
        programs[$program]=0
    done
}
_print_selected() {
    local -n categories=$1
    local -n categories_order=$2
    local -n programs=$3

    echo
    echo -e "\e[4m\e[1mSelected Programs: \e[0m\e[0m"
    for category in "${categories_order[@]}"; do
        echo -e "\e[1m$category: \e[0m"
        ## For each program in the category, check if it is selected, and if so, print it
        category_programs=${categories[$category]}

        # Ignore existing spaces, and split on commas, used for the loop below
        IFS='§' read -r -a category_programs <<<"$category_programs"

        # program_full being the full description, e.g fastp (All-in-one)
        for program_full in "${category_programs[@]}"; do
            # Take only the first word of the program, as that is the actual program name e.g fastp
            # TODO: Make it search / match lower case so to not worry about capitalisation mistakes
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
    # NOTE: For the above tood: maybe not, since sometimes we just want to
    # run a single program again
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
    echo "---------------------------"
    echo -e "\e[1mAre these choices correct? \e[0m"
    _print_selected $1 $2 $3
    select yn in "Yes" "No"; do
        case $yn in
        Yes) break ;;
        No)
            _reset_programs $3 $4
            category_chooser $1 $2 $3 $4
            break
            ;;
        esac
    done
}

## Get program arguments specified in ./arguments.xml config file
## Read arguments using 'xmlstarlet', choose arguments depending on program name
## The xml file contains files or directories as arguments
## Store strings and return
get_arguments() {
    local -n prog=$1

    while IFS= read -r line; do
        output+=("$line")
    done < <(xmlstarlet sel -t \
        --var prog_name="'$prog'" \
        -v '//programs/category[*]/program[@name = $prog_name]/arg' \
        -nl temp_arguments.xml)

    echo "${output[@]}"
}

## Read desription string from xml config file
get_description() {
    local -n prog=$1

    while IFS= read -r line; do
        output+=("$line")
    done < <(xmlstarlet sel -t \
        --var prog_name="'$prog'" \
        -v '//programs/category[*]/program[@name = $prog_name]/description' \
        -nl temp_arguments.xml)

    echo "${output[@]}"
}

## Check whether all necessary input arguments are set
## return 1 if all required arguments are given
## return 0 if argument is missing
arguments_complete() {
    argument_list=("$@")

    for argument in "${argument_list[@]}"; do
        EMPTY="_EMPTY"
        if [[ "$argument" == *"$EMPTY" ]]; then
            echo "$argument, please check input arguments"
            return 1 
        fi
    done

    return 0
}

## Read categories and asociated programs from XML file usinf xmlstarlet
## Store into associated array and return
get_categories() {
    local xml_file="$1"
    declare -A rna_categories

    while IFS= read -r line; do
        if [[ $line == *"<category"* ]]; then
            category_type=$(echo "$line" | awk -F '"' '/category/{print $2}')
        fi

        if [[ $line == *"<program"* ]]; then
            program_name=$(echo "$line" | awk -F '"' '/program/{print $2}')
            if [[ -z ${rna_categories["$category_type"]} ]]; then
                rna_categories["$category_type"]=$program_name
            else
                rna_categories["$category_type"]+="§$program_name"
            fi
        fi
    done < <(xmlstarlet sel -t -c "/programs/category" "$xml_file")

    declare -p rna_categories
}