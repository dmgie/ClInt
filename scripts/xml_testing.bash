#!/bin/bash

# XML-Dateipfad
XML_FILE="./arguments.xml"

# Array deklarieren
declare -A programs_array

# XML mit xmlstarlet analysieren und Array erstellen
while IFS= read -r line; do
    # Kategorie-Element überprüfen
    if [[ $line =~ \<category.*type=\"([^\"]+)\" ]]; then
        category_type="${BASH_REMATCH[1]}"
    fi

    # Programm-Element überprüfen
    if [[ $line =~ \<program.*name=\"([^\"]+)\" ]]; then
        program_name="${BASH_REMATCH[1]}"
        if [[ -z ${programs_array["$category_type"]} ]]; then
            programs_array["$category_type"]=$program_name
        else
            programs_array["$category_type"]+=",$program_name"
        fi
    fi
done < <(xmlstarlet sel -t -c "/programs/category" "$XML_FILE")

# Array ausgeben
for key in "${!programs_array[@]}"; do
    echo "Category: $key, Programs: ${programs_array[$key]}"
done
