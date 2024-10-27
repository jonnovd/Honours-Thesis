#!/bin/bash

# Help function
help() {
    echo "Usage: $0 <accession_file>"
    echo "Download SRA files from NCBI using fasterq-dump."
    echo
    echo "  accession_file: CSV file containing the list of accession numbers."
    echo
    echo "Example: $0 subsample.csv"
    exit 1
}

while getopts ":h" opt; do
    case ${opt} in
        h )
            help
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            help
            ;;
    esac
done
shift $((OPTIND -1))

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Error: Incorrect number of arguments."
    help
fi

# Assign the input file to a variable
accession_file="$1"

# Check if the file exists
if [ ! -f "$accession_file" ]; then
    echo "Error: File '$accession_file' not found."
    exit 1
fi

# Loop through the accessions in the CSV file and download using fasterq-dump
first_line=true
while IFS=, read -r accession; do
    # Skipping header lines
    if $first_line; then
        first_line=false
        continue
    fi

    echo "Downloading $accession..."
    fasterq-dump --split-files "$accession"
done < "$accession_file"

echo "Download complete."