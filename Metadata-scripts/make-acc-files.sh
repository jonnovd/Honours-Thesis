#!/bin/bash

# Ensure the script receives exactly one argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_csv_file>"
  exit 1
fi

# Input CSV file
input_file=$1

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "File not found: $input_file"
  exit 1
fi

# Read the CSV file and process it
tail -n +2 "$input_file" | while IFS=',' read -r i infant run; do
  # If the file doesn't exist, create it and add the header
  if [ ! -f "../srr/${infant}_accessions.txt" ]; then
    echo "acc" > "../srr/${infant}_accessions.txt"  # Add the header
  fi
  # Append the SRR accession to the Infant's file
  echo "$run" >> "../srr/${infant}_accessions.txt"
done

echo "SRR accession files created for each infant with 'acc' as the header."