#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 barcode_annotation.txt Output_Directory"
    exit 1
fi

barcode_file="$1"
output_dir="$2"

mkdir -p "$output_dir"

# Automatically detects column numbers in Sort_Population or sort_population (ignoring case).
sort_col=$(awk -F'\t' 'NR==1 {
    for(i=1; i<=NF; i++) {
        colname = tolower($i)
        if(colname == "sort_population") {
            print i
            exit
        }
    }
}' "$barcode_file")

echo "Sort_Population column number: $sort_col"

if [ -z "$sort_col" ]; then
    echo "Error: Sort_Population column not found"
    exit 1
fi

# Get all unique Sort Population values
sort_values=$(awk -F'\t' -v col="$sort_col" 'NR>1 {print $col}' "$barcode_file" | sort -u)

# Get the column number where XC_DNBPE is located.
col_num=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="XC_DNBPE") print i}' "$barcode_file")

# Ensure the XC_DNBPE column is found.
if [ -z "$col_num" ]; then
    echo "Error: XC_DNBPE column not found"
    exit 1
fi

# Iterate through each Sort Population value and perform data processing.
for value in $sort_values; do

    subset_file="$output_dir/${value}_subset.txt"
    output_file="$output_dir/${value}_subset_cellBClist.txt"

    # Extract the corresponding rows and store them in subset_file
    awk -F'\t' -v val="$value" -v col="$sort_col" '$col == val' "$barcode_file" > "$subset_file"

    # Extract the XC_DNBPE column and save it to the output file.
    awk -F'\t' -v col="$col_num" '{print $col}' "$subset_file" > "$output_file"

    echo "Generated: $output_file"
done
