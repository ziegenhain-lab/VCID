#!/bin/bash

# Input parameters
YAML_FILE=$1             # YAML config file path
OUT_DIR=$2               # Output directory
CHUNK_OUTPUT_DIR=$3     # input directory containing chunked output folders
MAKE_READ_COUNTS=$4 # Script path to generate count file


CELLTYPE1=$(grep "^celltype1:" "$YAML_FILE" | cut -d ':' -f2 | xargs)
CELLTYPE2=$(grep "^celltype2:" "$YAML_FILE" | cut -d ':' -f2 | xargs)

mkdir -p "$OUT_DIR"

# Merge all out_read_readID.txt files and filter out the corresponding readIDs 
# based on CELLTYPE1 and CELLTYPE2.
echo "Merge the out_read_readID.txt file and filter by readID..."

# Use GNU Parallel to process the out_read_readID.txt file for each subfolder.
find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read_readID.txt" | parallel "
  cat {} >> $OUT_DIR/all_out_read_readID.txt
"

# Filter out the corresponding readID based on CELLTYPE1 and save it to a text file
grep -w "$CELLTYPE1" "$OUT_DIR/all_out_read_readID.txt" | cut -f1 > "$OUT_DIR/${CELLTYPE1}_readID.txt"
grep -w "$CELLTYPE2" "$OUT_DIR/all_out_read_readID.txt" | cut -f1 > "$OUT_DIR/${CELLTYPE2}_readID.txt"

echo "Merging and filtering readIDs complete."
echo "Merge out_read.txt file..."

FINAL_OUT="$OUT_DIR/all_out_read.txt"
# Locate the first file and extract its header.
FIRST_FILE=$(find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read.txt" | head -n 1)
# Write header to the final output file
head -n 1 "$FIRST_FILE" > "$FINAL_OUT"

# Process all files in parallel (remove headers) and append them to the final output file.
find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read.txt" | grep -v "$FIRST_FILE" | \
  parallel "tail -n +2 {}" >> "$FINAL_OUT"

cp "$FINAL_OUT" "$OUT_DIR/all_out_read_backup.txt"

echo "Gene counts were performed on the merged out_read.txt file..."

Rscript $MAKE_READ_COUNTS "$FINAL_OUT" "$OUT_DIR"


echo "Gene counts have been completed, and the results are saved in [location]. $OUT_DIR/separated_counts.rds"

echo "All tasks completed.！"

