#!/bin/bash

# Usage: ./count_bc_from_bam.sh input.bam output.csv 8

input_bam="$1"
output_csv="$2"
threads="$3"

# Step 1: Extract reads from BAM with multi-threaded decompression
# Step 2: Extract BC tag
# Step 3: Count using awk
# Step 4: Output CSV

samtools view -@ "$threads" "$input_bam" \
  | awk '
    {
      bc = ""; 
      for (i=12; i<=NF; i++) {
        if ($i ~ /^BC:Z:/) {
          split($i, a, ":");
          bc = a[3];
          break;
        }
      }
      if (bc != "") {
        count[bc]++;
      }
    }
    END {
      print "Barcode,Read_Count";
      for (b in count) {
        print b "," count[b];
      }
    }
  ' > "$output_csv"

echo "✅ Barcode read count written to: $output_csv"
