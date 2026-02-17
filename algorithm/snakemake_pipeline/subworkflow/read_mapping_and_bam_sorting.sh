#!/bin/bash


# Two parameters: configuration file path and BAM file path
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <config_file_path> <bam_file_path>"
    exit 1
fi

CONFIG_FILE="$1"
ZUMIS_BAM="$2"

echo "configuration file path: $CONFIG_FILE"
echo "BAM file path: $BAM_FILE"

# 1. Read the config file
PROJECT_NAME=$(grep '^project:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
Program_folder=$(grep '^pipeline_path:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
algorithem_path="$Program_folder/algorithem/snakemake_pipeline"
CELLTYPE1=$(grep '^celltype1:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
CELLTYPE2=$(grep '^celltype2:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
CELLTYPE3=$(grep '^celltype3:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)

ANNOTATION=$(grep '^annotation:' $CONFIG_FILE | cut -d ':' -f 2 | sed 's/^ *//g' | sed 's/ *$//g')
Data_results_path=$(grep '^out_dir:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
THREADS=$(grep '^num_threads:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
ZUMIS_PATH=$(grep '^zUMIs_directory:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)

# 2.Output directories
OUT_DIR=$Data_results_path
ZUMIS_OUTPUT_DIR="$OUT_DIR/zUMIs_output"
SPLIT_BAM_OUTPUT_DIR="$OUT_DIR/split_bam_output"
READ_SEPARATION_RESULTS_DIR="$OUT_DIR/read_separation_results"

mkdir -p $OUT_DIR $ZUMIS_OUTPUT_DIR $SPLIT_BAM_OUTPUT_DIR $READ_SEPARATION_RESULTS_DIR

# 3. Extract file path
CELLBC_CELLTYPE1="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE1}_subset_cellBClist.txt"
CELLBC_CELLTYPE2="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE2}_subset_cellBClist.txt"
CELLBC_CELLTYPE3="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE3}_subset_cellBClist.txt"

SPLIT_BAM_CELLTYPE1="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE1}_subset.bam"
SPLIT_BAM_CELLTYPE2="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE2}_subset.bam"

# 4. Execute the script to extract the barcode
echo "Extracting barcodes and splitting BAM files..."
BARCODE_SCRIPT="$algorithem_path/SNP_calling/extract_barcodes_and_split_bam.sh"
chmod +x $BARCODE_SCRIPT
echo "BARCODE_SCRIPT: $BARCODE_SCRIPT"
echo "ANNOTATION: $ANNOTATION"
echo "SPLIT_BAM_OUTPUT_DIR: $SPLIT_BAM_OUTPUT_DIR"
bash $BARCODE_SCRIPT $ANNOTATION $SPLIT_BAM_OUTPUT_DIR


# 5. Split the BAM file
echo "Splitting BAM files for cell types..."
samtools view -D BC:$CELLBC_CELLTYPE1 -@ $THREADS -o $SPLIT_BAM_CELLTYPE1 $ZUMIS_BAM
samtools view -D BC:$CELLBC_CELLTYPE2 -@ $THREADS -o $SPLIT_BAM_CELLTYPE2 $ZUMIS_BAM

echo "Pipeline completed successfully."
