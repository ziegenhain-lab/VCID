#!/bin/bash


# 检查是否传入了两个参数：配置文件路径 和 BAM 文件路径
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <config_file_path> <bam_file_path>"
    exit 1
fi

# 从参数中读取配置文件路径和 BAM 文件路径
CONFIG_FILE="$1"
ZUMIS_BAM="$2"

echo "配置文件路径: $CONFIG_FILE"
echo "BAM 文件路径: $BAM_FILE"

# 1. 读取 config 文件
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

# 2. 创建输出目录
OUT_DIR=$Data_results_path
ZUMIS_OUTPUT_DIR="$OUT_DIR/zUMIs_output"
SPLIT_BAM_OUTPUT_DIR="$OUT_DIR/split_bam_output"
READ_SEPARATION_RESULTS_DIR="$OUT_DIR/read_separation_results"

mkdir -p $OUT_DIR $ZUMIS_OUTPUT_DIR $SPLIT_BAM_OUTPUT_DIR $READ_SEPARATION_RESULTS_DIR

# 3. 提取文件路径
CELLBC_CELLTYPE1="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE1}_subset_cellBClist.txt"
CELLBC_CELLTYPE2="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE2}_subset_cellBClist.txt"
CELLBC_CELLTYPE3="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE3}_subset_cellBClist.txt"

SPLIT_BAM_CELLTYPE1="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE1}_subset.bam"
SPLIT_BAM_CELLTYPE2="$SPLIT_BAM_OUTPUT_DIR/${CELLTYPE2}_subset.bam"

# 4. 执行提取条形码的脚本
echo "Extracting barcodes and splitting BAM files..."
BARCODE_SCRIPT="$algorithem_path/SNP_calling/extract_barcodes_and_split_bam.sh"
chmod +x $BARCODE_SCRIPT
echo "BARCODE_SCRIPT: $BARCODE_SCRIPT"
echo "ANNOTATION: $ANNOTATION"
echo "SPLIT_BAM_OUTPUT_DIR: $SPLIT_BAM_OUTPUT_DIR"
bash $BARCODE_SCRIPT $ANNOTATION $SPLIT_BAM_OUTPUT_DIR


# 现在不需要在此处运行 zUMIs了，该功能已经移动到了 snakemake pipeline 中
# 5. 运行 zUMIs
#zUMIs_script="$ZUMIS_PATH/zUMIs.sh"
#bash $zUMIs_script -c -y $CONFIG_FILE 

# 6. zUMIs BAM 索引
#echo "Indexing zUMIs BAM..."
#samtools index $ZUMIS_BAM

# 7. 分割 BAM 文件
echo "Splitting BAM files for cell types..."
samtools view -D BC:$CELLBC_CELLTYPE1 -@ $THREADS -o $SPLIT_BAM_CELLTYPE1 $ZUMIS_BAM
samtools view -D BC:$CELLBC_CELLTYPE2 -@ $THREADS -o $SPLIT_BAM_CELLTYPE2 $ZUMIS_BAM

echo "Pipeline completed successfully."
