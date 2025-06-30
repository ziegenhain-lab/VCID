#!/bin/bash

# 检查是否传入了配置文件路径
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

# 从命令行参数中获取配置文件路径
CONFIG_FILE="$1"

# 读取配置文件中的变量
PROJECT_NAME=$(grep '^project:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
Program_folder=$(grep '^pipeline_path:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
algorithem_path="$Program_folder/algorithem/snakemake_pipeline"
CELLTYPE1=$(grep '^celltype1:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
CELLTYPE2=$(grep '^celltype2:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
CELLTYPE3=$(grep '^celltype3:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
ANNOTATION=$(grep '^annotation:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
CELLBCs_FOR_zUMIs=$(grep '^zUMIs_barcode_list:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
CELLBCs_FOR_READSEPARATION=$(grep '^readseparation_barcode_list:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
Differetial_SNPs_from_other_source=$(grep '^Differetial_SNPs_from_other_source:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
BAM_file_from_other_source=$(grep '^BAM_file_from_other_source:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
DATA_RESULTS_PATH=$(grep '^out_dir:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)
THREADS=$(grep '^num_threads:' $CONFIG_FILE | cut -d ':' -f 2 | xargs)

# 创建输出目录
READ_SEPARATION_RESULTS_DIR="$DATA_RESULTS_PATH/read_separation_results"
SPLIT_CHUNK_OUTPUT="$DATA_RESULTS_PATH/read_separation_results/split_chunk_output"
CHUNK_COMBINED_OUTPUT="$DATA_RESULTS_PATH/read_separation_results/chunk_combined_output"
mkdir -p $READ_SEPARATION_RESULTS_DIR
mkdir -p $SPLIT_CHUNK_OUTPUT
mkdir -p $CHUNK_COMBINED_OUTPUT

echo "$DATA_RESULTS_PATH"

# 2. 根据是否提供了外部 SNP 文件路径来选择 VCF 文件路径
if [ -z "$Differetial_SNPs_from_other_source" ] || [ "$Differetial_SNPs_from_other_source" = "~" ]; then
    # 如果没有提供外部 SNP 文件路径，使用默认路径
    PARENT_DIR=$(dirname "$DATA_RESULTS_PATH")
    VCF_FILE="$PARENT_DIR/SNP_calling_output/${PROJECT_NAME}.processed.vcf.gz"
    echo "Default VCF file detected: $VCF_FILE"
else
    # 否则，使用提供的外部 SNP 文件路径
    VCF_FILE="$Differetial_SNPs_from_other_source"
    VCF_FILE_INDEX="${Differetial_SNPs_from_other_source}.tbi"
    echo "External SNP file detected: $Differetial_SNPs_from_other_source"
    if [ -f "$VCF_FILE_INDEX" ]; then
        echo "VCF index already exists: $VCF_FILE_INDEX"
    fi

fi

# 3. 根据是否提供了外部 BAM 文件路径来选择 BAM 文件路径
if [ -z "$BAM_file_from_other_source" ] || [ "$BAM_file_from_other_source" = "~" ]; then
    # 如果没有提供外部 BAM 文件路径，使用默认路径
    BAM_FILE="$DATA_RESULTS_PATH/${PROJECT_NAME}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam"
    echo "Default BAM file detected: $BAM_FILE"
else
    # 否则，使用提供的外部 BAM 文件路径
    echo "External BAM file detected"
    BAM_FILE="$BAM_file_from_other_source"
    echo "External BAM file locates in: $BAM_file_from_other_source"
fi


# 1. 执行 Doublet Separation 脚本
echo "Running Doublet Separation..."
RUN_BAM_FILE=""
Test_BAM_FILE="$DATA_RESULTS_PATH/${PROJECT_NAME}.filtered.Aligned.GeneTagged.UBcorrected.sorted.test.bam"

# 如果 readseparation list 和 zUMIs list 不同，则进行 BAM 文件筛选，筛选出我们想要的数据进行小批量的测试
if [ "$CELLBCs_FOR_READSEPARATION" != "$CELLBCs_FOR_zUMIs" ]; then
    echo "The test list has been detected and will be used to filter the BAM file for read separation."
    samtools view -D BC:$CELLBCs_FOR_READSEPARATION -@ $THREADS -o $Test_BAM_FILE $BAM_FILE
    RUN_BAM_FILE=$Test_BAM_FILE
else
    echo "The test list is the same as the zUMIs list, so the original BAM file will be used for read separation."
    RUN_BAM_FILE=$BAM_FILE
fi

# 检查RUN_BAM_FILE的bai索引文件是否存在
if [ ! -f "${RUN_BAM_FILE}.bai" ]; then
    echo "Index file for BAM file not found. Creating index..."
    samtools index $RUN_BAM_FILE
fi

# First, split the chunks of the original bam file & vcf file, and then use the split bam files for doublet separation
SPLIT_CHUNK_AND_READ_SEPARATION_SCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/split_chunk.sh"
READ_SEPARATION_RSCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/DoubletsSeparation_SNPbased_v2.0_separation.R"
READ_SEPARATION_COMBINING_SCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/read_separation_result_combining_and_counting.sh"
MAKE_COUNT_TABLE_SCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/make_read_counts.R"
SPLIT_CHUNK_OUTPUT="$DATA_RESULTS_PATH/SNP_calling_output"
chmod +x $SPLIT_CHUNK_AND_READ_SEPARATION_SCRIPT
chmod +x $READ_SEPARATION_COMBINING_SCRIPT
bash $SPLIT_CHUNK_AND_READ_SEPARATION_SCRIPT $RUN_BAM_FILE $VCF_FILE $SPLIT_CHUNK_OUTPUT $READ_SEPARATION_RSCRIPT
bash $READ_SEPARATION_COMBINING_SCRIPT $CONFIG_FILE $CHUNK_COMBINED_OUTPUT $SPLIT_CHUNK_OUTPUT $MAKE_COUNT_TABLE_SCRIPT

# 2. 执行 Doublet Separation Plotting 脚本
echo "Running Doublet Separation Plotting..."
echo "Maybe someday I will modify this script to use the new version of plotting script.... Please waite for the update"
#PLOTTING_SCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/DoubletsSeparation_SNPbased_v0.6.0_plotting.R"
#Rscript $PLOTTING_SCRIPT --input "$CHUNK_COMBINED_OUTPUT/all_out_read.txt" --annotation $ANNOTATION \
#    --outdir $READ_SEPARATION_RESULTS_DIR --yaml $CONFIG_FILE

# 3. 执行 BAM 文件分离操作
echo "Splitting BAM files for cell types..."

# 目前取消这些parallel extraction，因为split和parallel并不能提升效率，限速步骤在解压缩bam
#SEPARATE_BAM_SCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/parallel_extract_bam.sh"
#chmod +x $SEPARATE_BAM_SCRIPT
# 分离 CellType1 和 CellType2 的 BAM 文件
#$SEPARATE_BAM_SCRIPT "$CHUNK_COMBINED_OUTPUT/${CELLTYPE1}_readID.txt" $RUN_BAM_FILE $CHUNK_COMBINED_OUTPUT 30
#$SEPARATE_BAM_SCRIPT "$CHUNK_COMBINED_OUTPUT/${CELLTYPE2}_readID.txt" $RUN_BAM_FILE $CHUNK_COMBINED_OUTPUT 30
samtools view -@ 30 -N "$CHUNK_COMBINED_OUTPUT/${CELLTYPE1}_readID.txt" -b $RUN_BAM_FILE -o "$CHUNK_COMBINED_OUTPUT/${CELLTYPE1}_separated_bam.txt"
samtools view -@ 30 -N "$CHUNK_COMBINED_OUTPUT/${CELLTYPE2}_readID.txt" -b $RUN_BAM_FILE -o "$CHUNK_COMBINED_OUTPUT/${CELLTYPE2}_separated_bam.txt"

echo "Counting reads number for each cell before read separation..."
READCOUNT_SCRIPT="$algorithem_path/DoubletSeparation_SNPbased_v1.0/read_per_cell_from_bam.sh"
READCOUNT_CSV="$CHUNK_COMBINED_OUTPUT/read_count_per_cell_basedon_bam.csv"
bash $READCOUNT_SCRIPT $RUN_BAM_FILE $READCOUNT_CSV 30

echo "Pipeline completed successfully."
