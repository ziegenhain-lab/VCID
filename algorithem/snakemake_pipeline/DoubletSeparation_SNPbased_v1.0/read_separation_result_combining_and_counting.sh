#!/bin/bash
#! R脚本路径需要指明！

# 接收输入参数
YAML_FILE=$1             # YAML文件路径
OUT_DIR=$2               # 输出目录
CHUNK_OUTPUT_DIR=$3     # 输入目录，里面是各个chunk的输出结果
MAKE_READ_COUNTS=$4 # 生成计数文件的脚本路径

# 读取YAML文件中的CELLTYPE1和CELLTYPE2（假设你已经有了一个方法来从yaml读取这些信息）
CELLTYPE1=$(grep "^celltype1:" "$YAML_FILE" | cut -d ':' -f2 | xargs)
CELLTYPE2=$(grep "^celltype2:" "$YAML_FILE" | cut -d ':' -f2 | xargs)

# 创建输出目录
mkdir -p "$OUT_DIR"

# 合并所有的out_read_readID.txt，并根据CELLTYPE1和CELLTYPE2筛选出对应的readID
echo "合并 out_read_readID.txt 文件并筛选 readID..."

# 使用GNU Parallel并行化处理每个子文件夹的out_read_readID.txt
find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read_readID.txt" | parallel "
  cat {} >> $OUT_DIR/all_out_read_readID.txt
"

# 根据CELLTYPE1筛选出对应的readID并保存到文本文件
grep -w "$CELLTYPE1" "$OUT_DIR/all_out_read_readID.txt" | cut -f1 > "$OUT_DIR/${CELLTYPE1}_readID.txt"
grep -w "$CELLTYPE2" "$OUT_DIR/all_out_read_readID.txt" | cut -f1 > "$OUT_DIR/${CELLTYPE2}_readID.txt"

echo "合并并筛选 readID 完成."

# 合并所有的out_read.txt文件
echo "合并 out_read.txt 文件..."

# 设置最终输出文件路径
FINAL_OUT="$OUT_DIR/all_out_read.txt"
# 找到第一个文件，并提取其 header
FIRST_FILE=$(find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read.txt" | head -n 1)
# 写入 header 到最终输出文件
head -n 1 "$FIRST_FILE" > "$FINAL_OUT"

# 并行处理所有文件（去掉 header），并追加到最终输出文件
find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read.txt" | grep -v "$FIRST_FILE" | \
  parallel "tail -n +2 {}" >> "$FINAL_OUT"

cp "$FINAL_OUT" "$OUT_DIR/all_out_read_backup.txt"



# 找到第一个文件
#FIRST_FILE=$(find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read.txt" | head -n 1)

# 1. 把 header（以RG开头的行）提取并写入最终输出文件
#grep '^RG' "$FIRST_FILE" > "$OUT_DIR/all_out_read.txt"

# 2. 并行处理所有文件（包括第一个），去掉 header，然后追加到输出
#find "$CHUNK_OUTPUT_DIR" -type f -name "*_out_read.txt" | \
#  parallel "tail -n +2 {}" >> "$OUT_DIR/all_out_read.txt"



# 提供基因数量统计
echo "对合并后的 out_read.txt 进行基因数量统计..."

# 调用 R 脚本处理计数数据
Rscript $MAKE_READ_COUNTS "$FINAL_OUT" "$OUT_DIR"


# 输出基因数量统计结果
echo "基因数量统计完成，结果保存在 $OUT_DIR/separated_counts.rds"

# 统计输出完成
echo "所有任务完成！"

