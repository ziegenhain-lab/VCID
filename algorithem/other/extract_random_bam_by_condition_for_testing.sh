#!/bin/bash

# 参数检查
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <metadata.txt> <input.bam> <output_dir>"
  exit 1
fi

metadata_file="$1"
input_bam="$2"
output_dir="$3"

# 检查文件是否存在
if [[ ! -f "$metadata_file" ]]; then
  echo "Error: Metadata file '$metadata_file' not found."
  exit 1
fi

if [[ ! -f "$input_bam" ]]; then
  echo "Error: BAM file '$input_bam' not found."
  exit 1
fi

mkdir -p "$output_dir"

# 获取 sort_population 和 XC_DNBPE 的列号
cond_col=$(head -n1 "$metadata_file" | tr '\t' '\n' | grep -n "^sort_population$" | cut -d: -f1)
bc_col=$(head -n1 "$metadata_file" | tr '\t' '\n' | grep -n "^XC_DNBPE$" | cut -d: -f1)

if [[ -z "$cond_col" || -z "$bc_col" ]]; then
  echo "Error: Could not find 'sort_population' or 'XC_DNBPE' columns in the metadata file."
  exit 1
fi

# 获取所有唯一的 sort_population
sort_populations=$(awk -v c="$cond_col" 'BEGIN{FS="\t"} NR>1 {print $c}' "$metadata_file" | sort -u)

for cond in $sort_populations; do
  echo "Processing sort_population: $cond"

  # 先生成随机6位码
  randcode=$(tr -dc 'A-Za-z0-9' </dev/urandom | head -c6)

  # 生成barcode子集文件名
  barcode_file="${output_dir}/${cond}_${randcode}_subset_cellBClist.txt"
  # 输出bam路径
  out_bam="${output_dir}/${cond}_${randcode}.bam"

  # 获取该 sort_population 下所有 XC_DNBPE，随机抽10个写入文件
  awk -v FS="\t" -v c="$cond_col" -v b="$bc_col" -v val="$cond" 'NR>1 && $c==val {print $b}' "$metadata_file" | sort -u | shuf | head -n 10 > "$barcode_file"

  barcode_count=$(wc -l < "$barcode_file")
  if [[ "$barcode_count" -eq 0 ]]; then
    echo "Warning: No barcodes found for sort_population $cond. Skipping."
    continue
  fi

  echo "  Selected $barcode_count barcodes saved to $barcode_file"
  
  samtools view -D BC:$barcode_file -@ 40 -o $out_bam $input_bam

  echo "  Output BAM written to: $out_bam"
done
