#!/bin/bash
# 使用方法:
#   ./run_chunks.sh input.bam input.vcf /path/to/output_dir
#
# 要求:
#   1. 已安装 samtools、bcftools、GNU parallel 和 tabix（原VCF需建立索引）。

#
# 脚本说明：
#   1. 使用 samtools idxstats 获取所有染色体名称和长度，然后保存到输出路径下的 chrom_lengths.txt；
#   2. 对每条染色体按每 5,000,000 个碱基切分，若最后一段不足 5,000,000，则取实际长度；
#   3. 对 VCF 文件提取时，右侧区域扩展 2,000bp，但文件以及文件夹名称仍使用原始区间标识（例如 “1-500000”）；(由于目前这一阶段并不需要考虑read pair，所以2000bp完全足够了)
#   4. 只处理染色体名称为 1-23、X、Y（可支持“chr1”“chrX”等格式）的记录，其它染色体直接跳过；
#   5. 每个切分区域会在指定输出路径下创建一个独立的文件夹，并输出 BAM 与 VCF 文件；
#   6. 使用 GNU parallel 并行执行每个区域对应的任务。
#
# 如果暂时不需要调用 Rscript，可将 Rscript 命令行前加上注释（已示例）。


# -m 表示 monitor，开启作业监控模式，也就是 Bash 的作业控制机制。
# 开启后，Bash 能够更好地管理后台任务和进程组。
set -m

# 输入参数
BAM_FILE="$1"
VCF_FILE="$2"
OUTPUT_DIR="$3"
READ_SEPARATION_SCRIPT="$4"

# 检查输入参数
if [ -z "$BAM_FILE" ] || [ -z "$VCF_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$READ_SEPARATION_SCRIPT" ]; then
    echo "Usage: $0 input.bam input.vcf /path/to/output_dir READ_SEPARATION_SCRIPT.R"
    exit 1
fi

# 创建输出路径（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 每个 chunk 的长度（单位：bp）
CHUNK_SIZE=15000000
# VCF 提取区域右侧扩展 bp 数量
VCF_EXTEND=2000
# 并行运行的任务数，请根据实际情况调整
NPROC=25

# 在脚本开头记录开始时间
start_time=$(date +%s)

# 利用 samtools idxstats 获取染色体及其长度（过滤掉行首为 '*' 的行）
samtools idxstats "$BAM_FILE" | awk '$1!="*" {print $1"\t"$2}' > "${OUTPUT_DIR}/chrom_lengths.txt"

# 临时命令文件，用来存放所有子任务的命令
CMD_FILE_MKDIR="${OUTPUT_DIR}/chunk_commands_mkdir.sh"
> "$CMD_FILE_MKDIR"

CMD_FILE_SAM="${OUTPUT_DIR}/chunk_commands_samtools.sh"
> "$CMD_FILE_SAM"

CMD_FILE_VCF="${OUTPUT_DIR}/chunk_commands_bcftools.sh"
> "$CMD_FILE_VCF"

CMD_FILE_READSEPARATION="${OUTPUT_DIR}/chunk_commands_readseparation.sh"
> "$CMD_FILE_READSEPARATION"

# 定义一个函数判断染色体是否合法（符合 1-23, X, Y，不管是否带 chr 前缀）
is_valid_chrom() {
    # 去除可能存在的 "chr" 前缀
    CHROM_NAME=${1#chr}
    if [[ "$CHROM_NAME" =~ ^(1?[0-9]|2[0-3]|[1-9]|X|Y)$ ]]; then
        return 0
    else
        return 1
    fi
}

cleanup() {
    echo "捕获到终止信号，正在杀掉所有子进程..."
    kill 0
    exit 1
}

trap cleanup SIGTERM SIGINT


# 针对每条染色体进行区域切分，生成对应的命令（只处理合法的染色体，输出的目录也在指定输出路径下）
while read -r CHROM LENGTH; do
    # 过滤非指定染色体
    if ! is_valid_chrom "$CHROM"; then
        echo "跳过非指定染色体：$CHROM" 
        continue
    fi

    START=1
    while [ "$START" -le "$LENGTH" ]; do
        END=$((START + CHUNK_SIZE - 1))
        # 如果当前 chunk 超出染色体总长度，则取实际长度
        if [ "$END" -gt "$LENGTH" ]; then
            END="$LENGTH"
        fi

        # 计算双向扩展区域
        VCF_START=$((START - VCF_EXTEND))
        # 确保起始位置不小于 1
        if [ "$VCF_START" -lt 1 ]; then
            VCF_START=1
        fi

        # 针对 VCF 提取区域，右侧扩展 VCF_EXTEND bp
        VCF_END=$((END + VCF_EXTEND))

        # 确保终止位置不超过染色体总长
        if [ "$VCF_END" -gt "$LENGTH" ]; then
            VCF_END="$LENGTH"
        fi

        # 定义当前 chunk 的目录名称，命名格式为 {染色体}-{起始}-{终止}
        CHUNK_DIR="${OUTPUT_DIR}/${CHROM}-${START}-${END}"
        # 定义提取后的 BAM 与 VCF 文件名
        CHUNK_BAM="${CHROM}_${START}-${END}.bam"
        CHUNK_VCF="${CHROM}_${START}-${END}.vcf"

        # 将本次任务的命令写入命令文件
        echo "mkdir -p ${CHUNK_DIR}" >> "$CMD_FILE_MKDIR"
        echo "samtools view -b ${BAM_FILE} ${CHROM}:${START}-${END} > ${CHUNK_DIR}/${CHUNK_BAM}" >> "$CMD_FILE_SAM"
        echo "bcftools view -Ov -r ${CHROM}:${VCF_START}-${VCF_END} ${VCF_FILE} > ${CHUNK_DIR}/${CHUNK_VCF}" >> "$CMD_FILE_VCF"
        echo "if [ -f "${CHUNK_DIR}/${CHUNK_BAM}" ] && [ -f "${CHUNK_DIR}/${CHUNK_VCF}" ]; then \
            Rscript ${READ_SEPARATION_SCRIPT} \
                -b ${CHUNK_DIR}/${CHUNK_BAM} \
                -s ${CHUNK_DIR}/${CHUNK_VCF}; \
            fi" >> "$CMD_FILE_READSEPARATION"

        # 更新下一区间的起始位置
        START=$((END + 1))
    done
done < "${OUTPUT_DIR}/chrom_lengths.txt"

# 利用 GNU parallel 并行执行所有任务
# 增加halt用于杀死所有子进程————当存在10子进程错误的时候,立即停止所有子进程
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_MKDIR"
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_SAM"
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_VCF"
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_READSEPARATION"

# 并行任务执行完毕后记录结束时间
end_time=$(date +%s)
elapsed_seconds=$(( end_time - start_time ))
elapsed_hours=$(echo "scale=2; $elapsed_seconds/3600" | bc)
echo "Total runtime: ${elapsed_hours} hours"
