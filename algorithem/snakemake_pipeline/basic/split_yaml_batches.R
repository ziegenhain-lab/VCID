
# 加载依赖
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "输入的 YAML 配置文件路径"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".", help = "输出目录")
)

# 解析参数
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) {
  stop("请提供输入配置文件路径，用法: Rscript split_yaml_batches.R -i config.yaml -o output_dir")
}

# 加载 YAML
config <- yaml.load_file(opt$input)

# 获取批次数量
num_batches <- as.integer(config$Multi_batches)




# 遍历每个 batch
for (i in 1:num_batches) {
  new_config <- config
  # 注释：因为本pipeline加入了直接基于已有的bam文件和vcf文件进行read separation的功能，
  # 因此这一步中增加了对外来bam文件建立索引的步骤
  # 构建索引文件路径
  bai_file <- paste0(new_config$BAM_file_from_other_source, ".bai")

  # 判断索引文件是否存在
  if (!file.exists(bai_file)) {
    message("External BAM index file not found. Building index with samtools...")
  
    # 构建索引
    cmd <- paste("samtools index", shQuote(new_config$BAM_file_from_other_source))
    system(cmd)
  }
  
  # 替换 batch-specific 配置
  new_config$sequence_files <- config[[paste0("sequence_files", i)]]
  new_config$annotation <- config[[paste0("annotation", i)]]
  
  # zUMIs_barcode_list 配置
  zUMIs_key <- paste0("zUMIs_barcode_list", i)
  new_config$zUMIs_barcode_list <- config[[zUMIs_key]]
  new_config$barcodes$barcode_file <- config[[zUMIs_key]]
  
  # 读取 readseparation_barcode_list，如果为空则赋值为 zUMIs_barcode_list
  rs_key <- paste0("readseparation_barcode_list", i)
  new_config$readseparation_barcode_list <- if (!is.null(config[[rs_key]])) {
    config[[rs_key]]
  } else {
    config[[zUMIs_key]]  # 如果为空则用 zUMIs_barcode_list
  }

  # 修改 out_dir 为每个 batch 的路径
  new_config$out_dir <- file.path(config$out_dir, paste0("batch", i, "_results"))
  # 修改 project, 增加 batch 的后缀
  # new_config$project <- paste0(config$project, "_batch", i) # 为了简化后续的wildcard替换，暂时不加后缀
  
  # 清理原始列表字段（确保 batch-specific 配置不会混淆）
  for (j in 1:num_batches) {
    new_config[[paste0("sequence_files", j)]] <- NULL
    new_config[[paste0("annotation", j)]] <- NULL
    new_config[[paste0("zUMIs_barcode_list", j)]] <- NULL
    new_config[[paste0("readseparation_barcode_list", j)]] <- NULL
  }
  
  # 这里可以保存或使用新的 new_config 进行其他操作，如写入文件
  # 写入输出 YAML
  output_file <- file.path(opt$output_dir, paste0("batch", i, "_config.yaml"))
  writeLines(as.yaml(new_config), con = output_file)
  cat("✅ 配置已保存:", output_file, "\n")
}
