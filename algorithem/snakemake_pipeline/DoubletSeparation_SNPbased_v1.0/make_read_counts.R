#!/usr/bin/env Rscript

# 加载所需的库
library(data.table)
library(reshape2)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 输入参数：输入的all_out_read.txt文件、输出目录
input_file <- args[1]   # all_out_read.txt 文件路径
out_dir <- args[2]      # 输出目录
 

# 读取数据
print(paste("读取文件:", input_file))
allele_dat <- fread(input_file, header = TRUE, sep = "\t")

# 函数
makeWide <- function(allele_dat, metric = c("reads","UMIs"), level = c("exon","intron","exonintron")){
  dat <- allele_dat
  if(metric == "UMIs"){
    dat <- allele_dat[!UB==""]
  }
  if(level == "exon"){
    dat[, GeneID := GE]
  }else if(level == "intron"){
    dat[, GeneID := GEin]
  }else{
    dat[, GeneID := GE][
      is.na(GeneID), GeneID := GEin]
  }
  dat <- dat[!is.na(GeneID)]
  dat[, c("GE","GEin") := NULL]
  
  cell_lines <- unique(dat$readcall)
  cell_lines <- cell_lines[cell_lines != "other"]
  
  dat_summarized <- dat[, .(reads = length(unique(readID)), UMIs = length(unique(UB))), by =c("RG","GeneID","readcall")]
  
  out_tabs <- list(
    dcast(dat_summarized[readcall == cell_lines[1]], formula = GeneID ~ RG, fill=0, value.var = metric),
    dcast(dat_summarized[readcall == cell_lines[2]], formula = GeneID ~ RG, fill=0, value.var = metric)
  )
  names(out_tabs) <- cell_lines
  return(out_tabs)
}

# 计算并存储结果
print("生成不同的计数表...")
final <- list(
  "counts" = list(
    "exon" = test <- makeWide(allele_dat, metric = "reads", level = "exon"),
    "intron" = test <- makeWide(allele_dat, metric = "reads", level = "intron"),
    "exonintron" = test <- makeWide(allele_dat, metric = "reads", level = "exonintron")
  ),
  
  "umicounts" = list(
    "exon" = test <- makeWide(allele_dat, metric = "UMIs", level = "exon"),
    "intron" = test <- makeWide(allele_dat, metric = "UMIs", level = "intron"),
    "exonintron" = test <- makeWide(allele_dat, metric = "UMIs", level = "exonintron")
  )
)

print("完成计数表的计算.")

# 输出结果到RDS文件
output_file <- paste0(out_dir, "/separated_counts.rds")
print(paste("保存输出结果到:", output_file))
saveRDS(final, file = output_file, compress = FALSE)

print("所有任务完成！")