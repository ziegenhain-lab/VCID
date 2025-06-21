#!/usr/bin/env Rscript

# DoubletSeparation_SNPbased_v0.5_chunked.R
# 适用于分块处理的版本，输入为一个 chunk 的 BAM 与 VCF 文件
# 仅处理一个 chunk，无需并行与染色体循环


# -----------------------------------------------------------------------------
## 1.基础配置模块
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
  library(GenomicAlignments)
})

option_list <- list(
  make_option(c("-b", "--bam"), type = "character", 
              help = "Input BAM file (chunked)"),
  make_option(c("-s", "--snps"), type = "character", 
              help = "Input SNP VCF file (chunked)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 创建输出目录
if (any(is.null(opt$snps),is.null(opt$bam))) {
  stop("All mandatory parameters must be provided. See script usage (--help)")
}


# -----------------------------------------------------------------------------
## 2.函数模块
# 函数：从 bam 中加载 reads
load_reads_from_bam <- function(path_bam) {
  taglist <- c("BC", "UB", "GE", "GI")
  whatlist <- c("qname", "pos", "cigar", "seq")

  parms <- ScanBamParam(
            tag = taglist,
            what = whatlist
            )

  dat <- scanBam(file = path_bam, param = parms)
  # 获取第一个 read 的信息
  first_read <- dat[[1]]
  # 提取第一个 read 的染色体名称（rname 对应染色体名称）
  chrom <- first_read$rname[1]

  dat <- dat[[1]]
  # 检查所有的标签内容，看看是否有 GE 和其他标签
  #print(names(dat$tag))  # 打印所有的标签名

    # 检查 GE 和 GI 是否存在，并且是否为空
if ( (!"GE" %in% names(dat$tag) || all(is.na(dat$tag$GE))) &&
     (!"GI" %in% names(dat$tag) || all(is.na(dat$tag$GI))) ) {
  stop("Both GE and GI in all reads of this chunk are empty. Ending the program.")
}

  dt <- data.table(
    readID = dat$qname,
    chrom  = chrom,
    pos    = dat$pos,
    cigar  = dat$cigar,
    seq    = as.character(dat$seq),
    RG     = dat$tag$BC,
    UB     = dat$tag$UB,
    GE     = if ("GE" %in% names(dat$tag)) dat$tag$GE else rep(NA_character_, length(dat$qname)),
    GEin   = if ("GI" %in% names(dat$tag)) dat$tag$GI else rep(NA_character_, length(dat$qname))
  )

  # rm删除变量，但 R 的内存回收机制并不会立即清理已删除对象占用的内存空间。
  # 这时需要使用 gc() 来强制进行垃圾回收，释放内存。
  #rm(dat)
  #gc(verbose = F)
  #dt <- dt[ !(is.na(GE) & is.na(GEin)) ]
  return(dt)
}


variant_parsing <- function(reads, variant_positions){
  #parse all cigars to reference seq
  ops <- c("M", "=", "X")
  ranges_on_ref <- cigarRangesAlongReferenceSpace(reads$cigar, pos=reads$pos, ops=ops)
  ranges_on_query <- cigarRangesAlongQuerySpace(reads$cigar, ops=ops)
  gc(verbose = F)
  range_group <- togroup(PartitioningByWidth(ranges_on_ref))
  ranges_on_ref <- unlist(ranges_on_ref, use.names=FALSE)
  ranges_on_query <- unlist(ranges_on_query, use.names=FALSE)
  query2ref_shift <- start(ranges_on_ref) - start(ranges_on_query)
  
  var_pos <- variant_positions
  hits <- findOverlaps(var_pos, ranges_on_ref)
  hits_at_in_x <- var_pos[queryHits(hits)] - query2ref_shift[subjectHits(hits)]
  hits_group <- range_group[subjectHits(hits)]
  fetched_bases <- subseq(reads[hits_group,]$seq, start=hits_at_in_x, width=1L)
  
  #now add everything together in the output data.table
  out_vars <- data.table(
    obs_base = fetched_bases,
    pos = var_pos[queryHits(hits)]
  )
  out_vars[, c("RG","UB","GE","GEin","readID") := reads[hits_group, c("RG","UB", "GE","GEin", "readID"), with = F] ]
  
  out_vars <- out_vars[obs_base %in% c("A","C","G","T") ]
  setnames(out_vars,"pos","POS")
  
  if (nrow(out_vars) == 0) {
      message("No observations in out_vars, stopping script.")
      quit(save = "no", status = 1)
  }

  return(out_vars)

}

read_decision <- function(basecalls){
  if(length(basecalls) == 1){
    return(basecalls)
  }else{
    ux <- unique(basecalls)
    basecall_summary <- tabulate(match(basecalls, ux))
    names(basecall_summary) <- ux
    majority_basecall <- ux[which.max(basecall_summary)]
    if(basecall_summary[majority_basecall]/sum(basecall_summary) >= 0.66){
      return(majority_basecall)
    }else{
      return("other")
    }
  }
}

SNPdeconv <- function(vcf_chunk, reads_chunk){
  
  #get the bases that were observed in our reads according to the positions in VCF
  out_vars <- variant_parsing(reads = reads_chunk, variant_positions = vcf_chunk$POS)

  #combine VCF with observations
  out_vars <- merge(out_vars,vcf_chunk,by = "POS" )
 
  #categorize:
  out_vars[              , basecall := "other"][
    obs_base == REF, basecall := "REF"][
      obs_base == ALT, basecall := "ALT"]
  #注释：也许我们可以省掉这个输出，太多输出没有必要
  #fwrite(out_vars, file = paste0(out_prefix, "_out_vars.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  #translate ref/alt into cell line:
  cell_lines <- tail(colnames(vcf_chunk),2)
  
  # The default setting is "other"
  out_vars[, basecall2 := "other"]
  
  for (i in seq_len(nrow(out_vars))) {
    # Get the basecall and cell_lines values ​​for the current row
    current_basecall <- out_vars[i, basecall]
    genotype_1 <- out_vars[i, get(cell_lines[1])]
    genotype_2 <- out_vars[i, get(cell_lines[2])]
    
    # Logic of Decision Tree 
    if (current_basecall == "REF") {
      if (genotype_1 %in% c("./.", ".|.", "0/0", "0|0") & genotype_2 %in% c("1/1", "1|1")) {
        out_vars[i, basecall2 := cell_lines[1]]
      } else if (genotype_1 %in% c("1/1", "1|1") & genotype_2 %in% c("./.", ".|.", "0/0", "0|0")) {
        out_vars[i, basecall2 := cell_lines[2]]
      } else if (genotype_1 %in% c("1/0", "0/1", "0|1", "1|0") & genotype_2 %in% c("1/1", "1|1")) {
        out_vars[i, basecall2 := cell_lines[1]]
      } else if (genotype_1 %in% c("1/1", "1|1") & genotype_2 %in% c("1/0", "0/1", "0|1", "1|0")) {
        out_vars[i, basecall2 := cell_lines[2]]
      }
    } else if (current_basecall == "ALT") {
      if (genotype_1 %in% c("./.", ".|.", "0/0", "0|0") & genotype_2 %in% c("1/1", "1|1")) {
        out_vars[i, basecall2 := cell_lines[2]]
      } else if (genotype_1 %in% c("1/1", "1|1") & genotype_2 %in% c("./.", ".|.", "0/0", "0|0")) {
        out_vars[i, basecall2 := cell_lines[1]]
      } else if (genotype_1 %in% c("./.", ".|.", "0/0", "0|0") & genotype_2 %in% c("1/0", "0/1", "0|1", "1|0")) {
        out_vars[i, basecall2 := cell_lines[2]]
      } else if (genotype_1 %in% c("1/0", "0/1", "0|1", "1|0") & genotype_2 %in% c("./.", ".|.", "0/0", "0|0")) {
        out_vars[i, basecall2 := cell_lines[1]]
      }
    }
  }
  
  out_reads <- out_vars[, .(readcall = read_decision(basecall2), nVars =length(unique(POS))), by = c("RG","UB","GE","GEin","readID")]
  
  return(out_reads)
}



# -----------------------------------------------------------------------------
## 3.主流程 
out_prefix <- sub("\\.bam$", "", opt$bam)

# 加载 VCF（chunk 区域）
message("Reading VCF chunk...")
if (grepl(opt$snps, pattern = ".gz$")) {
  vcf <- fread(cmd = paste("zcat", opt$snps, "| grep -v '^##'"), header = TRUE)
} else {
  vcf <- fread(cmd = paste("grep -v '^##'", opt$snps), header = TRUE)
}
setnames(vcf, c("CHROM", colnames(vcf)[2:length(colnames(vcf))]))
vcf <- vcf[nchar(REF) == 1 & nchar(ALT) == 1]
message("Done reading VCF chunk.")

# 加载 BAM（chunk 区域）
message("Loading reads from BAM...")
reads_chunk <- load_reads_from_bam(path_bam = opt$bam)
#str(reads_chunk)
#print(reads_chunk$GE)
#print(reads_chunk$GEin)

# 假设 reads_chunk 是你处理的数据
if (nrow(reads_chunk) == 0) {
  cat("No data in reads_chunk. Exiting the program...\n")
  stop("Exiting due to empty reads_chunk")
}

message(sprintf("Loaded %d reads.", nrow(reads_chunk)))

message("Performing SNP-based read separation...")
out_reads <- SNPdeconv(vcf_chunk = vcf, reads_chunk = reads_chunk)

# 写入结果
fwrite(out_reads, file = paste0(out_prefix, "_out_read.txt"), sep = "\t", quote = FALSE)
fwrite(out_reads[, .(readID, readcall)], file = paste0(out_prefix, "_out_read_readID.txt"), sep = "\t", quote = FALSE)

message("Chunk processing complete.")




