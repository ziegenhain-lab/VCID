#!/usr/bin/env Rscript

# DoubletSeparation_SNPbased_v0.5_chunked.R
# Version for chunked processing, input is a single chunk of BAM and VCF file
# Processes only one chunk, no parallelism or chromosome looping required


# -----------------------------------------------------------------------------
## 1.Basic configuration module
#install.packages(c("tictoc", "optparse"), repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(tictoc)
})

option_list <- list(
  make_option(c("-b", "--bam"), type = "character", 
              help = "Input BAM file (chunked)"),
  make_option(c("-s", "--snps"), type = "character", 
              help = "Input SNP VCF file (chunked)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Output
if (any(is.null(opt$snps),is.null(opt$bam))) {
  stop("All mandatory parameters must be provided. See script usage (--help)")
}


# -----------------------------------------------------------------------------
## 2.function module
# Function: Load reads from bam
load_reads_from_bam <- function(path_bam) {
  taglist <- c("BC", "UB", "GE", "GI")
  whatlist <- c("qname", "pos", "cigar", "seq")

  parms <- ScanBamParam(
            tag = taglist,
            what = whatlist
            )

  dat <- scanBam(file = path_bam, param = parms)
  # Get the information from the first read
  first_read <- dat[[1]]
  # Extract the chromosome name from the first read (rname corresponds to the chromosome name).
  chrom <- first_read$rname[1]

  dat <- dat[[1]]
  # Check all the label content to see if there are any GE or other labels.
  #print(names(dat$tag))  

  # Check if GE and GI exist and if they are empty.
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
  return(dt)
}

variant_parsing_new <- function(reads, variant_positions){
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
  
      # batch base extraction using subseq for DNAStringSet (vectorized)
  tic("Fetching bases from reads...")
  read_seqs <- Biostrings::DNAStringSet(reads$seq[hits_group])
  fetched_bases <- as.character(
    Biostrings::subseq(read_seqs,
                       start = hits_at_in_x,
                       width = 1L)
  )
  toc()
  # assemble output data.table
  out_vars <- data.table(
    obs_base = fetched_bases,
    pos      = var_pos[queryHits(hits)]
  )
  out_vars[, c("RG","UB","GE","GEin","readID") := 
             reads[hits_group, .(RG, UB, GE, GEin, readID)]]
  
  # filter standard bases and rename
  out_vars <- out_vars[obs_base %in% c("A","C","G","T")]
  setnames(out_vars, "pos", "POS")
  
  if (nrow(out_vars) == 0) {
    message("No observations in out_vars, stopping script.")
    quit(save = "no", status = 1)
  }

  return(out_vars)
}


SNPdeconv_new <- function(vcf_chunk, reads_chunk){
  tic("Subprocess, variant_parsing time")
  out_vars <- variant_parsing_new(reads = reads_chunk, variant_positions = vcf_chunk$POS)
  toc()
  
  # Merge VCF information
  out_vars <- merge(out_vars, vcf_chunk, by = "POS")
  
  # Lable REF/ALT
  out_vars[, basecall := fcase(
    obs_base == REF, "REF",
    obs_base == ALT, "ALT",
    default = "other"
  )]
  
  # Get cell line column name
  cell_lines <- tail(colnames(vcf_chunk), 2)
  cl1 <- cell_lines[1]; cl2 <- cell_lines[2]
  
  # Vectorized decision-making: data.table::fcase supports multiple conditions.
  tic("Subprocess, time consumption of vectorized basecall2")
  out_vars[, basecall2 := fcase(
    # REF Conditions
    basecall == "REF" & get(cl1) %in% c("./.", ".|.", "0/0", "0|0") & get(cl2) %in% c("1/1", "1|1"), cl1,
    basecall == "REF" & get(cl1) %in% c("1/1", "1|1") & get(cl2) %in% c("./.", ".|.", "0/0", "0|0"), cl2,
    basecall == "REF" & get(cl1) %in% c("1/0", "0/1", "0|1", "1|0") & get(cl2) %in% c("1/1", "1|1"), cl1,
    basecall == "REF" & get(cl1) %in% c("1/1", "1|1") & get(cl2) %in% c("1/0", "0/1", "0|1", "1|0"), cl2,
    # ALT Conditions
    basecall == "ALT" & get(cl1) %in% c("./.", ".|.", "0/0", "0|0") & get(cl2) %in% c("1/1", "1|1"), cl2,
    basecall == "ALT" & get(cl1) %in% c("1/1", "1|1") & get(cl2) %in% c("./.", ".|.", "0/0", "0|0"), cl1,
    basecall == "ALT" & get(cl1) %in% c("./.", ".|.", "0/0", "0|0") & get(cl2) %in% c("1/0", "0/1", "0|1", "1|0"), cl2,
    basecall == "ALT" & get(cl1) %in% c("1/0", "0/1", "0|1", "1|0") & get(cl2) %in% c("./.", ".|.", "0/0", "0|0"), cl1,
    default = "other"
  )]
  toc()
  
  # Finally, a majority vote is conducted by grouping by read.
  tic("Sub-process, read_decision & summarizing time consumption")
  out_reads <- out_vars[, .(
    readcall = {
      tb <- tabulate(match(basecall2, unique(basecall2)))
      ux <- unique(basecall2)
      maj <- ux[which.max(tb)]
      if (max(tb)/sum(tb) >= 0.66) maj else "other"
    },
    nVars = uniqueN(POS)
  ), by = .(RG, UB, GE, GEin, readID)]
  toc()
  
  return(out_reads)
}



# -----------------------------------------------------------------------------
## 3. Main process 
out_prefix <- sub("\\.bam$", "", opt$bam)

# ----------------------------
tic("Read VCF chunk timing")
message("Reading VCF chunk...")
if (grepl(opt$snps, pattern = ".gz$")) {
  vcf <- fread(cmd = paste("zcat", opt$snps, "| grep -v '^##'"), header = TRUE)
} else {
  vcf <- fread(cmd = paste("grep -v '^##'", opt$snps), header = TRUE)
}
setnames(vcf, c("CHROM", colnames(vcf)[2:length(colnames(vcf))]))
vcf <- vcf[nchar(REF) == 1 & nchar(ALT) == 1]
toc()
message("Done reading VCF chunk.")

# ----------------------------
tic("Load BAM reads:")
message("Loading reads from BAM...")
reads_chunk <- load_reads_from_bam(path_bam = opt$bam)
toc()

# Check if it is empty
if (nrow(reads_chunk) == 0) {
  cat("No data in reads_chunk. Exiting the program...\n")
  stop("Exiting due to empty reads_chunk")
}
message(sprintf("Loaded %d reads.", nrow(reads_chunk)))

# ----------------------------
tic("SNP-based read separation:")
message("Performing SNP-based read separation...")
out_reads <- SNPdeconv_new(vcf_chunk = vcf, reads_chunk = reads_chunk)
toc()

# ----------------------------
tic("Write output files")
# Output
fwrite(out_reads, file = paste0(out_prefix, "_out_read.txt"), sep = "\t", quote = FALSE)
fwrite(out_reads[, .(readID, readcall)], file = paste0(out_prefix, "_out_read_readID.txt"), sep = "\t", quote = FALSE)
toc()

message("Chunk processing complete.")




