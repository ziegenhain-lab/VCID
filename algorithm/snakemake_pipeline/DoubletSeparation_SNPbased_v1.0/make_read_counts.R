#!/usr/bin/env Rscript

library(data.table)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]   # Path to the all_out_read.txt file
out_dir <- args[2]     
 
print(paste("Reading file:", input_file))
allele_dat <- fread(input_file, header = TRUE, sep = "\t")

# Function to create wide format count tables
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

# Calculate and store the results
print("Generate different counting tables...")
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

print("Complete the calculations on the counting table..")

# Output
output_file <- paste0(out_dir, "/separated_counts.rds")
print(paste("Save the output to:", output_file))
saveRDS(final, file = output_file, compress = FALSE)

print("All tasks completed.！")