
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to the input YAML configuration file"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) {
  stop("Please provide the path to the input configuration file. Usage: Rscript split_yaml_batches.R -i config.yaml -o output_dir")
}

config <- yaml.load_file(opt$input)

num_batches <- as.integer(config$Multi_batches)

for (i in 1:num_batches) {
  new_config <- config
  # Note: Because this pipeline includes the ability to perform read separation directly based on existing BAM and VCF files,
  # this step adds a step to build an index for external BAM files.
  # Build index file path
  bai_file <- paste0(new_config$BAM_file_from_other_source, ".bai")

  if (!file.exists(bai_file)) {
    message("External BAM index file not found. Building index with samtools...")
  
    cmd <- paste("samtools index", shQuote(new_config$BAM_file_from_other_source))
    system(cmd)
  }
  
  # Replace batch-specific configuration
  new_config$sequence_files <- config[[paste0("sequence_files", i)]]
  new_config$annotation <- config[[paste0("annotation", i)]]
  
  # zUMIs_barcode_list configuration
  zUMIs_key <- paste0("zUMIs_barcode_list", i)
  new_config$zUMIs_barcode_list <- config[[zUMIs_key]]
  new_config$barcodes$barcode_file <- config[[zUMIs_key]]
  
  # Read `readseparation_barcode_list`, and if it is empty, assign it to `zUMIs_barcode_list`.
  rs_key <- paste0("readseparation_barcode_list", i)
  new_config$readseparation_barcode_list <- if (!is.null(config[[rs_key]])) {
    config[[rs_key]]
  } else {
    config[[zUMIs_key]]  # If empty, use zUMIs_barcode_list
  }

  # Modify out_dir to the path for each batch
  new_config$out_dir <- file.path(config$out_dir, paste0("batch", i, "_results"))
  # Modify the project name and add the "batch" suffix.
  # To simplify subsequent wildcard replacements, the suffix will not be added for now.
  # new_config$project <- paste0(config$project, "_batch", i) 
  
  # Clean up the original list fields (to ensure batch-specific configurations are not obfuscated)
  for (j in 1:num_batches) {
    new_config[[paste0("sequence_files", j)]] <- NULL
    new_config[[paste0("annotation", j)]] <- NULL
    new_config[[paste0("zUMIs_barcode_list", j)]] <- NULL
    new_config[[paste0("readseparation_barcode_list", j)]] <- NULL
  }
  
  # Output
  output_file <- file.path(opt$output_dir, paste0("batch", i, "_config.yaml"))
  writeLines(as.yaml(new_config), con = output_file)
  cat("✅ Configuration saved:", output_file, "\n")
}
