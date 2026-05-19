#!/usr/bin/env Rscript

# Load required libraries
library(optparse)

# Command line argument parsing
option_list <- list(
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="VCF file to validate", metavar="character"),
  make_option(c("-m", "--mappingFile"), type="character", default=NULL, 
              help="Original DArT SNP mapping file (CSV)", metavar="character"),
  make_option(c("-s", "--samples"), type="integer", default=10, 
              help="Number of random samples to check [default= %default]", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$vcf) || is.null(opt$mappingFile)) {
  print_help(opt_parser)
  stop("Both VCF file and original mapping file must be specified", call.=FALSE)
}

if (!file.exists(opt$vcf)) {
  stop("VCF file does not exist", call.=FALSE)
}

if (!file.exists(opt$mappingFile)) {
  stop("Mapping file does not exist", call.=FALSE)
}

validate_vcf <- function(vcf_file, mapping_file, num_samples = 10) {
  
  cat("=== VCF Validation Report ===\n")
  cat("VCF file:", vcf_file, "\n")
  cat("Original mapping file:", mapping_file, "\n\n")
  
  # Read VCF file
  cat("Reading VCF file...\n")
  vcf_lines <- readLines(vcf_file)
  
  # Find header and data lines
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]
  column_header_line <- vcf_lines[grepl("^#CHROM", vcf_lines)]
  data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
  
  cat("VCF contains:", length(header_lines), "header lines,", length(data_lines), "data lines\n")
  
  # Parse column header
  vcf_header <- strsplit(column_header_line, "\t")[[1]]
  vcf_samples <- vcf_header[10:length(vcf_header)]  # Samples start at column 10
  
  # Read original mapping file
  cat("Reading original mapping file...\n")
  library(readr)
  mapping_data <- suppressMessages(read_csv(mapping_file, skip = 6, show_col_types = FALSE))
  mapping_data <- as.data.frame(mapping_data, check.names = FALSE)
  
  # Identify sample columns in mapping file
  metadata_cols <- c("AlleleID", "CloneID", "AlleleSequenceRef", "AlleleSequenceSnp", 
                     "TrimmedSequenceRef", "TrimmedSequenceSnp", 
                     "Chrom_Cassava_v7Phyto", "ChromPosTag_Cassava_v7Phyto", 
                     "ChromPosSnp_Cassava_v7Phyto", "AlnCnt_Cassava_v7Phyto", 
                     "AlnEvalue_Cassava_v7Phyto", "Strand_Cassava_v7Phyto", 
                     "SNP", "SnpPosition", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", 
                     "AvgPIC", "AvgCountRef", "AvgCountSnp", "RepAvg")
  
  mapping_samples <- names(mapping_data)[!(names(mapping_data) %in% metadata_cols)]
  
  cat("Original mapping file contains:", nrow(mapping_data), "markers,", length(mapping_samples), "samples\n")
  
  # Check 1: Sample consistency
  cat("\n=== Check 1: Sample Consistency ===\n")
  if (length(vcf_samples) == length(mapping_samples)) {
    cat("✓ Sample count matches:", length(vcf_samples), "samples\n")
  } else {
    cat("✗ Sample count mismatch! VCF:", length(vcf_samples), "vs Mapping:", length(mapping_samples), "\n")
  }
  
  if (all(vcf_samples %in% mapping_samples)) {
    cat("✓ All VCF samples found in mapping file\n")
  } else {
    missing <- setdiff(vcf_samples, mapping_samples)
    cat("✗ VCF samples not in mapping file:", paste(head(missing, 5), collapse = ", "), "\n")
  }
  
  # Check 2: Marker consistency
  cat("\n=== Check 2: Marker Consistency ===\n")
  if (length(data_lines) == nrow(mapping_data)) {
    cat("✓ Marker count matches:", length(data_lines), "markers\n")
  } else {
    cat("✗ Marker count mismatch! VCF:", length(data_lines), "vs Mapping:", nrow(mapping_data), "\n")
  }
  
  # Check 3: Genotype conversion accuracy
  cat("\n=== Check 3: Genotype Conversion Accuracy ===\n")
  
  # Function to convert DArT code to expected VCF genotype
  dart_to_vcf <- function(dart_code) {
    if (is.na(dart_code) || dart_code == "" || dart_code == "-") return("./.")
    dart_code <- as.character(dart_code)
    if (dart_code == "0") return("0/0")      # Reference homozygote
    else if (dart_code == "1") return("1/1") # Alternate homozygote
    else if (dart_code == "2") return("0/1") # Heterozygote
    else return("./.")
  }
  
  # Check random samples
  sample_indices <- sample(1:min(length(data_lines), nrow(mapping_data)), num_samples)
  
  errors <- 0
  total_checks <- 0
  
  for (i in sample_indices) {
    vcf_line <- strsplit(data_lines[i], "\t")[[1]]
    vcf_marker_id <- vcf_line[3]  # ID column
    vcf_genotypes <- vcf_line[10:length(vcf_line)]  # Sample genotypes
    
    # Find matching marker in mapping data
    mapping_row <- which(mapping_data$AlleleID == vcf_marker_id)
    
    if (length(mapping_row) == 0) {
      cat("✗ Marker", vcf_marker_id, "not found in mapping file\n")
      errors <- errors + 1
      next
    }
    
    # Check a few samples for this marker
    check_samples <- sample(1:min(length(vcf_samples), length(mapping_samples)), 3)
    
    for (j in check_samples) {
      sample_name <- vcf_samples[j]
      
      if (sample_name %in% mapping_samples) {
        dart_code <- mapping_data[mapping_row, sample_name]
        expected_gt <- dart_to_vcf(dart_code)
        actual_gt <- vcf_genotypes[j]
        
        total_checks <- total_checks + 1
        
        if (expected_gt != actual_gt) {
          cat("✗ Marker", vcf_marker_id, "Sample", sample_name, 
              "- Expected:", expected_gt, "Got:", actual_gt, 
              "(DArT code:", dart_code, ")\n")
          errors <- errors + 1
        }
      }
    }
  }
  
  success_rate <- (total_checks - errors) / total_checks * 100
  cat("Genotype conversion accuracy:", sprintf("%.1f%%", success_rate), 
      "(", total_checks - errors, "/", total_checks, "correct )\n")
  
  # Check 4: VCF format validation
  cat("\n=== Check 4: VCF Format Validation ===\n")
  
  # Check header format
  if (any(grepl("##fileformat=VCF", header_lines))) {
    cat("✓ VCF format declaration found\n")
  } else {
    cat("✗ Missing VCF format declaration\n")
  }
  
  # Check required columns
  required_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  if (all(required_cols %in% vcf_header[1:9])) {
    cat("✓ All required VCF columns present\n")
  } else {
    missing_cols <- setdiff(required_cols, vcf_header[1:9])
    cat("✗ Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  }
  
  # Check a few data lines for format issues
  cat("\n=== Check 5: Data Format Issues ===\n")
  format_errors <- 0
  
  for (i in 1:min(10, length(data_lines))) {
    line_parts <- strsplit(data_lines[i], "\t")[[1]]
    
    # Check correct number of columns
    if (length(line_parts) != length(vcf_header)) {
      cat("✗ Line", i, "has", length(line_parts), "columns, expected", length(vcf_header), "\n")
      format_errors <- format_errors + 1
    }
    
    # Check genotype format
    genotypes <- line_parts[10:length(line_parts)]
    invalid_gt <- !grepl("^(\\d/\\d|\\./\\.)$", genotypes)
    if (any(invalid_gt)) {
      cat("✗ Line", i, "has invalid genotype formats:", 
          paste(head(genotypes[invalid_gt], 3), collapse = ", "), "\n")
      format_errors <- format_errors + 1
    }
  }
  
  if (format_errors == 0) {
    cat("✓ No format errors detected in sample data lines\n")
  }
  
  # Summary
  cat("\n=== VALIDATION SUMMARY ===\n")
  if (errors == 0 && format_errors == 0) {
    cat("🎉 VCF file appears to be correctly converted!\n")
  } else {
    cat("⚠️  Issues detected:", errors, "genotype errors,", format_errors, "format errors\n")
  }
  
  cat("Total genotype checks performed:", total_checks, "\n")
  cat("Conversion accuracy:", sprintf("%.1f%%", success_rate), "\n")
}

# Run validation
validate_vcf(opt$vcf, opt$mappingFile, opt$samples)