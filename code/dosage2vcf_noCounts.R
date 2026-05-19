#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(readr)

# Command line argument parsing
option_list <- list(
  make_option(c("-d", "--dosage"), type="character", default=NULL, 
              help="DArT dosage report file (CSV)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output", 
              help="Output file base name (without .vcf extension) [default= %default]", metavar="character"),
  make_option(c("-p", "--ploidy"), type="integer", default=2, 
              help="Species ploidy level [default= %default]", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$dosage)) {
  print_help(opt_parser)
  stop("Dosage file must be specified", call.=FALSE)
}

if (!file.exists(opt$dosage)) {
  stop("Dosage file does not exist", call.=FALSE)
}

# Clean dosage2vcf function - for no counts
dosage2vcf_clean <- function(dart.report, output.file, ploidy = 2) {
  
  cat("Reading dosage report file:", dart.report, "\n")
  
  # Read the dosage report (skip first 6 rows)
  raw_report <- suppressMessages(readr::read_csv(dart.report, skip = 6, show_col_types = FALSE))
  raw_report <- as.data.frame(raw_report, check.names = FALSE)
  
  # Check if we have the basic columns we need
  if (!"AlleleID" %in% names(raw_report)) {
    stop("File must contain AlleleID column")
  }
  
  cat("File contains AlleleID column - proceeding with conversion...\n")
  
  # Use raw_report directly (no need to manipulate it)
  report_data <- raw_report
  
  # Identify metadata and sample columns
  metadata_cols <- c("AlleleID", "CloneID", "AlleleSequenceRef", "AlleleSequenceSnp", 
                     "TrimmedSequenceRef", "TrimmedSequenceSnp", 
                     "Chrom_Cassava_v7Phyto", "ChromPosTag_Cassava_v7Phyto", 
                     "ChromPosSnp_Cassava_v7Phyto", "AlnCnt_Cassava_v7Phyto", 
                     "AlnEvalue_Cassava_v7Phyto", "Strand_Cassava_v7Phyto", 
                     "SNP", "SnpPosition", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", 
                     "AvgPIC", "AvgCountRef", "AvgCountSnp", "RepAvg")
  
  sample_cols <- names(report_data)[!(names(report_data) %in% metadata_cols)]
  
  cat("Found", length(sample_cols), "samples and", nrow(report_data), "markers\n")
  
  # Convert dosage to genotype function
  convert_dosage_to_gt <- function(dosage_val, ploidy) {
    if (is.na(dosage_val) || dosage_val == "" || dosage_val == "-") {
      return("./.")
    }
    
    dosage_val <- as.numeric(dosage_val)
    if (is.na(dosage_val) || dosage_val < 0 || dosage_val > ploidy) {
      return("./.")
    }
    
    # Create genotype string (dosage_val = number of reference alleles)
    ref_count <- dosage_val
    alt_count <- ploidy - dosage_val
    
    # Build genotype: 0 = reference, 1 = alternate
    genotype <- c(rep("0", ref_count), rep("1", alt_count))
    return(paste(genotype, collapse = "/"))
  }
  
  # Extract REF/ALT from SNP column
  get_ref_alt <- function(snp_info) {
    if (!is.na(snp_info) && snp_info != "" && grepl(">", snp_info)) {
      parts <- strsplit(as.character(snp_info), ">")[[1]]
      if (length(parts) == 2) {
        return(list(ref = parts[1], alt = parts[2]))
      }
    }
    return(list(ref = "A", alt = "T"))  # Default
  }
  
  # Create minimal VCF header
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    paste0("##fileDate=", Sys.Date()),
    "##source=dosage2vcf_clean",
    "##reference=Cassava_v7Phyto",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DS,Number=1,Type=Integer,Description="Dosage of reference allele">'
  )
  
  output_vcf <- paste0(output.file, ".vcf")
  
  cat("Converting dosage to genotypes...\n")
  
  vcf_lines <- vcf_header
  
  # Add sample header
  header_line <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_cols), collapse = "\t")
  vcf_lines <- c(vcf_lines, header_line)
  
  # Process each marker
  for (i in 1:nrow(report_data)) {
    if (i %% 1000 == 0) cat("Processing marker", i, "of", nrow(report_data), "\n")
    
    # Use your actual column names
    marker_id <- report_data$AlleleID[i]
    chrom <- report_data$`Chrom_Cassava_v7Phyto`[i]
    pos <- report_data$`ChromPosSnp_Cassava_v7Phyto`[i]
    snp_info <- report_data$SNP[i]
    
    # Handle missing values
    if (is.na(chrom) || chrom == "") chrom <- "Unknown"
    if (is.na(pos) || pos == "") pos <- i
    
    # Get REF/ALT
    ref_alt <- get_ref_alt(snp_info)
    
    # Create genotype strings for each sample
    genotypes <- c()
    for (sample in sample_cols) {
      dosage_val <- report_data[[sample]][i]
      
      gt <- convert_dosage_to_gt(dosage_val, ploidy)
      ds <- if (is.na(dosage_val) || dosage_val == "" || dosage_val == "-") "." else as.character(dosage_val)
      
      # Format: GT:DS
      genotype_info <- paste(gt, ds, sep = ":")
      genotypes <- c(genotypes, genotype_info)
    }
    
    # Create VCF line
    vcf_line <- paste(c(
      chrom,
      pos, 
      marker_id,
      ref_alt$ref,
      ref_alt$alt,
      ".",           # QUAL
      "PASS",        # FILTER
      ".",           # INFO
      "GT:DS",       # FORMAT
      genotypes
    ), collapse = "\t")
    
    vcf_lines <- c(vcf_lines, vcf_line)
  }
  
  # Write VCF
  cat("Writing VCF file:", output_vcf, "\n")
  writeLines(vcf_lines, output_vcf)
  
  cat("Conversion complete!\n")
  return(output_vcf)
}

# Run the conversion
cat("Starting clean dosage-only VCF conversion...\n")
result <- dosage2vcf_clean(opt$dosage, opt$output, opt$ploidy)
cat("Clean VCF created:", result, "\n")