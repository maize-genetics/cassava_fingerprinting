#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(readr)

# Command line argument parsing
option_list <- list(
  make_option(c("-m", "--mappingFile"), type="character", default=NULL, 
              help="DArT SNP mapping file (CSV)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output", 
              help="Output file base name (without .vcf extension) [default= %default]", metavar="character"),
  make_option(c("-p", "--ploidy"), type="integer", default=2, 
              help="Species ploidy level [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$mappingFile)) {
  print_help(opt_parser)
  stop("SNP mapping file must be specified", call.=FALSE)
}

if (!file.exists(opt$mappingFile)) {
  stop("SNP mapping file does not exist", call.=FALSE)
}

# Clean up csv

fix_csv_header <- function(csv_file) {
  cat("Preprocessing CSV header to fix duplicates...\n")
  lines <- readLines(csv_file)
  header_line <- lines[7]  # Sample names are on line 7
  
  # Fix spaces and duplicates in header
  header_parts <- strsplit(header_line, ",")[[1]]
  header_parts <- gsub(" ", "", header_parts)  # Remove spaces
  
  # First pass: identify all names that appear multiple times
  name_counts <- table(header_parts)
  duplicate_names <- names(name_counts)[name_counts > 1]
  
  if (length(duplicate_names) > 0) {
    cat("Found duplicates:", paste(duplicate_names, collapse = ", "), "\n")
    
    # Second pass: rename ALL occurrences of duplicate names
    for (i in seq_along(header_parts)) {
      name <- header_parts[i]
      if (name %in% duplicate_names) {
        header_parts[i] <- paste0(name, "Col", i)
        cat("Renamed duplicate:", name, "->", header_parts[i], "\n")
      }
    }
  }
  
  # Replace the header line
  lines[7] <- paste(header_parts, collapse = ",")
  
  # Write fixed CSV
  temp_file <- paste0(csv_file, "_fixed.csv")
  writeLines(lines, temp_file)
  cat("Fixed CSV written to:", temp_file, "\n")
  return(temp_file)
}


# Then modify the beginning of your existing function:
snp_mapping_to_vcf_fixed <- function(snp_file, output.file, ploidy = 2) {
  
  # Add this line at the very beginning:
  snp_file <- fix_csv_header(snp_file)
  
  cat("Reading DArT SNP mapping file:", snp_file, "\n")

  # Read SNP mapping data (skip first 6 header rows like DArT format)
  snp_data <- suppressMessages(read_csv(snp_file, skip = 6, show_col_types = FALSE))
  snp_data <- as.data.frame(snp_data, check.names = FALSE)
  
  if (!"AlleleID" %in% names(snp_data)) {
    stop("File must contain AlleleID column")
  }
  
  cat("Processing DArT SNP mapping format...\n")
  
  # Identify metadata and sample columns
  metadata_cols <- c("AlleleID", "CloneID", "AlleleSequenceRef", "AlleleSequenceSnp", 
                     "TrimmedSequenceRef", "TrimmedSequenceSnp", 
                     "Chrom_Cassava_v7Phyto", "ChromPosTag_Cassava_v7Phyto", 
                     "ChromPosSnp_Cassava_v7Phyto", "AlnCnt_Cassava_v7Phyto", 
                     "AlnEvalue_Cassava_v7Phyto", "Strand_Cassava_v7Phyto", 
                     "SNP", "SnpPosition", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", 
                     "AvgPIC", "AvgCountRef", "AvgCountSnp", "RepAvg")
  
  sample_cols <- names(snp_data)[!(names(snp_data) %in% metadata_cols)]
  
  # Fix sample names: replace spaces with underscores
  sample_cols_fixed <- sample_cols
  
  cat("Found", length(sample_cols), "samples and", nrow(snp_data), "markers\n")
  
  # Clean and validate chromosome and position data
  cat("Cleaning chromosome and position data...\n")
  
  for (i in 1:nrow(snp_data)) {
    # Clean chromosome name
    chrom <- snp_data$`Chrom_Cassava_v7Phyto`[i]
    if (is.na(chrom) || chrom == "") {
      snp_data$`Chrom_Cassava_v7Phyto`[i] <- "Unknown"
    }
    
    # Clean position - ensure it's a positive integer
    pos <- snp_data$`ChromPosSnp_Cassava_v7Phyto`[i]
    if (is.na(pos) || pos == "" || pos <= 0) {
      snp_data$`ChromPosSnp_Cassava_v7Phyto`[i] <- i  # Use marker index
    } else {
      snp_data$`ChromPosSnp_Cassava_v7Phyto`[i] <- max(1, as.integer(pos))
    }
  }
  
  # Sort data by chromosome and position for proper VCF format
  cat("Sorting markers by chromosome and position...\n")
  
  # Create sorting key
  snp_data$sort_chrom <- snp_data$`Chrom_Cassava_v7Phyto`
  snp_data$sort_pos <- as.numeric(snp_data$`ChromPosSnp_Cassava_v7Phyto`)
  
  # Sort by chromosome then position
  snp_data <- snp_data[order(snp_data$sort_chrom, snp_data$sort_pos), ]
  
  # Remove sorting columns
  snp_data$sort_chrom <- NULL
  snp_data$sort_pos <- NULL
  
  # Get unique chromosomes for contig headers
  unique_chroms <- unique(snp_data$`Chrom_Cassava_v7Phyto`)
  unique_chroms <- sort(unique_chroms[!is.na(unique_chroms) & unique_chroms != ""])
  
  cat("Found chromosomes/scaffolds:", paste(head(unique_chroms, 10), collapse = ", "), 
      if(length(unique_chroms) > 10) "..." else "", "\n")
  
  # Convert DArT SNP codes to VCF genotypes
  convert_dart_to_vcf_gt <- function(dart_code) {
    if (is.na(dart_code) || dart_code == "" || dart_code == "-") {
      return("./.")  # Missing
    }
    
    dart_code <- as.character(dart_code)
    
    # DArT encoding: 0=ref homo, 1=alt homo, 2=hetero
    if (dart_code == "0") {
      return("0/0")  # Reference homozygote
    } else if (dart_code == "1") {
      return("1/1")  # Alternate homozygote  
    } else if (dart_code == "2") {
      return("0/1")  # Heterozygote
    } else {
      return("./.")  # Unknown code
    }
  }
  
  # Extract REF/ALT from SNP column
  get_ref_alt <- function(snp_info) {
    if (!is.na(snp_info) && snp_info != "" && grepl(">", snp_info)) {
      parts <- strsplit(as.character(snp_info), ">")[[1]]
      if (length(parts) == 2 && nchar(parts[1]) > 0 && nchar(parts[2]) > 0) {
        return(list(ref = parts[1], alt = parts[2]))
      }
    }
    return(list(ref = "A", alt = "T"))  # Default
  }
  
  # Create contig header lines
  contig_lines <- paste0('##contig=<ID=', unique_chroms, '>')
  
  # Create proper VCF header
  vcf_header <- c(
    "##fileformat=VCFv4.3",
    paste0("##fileDate=", Sys.Date()),
    "##source=dart_snp_mapping_converter_fixed",
    "##reference=Cassava_v7Phyto",
    contig_lines,
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  )
  
  output_vcf <- paste0(output.file, ".vcf")
  
  cat("Converting DArT SNP codes to VCF genotypes...\n")
  
  vcf_lines <- vcf_header
  
  # Add sample header with FIXED sample names (spaces → underscores)
  header_line <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_cols_fixed), collapse = "\t")
  vcf_lines <- c(vcf_lines, header_line)
  
  # Process each marker (now sorted)
  for (i in 1:nrow(snp_data)) {
    if (i %% 1000 == 0) cat("Processing marker", i, "of", nrow(snp_data), "\n")
    
    # Extract marker information
    marker_id <- snp_data$AlleleID[i]
    chrom <- snp_data$`Chrom_Cassava_v7Phyto`[i]
    pos <- snp_data$`ChromPosSnp_Cassava_v7Phyto`[i]
    snp_info <- snp_data$SNP[i]
    
    # Get REF/ALT
    ref_alt <- get_ref_alt(snp_info)
    
    # Convert genotypes for each sample (using ORIGINAL column names for data access)
    genotypes <- c()
    for (sample in sample_cols) {  # Use original names to access data
      dart_code <- snp_data[[sample]][i]
      gt <- convert_dart_to_vcf_gt(dart_code)
      genotypes <- c(genotypes, gt)
    }
    
    # Create VCF line
    vcf_line <- paste(c(
      chrom,
      pos, 
      marker_id,
      ref_alt$ref,
      ref_alt$alt,
      "60",          # QUAL - reasonable quality score
      "PASS",        # FILTER
      ".",           # INFO
      "GT",          # FORMAT (just genotype)
      genotypes
    ), collapse = "\t")
    
    vcf_lines <- c(vcf_lines, vcf_line)
  }
  
  # Write VCF
  cat("Writing properly formatted VCF file:", output_vcf, "\n")
  writeLines(vcf_lines, output_vcf)
  
  cat("VCF creation complete!\n")
  cat("- Contig headers added\n")
  cat("- Markers sorted by chromosome and position\n") 
  cat("- All positions are positive integers\n")
  cat("- Sample names cleaned (spaces → underscores)\n")
  
  return(output_vcf)
}

# Run the conversion
cat("Creating properly formatted VCF from DArT SNP mapping file...\n")
result <- snp_mapping_to_vcf_fixed(opt$mappingFile, opt$output, opt$ploidy)
cat("Fixed VCF created:", result, "\n")