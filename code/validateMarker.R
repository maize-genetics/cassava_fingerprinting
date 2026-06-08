# Find examples of each genotype type
find_genotype_examples <- function(marker_id, mapping_file, vcf_file, dosage_file) {
  
  library(readr)
  
  # Read files 
  mapping_data <- suppressMessages(read_csv(mapping_file, skip = 6, show_col_types = FALSE))
  mapping_data <- as.data.frame(mapping_data, check.names = FALSE)
  
  vcf_lines <- readLines(vcf_file)
  header_line <- vcf_lines[grepl("^#CHROM", vcf_lines)]
  vcf_header_parts <- strsplit(header_line, "\t")[[1]]
  vcf_sample_names <- vcf_header_parts[10:length(vcf_header_parts)]
  
  dosage_data <- read.table(dosage_file, header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  # Get marker data
  marker_row <- which(mapping_data$AlleleID == marker_id)
  data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
  marker_line <- NULL
  for (i in 1:length(data_lines)) {
    if (grepl(marker_id, data_lines[i], fixed = TRUE)) {
      marker_line <- data_lines[i]
      break
    }
  }
  vcf_parts <- strsplit(marker_line, "\t")[[1]]
  vcf_genotypes <- vcf_parts[10:length(vcf_parts)]
  
  # Get sample columns
  metadata_cols <- c("AlleleID", "CloneID", "AlleleSequenceRef", "AlleleSequenceSnp", 
                     "TrimmedSequenceRef", "TrimmedSequenceSnp", 
                     "Chrom_Cassava_v7Phyto", "ChromPosTag_Cassava_v7Phyto", 
                     "ChromPosSnp_Cassava_v7Phyto", "AlnCnt_Cassava_v7Phyto", 
                     "AlnEvalue_Cassava_v7Phyto", "Strand_Cassava_v7Phyto", 
                     "SNP", "SnpPosition", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", 
                     "AvgPIC", "AvgCountRef", "AvgCountSnp", "RepAvg")
  
  mapping_sample_cols <- names(mapping_data)[!(names(mapping_data) %in% metadata_cols)]
  
  # Find examples of each genotype
  examples_00 <- c()
  examples_01 <- c()  
  examples_11 <- c()
  examples_missing <- c()
  
  for (i in 1:length(vcf_genotypes)) {
    if (vcf_genotypes[i] == "0/0" && length(examples_00) < 2) {
      examples_00 <- c(examples_00, i)
    } else if (vcf_genotypes[i] == "0/1" && length(examples_01) < 2) {
      examples_01 <- c(examples_01, i)
    } else if (vcf_genotypes[i] == "1/1" && length(examples_11) < 2) {
      examples_11 <- c(examples_11, i)
    } else if (vcf_genotypes[i] == "./." && length(examples_missing) < 2) {
      examples_missing <- c(examples_missing, i)
    }
  }
  
  cat("=== GENOTYPE EXAMPLES FOR MARKER", marker_id, "===\n\n")
  
  # Show 0/0 examples
  cat("0/0 EXAMPLES (should be dosage 0 - zero alternate alleles):\n")
  for (pos in examples_00) {
    sample_name <- mapping_sample_cols[pos]
    mapping_code <- mapping_data[[sample_name]][marker_row]
    dosage_val <- dosage_data[[marker_id]][pos]
    cat(sprintf("Sample %s: Mapping=%s, VCF=0/0, Dosage=%s\n", 
                sample_name, mapping_code, dosage_val))
  }
  
  cat("\n0/1 EXAMPLES (should be dosage 1 - one alternate allele):\n")
  for (pos in examples_01) {
    sample_name <- mapping_sample_cols[pos]
    mapping_code <- mapping_data[[sample_name]][marker_row] 
    dosage_val <- dosage_data[[marker_id]][pos]
    cat(sprintf("Sample %s: Mapping=%s, VCF=0/1, Dosage=%s\n", 
                sample_name, mapping_code, dosage_val))
  }
  
  cat("\n1/1 EXAMPLES (should be dosage 2 - two alternate alleles):\n")
  for (pos in examples_11) {
    sample_name <- mapping_sample_cols[pos]
    mapping_code <- mapping_data[[sample_name]][marker_row]
    dosage_val <- dosage_data[[marker_id]][pos]
    cat(sprintf("Sample %s: Mapping=%s, VCF=1/1, Dosage=%s\n", 
                sample_name, mapping_code, dosage_val))
  }
  
  cat("\n./. EXAMPLES (should be dosage NA - missing):\n")
  for (pos in examples_missing) {
    sample_name <- mapping_sample_cols[pos]
    mapping_code <- mapping_data[[sample_name]][marker_row]
    dosage_val <- dosage_data[[marker_id]][pos]
    cat(sprintf("Sample %s: Mapping=%s, VCF=./., Dosage=%s\n", 
                sample_name, mapping_code, dosage_val))
  }
}

# Run it
find_genotype_examples(
  marker_id = "15484497|F|0-51:A>C-51:A>C",
  mapping_file = "data/OrderAppendix_2_DCas22-7517/Report_DCas22-7517_SNP_mapping_2.csv",
  vcf_file = "output/Report_DCas22-7517_SNP_mapping_2.vcf", 
  dosage_file = "output/Report_DCas22-7517_SNP_mapping_2plink_dosageMatrix.txt"
)