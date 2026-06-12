# Convert dosage file from PLINK .traw file into a matrix 

library(optparse)

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="PLINK .traw file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="dosage_matrix.txt",
              help="Output dosage file [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input PLINK .traw file must be specified", call.=FALSE)
}

format_plink_to_dosage_matrix <- function(traw_file, output_file) {
  
  cat("Reading PLINK .traw file:", traw_file, "\n")
  
  # Read PLINK .traw file
  traw_data <- read.table(traw_file, header = TRUE, stringsAsFactors = FALSE, 
                          sep = "\t", check.names = FALSE)
  
  cat("Found", nrow(traw_data), "markers and", (ncol(traw_data) - 6), "samples\n")
  
  # Extract sample data (columns 7 onwards in .traw format)
  sample_data <- traw_data[, 7:ncol(traw_data)]
  
  # Get sample names and clean them up
  sample_names <- names(sample_data)
  
  # PLINK duplicates sample names as SAMPLE_SAMPLE
  # Check if splitting at the last underscore gives us identical parts
  sample_names_clean <- sapply(sample_names, function(name) {
    # Find the last underscore
    last_underscore <- regexpr("_[^_]*$", name)
    
    if (last_underscore > 0) {
      # Split at the last underscore
      first_part <- substr(name, 1, last_underscore - 1)
      second_part <- substr(name, last_underscore + 1, nchar(name))
      
      # If both parts are identical, it's PLINK's duplication pattern
      if (first_part == second_part) {
        return(first_part)
      }
    }
    
    # If no duplication pattern found, return original name
    return(name)
  }, USE.NAMES = FALSE)
  
  cat("Sample name cleaning examples:\n")
  for (i in 1:min(5, length(sample_names))) {
    if (sample_names[i] != sample_names_clean[i]) {
      cat("  ", sample_names[i], "→", sample_names_clean[i], "\n")
    } else {
      cat("  ", sample_names[i], "(no change - legitimate underscores)\n")
    }
  }
  
  # Use original marker IDs from the SNP column (column 2 in .traw)
  marker_names <- traw_data$SNP
  cat("Using original marker IDs\n")
  cat("Example marker IDs:", paste(head(marker_names, 3), collapse = ", "), "...\n")
  
  # Transpose so samples are rows, markers are columns (for find_parentage function)
  dosage_matrix <- t(sample_data)
  
  # Set proper names  
  colnames(dosage_matrix) <- marker_names
  rownames(dosage_matrix) <- sample_names_clean
  
  # Convert to data frame and add ID column for sample names
  final_matrix <- data.frame(
    id = sample_names_clean,  # Sample names in the id column
    dosage_matrix,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("Writing formatted dosage matrix:", output_file, "\n")
  
  # Write with tab separation
  write.table(final_matrix, file = output_file, sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
  
  cat("Dosage matrix created successfully!\n")
  cat("Format: id column + markers as columns\n")  
  cat("Samples as rows (", nrow(final_matrix), "), markers as columns (", ncol(final_matrix)-1, ")\n")
  
  return(output_file)
}

# Run formatting
result <- format_plink_to_dosage_matrix(opt$input, opt$output)
cat("Final dosage matrix:", result, "\n")