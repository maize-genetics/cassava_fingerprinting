# Load libraries ----
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set working directory ----
setwd("~/Documents/cassava_fingerprinting/")

# ==============================================================================
# 1. READ AND PROCESS DATA
# ==============================================================================

# Read KING kinship data ----
king_data <- read.table("output/plinkAndKing_geno0.2_mind0.2_maf0.01.kin0", 
                        header = TRUE, comment.char = "", stringsAsFactors = FALSE)
colnames(king_data) <- c("FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP")

# Read PLINK IBD data ----
plink_data <- read.table("output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01.genome", 
                         header = TRUE, stringsAsFactors = FALSE)

# Define reference varieties ----
reference_varieties <- c(
  "Nase1", "NASE11", "UG120024", "IITA-TMS-MM960608", "KABWA", "MUWOGO-MUMYUFU", 
  "MERCURY", "NASE2", "NASE12", "UG120183", "IITA-TMS-IBA120067", "UG110310", 
  "MASAKA_LOCAL-2", "NJULE-WHITE", "Nase3", "Nase13", "UG120156", "UG110052", 
  "BALICol1998", "NJULE_RED", "MASAKA_LOCAL-1", "Nase4", "NASE14", "UG120193", 
  "UG110309", "KITENGA", "NYARABOKE", "Nase5", "Nase16", "UG110164", "UG110114", 
  "MAGANA", "MUREFU", "NASE6", "Nase19", "Mkumba", "UG110304", "UG110306", 
  "MACHUNDE", "NASE8", "Narocas1", "BALICol2021", "KWATAMUMPALE", "BAO", 
  "LYAHOROLE", "Nase9", "NAROCASS2", "TMEB14", "EDYAL", "BUKALASA-11", "OFUMBA_CHAI"
)

# Merge KING and PLINK data ----
king_data$ID_pair <- paste(pmin(king_data$IID1, king_data$IID2), 
                           pmax(king_data$IID1, king_data$IID2), sep = "_")

plink_data$ID_pair <- paste(pmin(plink_data$IID1, plink_data$IID2), 
                            pmax(plink_data$IID1, plink_data$IID2), sep = "_")

combined_data <- merge(
  king_data[, c("ID_pair", "IID1", "IID2", "KINSHIP", "NSNP", "IBS0")], 
  plink_data[, c("ID_pair", "Z0", "Z1", "Z2", "PI_HAT")], 
  by = "ID_pair"
)

# ==============================================================================
# 2. CLASSIFY RELATIONSHIPS
# ==============================================================================

classify_relationships <- function(ibs0, kinship) {
  case_when(
    kinship >= 0.36 ~ "Vegetative clones",
    kinship >= 0.19 & kinship < 0.36 ~ "First degree",
    kinship >= 0.088 & kinship < 0.19 ~ "Second degree",
    kinship >= 0.044 & kinship < 0.088 ~ "Third degree",
    kinship < 0.044 ~ "Fourth degree and unrelated",
    TRUE ~ "Unclassified"
  )
}

combined_data$Relationship <- classify_relationships(
  ibs0 = combined_data$IBS0,
  kinship = combined_data$KINSHIP
)

# ==============================================================================
# 3a. CREATE MAX_REF_RELATIONSHIPS.CSV (Farm x Reference only)
# ==============================================================================

# Get ONLY farm samples and their strongest reference relationship ----
max_ref_relationships <- combined_data %>%
  filter(IID1 %in% reference_varieties | IID2 %in% reference_varieties) %>%
  # Filter OUT reference-reference pairs
  filter(!(IID1 %in% reference_varieties & IID2 %in% reference_varieties)) %>%
  mutate(
    farm_sample = case_when(
      IID1 %in% reference_varieties ~ IID2,
      IID2 %in% reference_varieties ~ IID1,
      TRUE ~ NA_character_
    ),
    ref_partner = case_when(
      IID1 %in% reference_varieties ~ IID1,
      IID2 %in% reference_varieties ~ IID2,  
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(farm_sample)) %>%
  group_by(farm_sample) %>%
  slice_min(abs(KINSHIP - 0.5), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample = farm_sample, ref_partner, KINSHIP, IBS0, Relationship, NSNP) %>%
  mutate(relationship_category = "Farm-Reference")

# Save ----
write.csv(max_ref_relationships, "max_ref_relationships.csv", row.names = FALSE)

# ==============================================================================
# 3b. CREATE REF_REF_RELATIONSHIPS.CSV
# ==============================================================================

# Get reference-reference relationships (ALL relationships, not just strongest) ----
ref_ref_relationships <- combined_data %>%
  # Only reference-reference pairs
  filter(IID1 %in% reference_varieties & IID2 %in% reference_varieties) %>%
  # Remove self-comparisons (a reference to itself)
  filter(IID1 != IID2) %>%
  select(
    sample = IID1, 
    ref_partner = IID2, 
    KINSHIP, 
    IBS0, 
    Relationship, 
    NSNP
  ) %>%
  mutate(relationship_category = "Reference-Reference")

# Save reference-reference relationships ----
write.csv(ref_ref_relationships, "ref_ref_relationships.csv", row.names = FALSE)

cat("=== REF-REF RELATIONSHIPS SAVED ===\n")
cat("Total reference-reference relationships:", nrow(ref_ref_relationships), "\n")

# Summary of ref-ref relationships
cat("\nReference-Reference relationship types:\n")
print(table(ref_ref_relationships$Relationship))

cat("\nSample ref-ref relationships:\n")
print(head(ref_ref_relationships))

# ==============================================================================
# 4a. CREATE MAXIMUM REFERENCE RELATIONSHIPS PLOT
# ==============================================================================

# Set up colors and shapes ----
relationship_colors <- c(
  "Vegetative clones" = "#FF0000",
  "First degree" = "#0066CC", 
  "Second degree" = "#FF69B4",
  "Third degree" = "#8A2BE2",
  "Fourth degree and unrelated" = "#CCCCCC"
)

# Create the plot ----
p_max_ref <- ggplot(max_ref_relationships, aes(x = KINSHIP, y = IBS0, 
                                               color = Relationship, 
                                               shape = relationship_category)) +
  geom_point(size = 2, alpha = 0.7) +
  
  # KING threshold lines
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = 0.05, 
           label = c("3rd", "2nd", "1st", "Clone"), 
           color = "gray30", size = 3, angle = 90, vjust = -0.5) +
  
  scale_color_manual(values = relationship_colors, name = "Relationship Type") +
  scale_shape_manual(values = c("Farm-Reference" = 17, 
                                "Reference-Reference" = 16),
                     name = "Pairing Type") +
  
  labs(
    title = "Maximum reference x farm relationships",
    subtitle = paste("Strongest reference connection for", nrow(max_ref_relationships), "samples"),
    x = "KING Kinship Coefficient (Φ)",
    y = "IBS0 Coefficient"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(0, 0.52), ylim = c(0, 0.1))

# Display and save plot ----
print(p_max_ref)
ggsave("maximum_reference_relationships_plot.png", p_max_ref, 
       width = 12, height = 8, dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. max_ref_relationships.csv -", nrow(max_ref_relationships), "samples\n")
cat("2. maximum_reference_relationships_plot.png\n")

# ==============================================================================
# 4b. Reference x reference figure (King coefficient x IBS0) 
# ==============================================================================

# Filter for strong relationships only in ref-ref relationships
strong_ref_ref_relationships <- ref_ref_relationships %>%
  filter(Relationship %in% c("Vegetative clones", "First degree", "Second degree", "Third degree"))

cat("=== STRONG REF-REF ANALYSIS ===\n")
cat("Total strong ref-ref relationships found:", nrow(strong_ref_ref_relationships), "\n")

# Show breakdown by relationship type
if(nrow(strong_ref_ref_relationships) > 0) {
  cat("\nStrong relationship types between references:\n")
  print(table(strong_ref_ref_relationships$Relationship))
  
  cat("\nTop 10 strongest ref-ref relationships:\n")
  print(strong_ref_ref_relationships %>% 
          arrange(desc(KINSHIP)) %>% 
          select(sample, ref_partner, KINSHIP, Relationship) %>%
          head(10))
  
  # Create the strong ref-ref plot
  p_strong_ref_ref <- ggplot(strong_ref_ref_relationships, aes(x = KINSHIP, y = IBS0, color = Relationship)) +
    geom_point(size = 3, alpha = 0.8) +
    
    # Use ggrepel for non-overlapping labels
    geom_text_repel(
      aes(label = paste0(sample, " → ", ref_partner)),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      max.overlaps = 20,
      force = 2,
      show.legend = FALSE
    ) +
    
    # KING threshold lines
    geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
               linetype = "dashed", alpha = 0.6, color = "gray30") +
    annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = 0.045, 
             label = c("3rd", "2nd", "1st", "Clone"), 
             color = "gray30", size = 4, angle = 90, vjust = -0.5) +
    
    scale_color_manual(values = relationship_colors, name = "Relationship Type") +
    
    labs(
      title = "Strongest relationships between reference varieties",
      subtitle = paste0("Strongest relationships among reference varieties"),
      x = "KING Kinship Coefficient (Φ)",
      y = "IBS0 Coefficient"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    ) +
    coord_cartesian(xlim = c(0, 0.52), ylim = c(0, 0.1))
  
  # Display and save
  print(p_strong_ref_ref)
  ggsave("reference_strong_relationships.png", p_strong_ref_ref, 
         width = 14, height = 10, dpi = 300)
  
} else {
  cat("No strong relationships found between reference varieties.\n")
  
  # Show the overall distribution instead
  cat("\nOverall ref-ref relationship distribution:\n")
  print(table(ref_ref_relationships$Relationship))
}

# ==============================================================================
# 5. REFERENCE VARIETY ANALYSIS TABLE
# ==============================================================================

# Create complete table for ALL references with ALL relationship stats
complete_reference_stats <- data.frame(reference_variety = reference_varieties) %>%
  left_join(
    # Get vegetative clone counts
    max_ref_relationships %>%
      filter(Relationship == "Vegetative clones") %>%
      group_by(ref_partner) %>%
      summarise(
        clone_count = n(),
        clone_avg_kinship = round(mean(KINSHIP, na.rm = TRUE), 3),
        .groups = 'drop'
      ),
    by = c("reference_variety" = "ref_partner")
  ) %>%
  left_join(
    # Get ALL relationship stats (including clones, 1st, 2nd, 3rd degree, etc.)
    max_ref_relationships %>%
      group_by(ref_partner) %>%
      summarise(
        total_connections = n(),
        overall_avg_kinship = round(mean(KINSHIP, na.rm = TRUE), 3),
        overall_min_kinship = round(min(KINSHIP, na.rm = TRUE), 3),
        overall_max_kinship = round(max(KINSHIP, na.rm = TRUE), 3),
        # Count each relationship type
        first_degree = sum(Relationship == "First degree"),
        second_degree = sum(Relationship == "Second degree"),
        third_degree = sum(Relationship == "Third degree"),
        vegetative_clones = sum(Relationship == "Vegetative clones"),
        fourth_degree = sum(Relationship == "Fourth degree and unrelated"),
        unrelated = sum(Relationship == "unrelated"),
        .groups = 'drop'
      ),
    by = c("reference_variety" = "ref_partner")
  ) %>%
  # Clean up NAs
  mutate(
    number_of_clones = ifelse(is.na(clone_count), 0, clone_count),
    clone_avg_kinship = ifelse(is.na(clone_avg_kinship), 0, clone_avg_kinship),
    total_connections = ifelse(is.na(total_connections), 0, total_connections),
    overall_avg_kinship = ifelse(is.na(overall_avg_kinship), 0, overall_avg_kinship),
    overall_min_kinship = ifelse(is.na(overall_min_kinship), 0, overall_min_kinship),
    overall_max_kinship = ifelse(is.na(overall_max_kinship), 0, overall_max_kinship),
    # Replace NAs in relationship counts with 0
    across(first_degree:unrelated, ~ifelse(is.na(.), 0, .))
  ) %>%
  # Reorder columns for clarity
  select(reference_variety, number_of_clones, clone_avg_kinship, total_connections, 
         overall_avg_kinship, overall_min_kinship, overall_max_kinship,
         first_degree, second_degree, third_degree, fourth_degree, unrelated) %>%
  arrange(desc(number_of_clones), desc(total_connections))

# Display the complete table
print(complete_reference_stats)
write.csv(complete_reference_stats, "complete_reference_stats.csv", row.names = FALSE)

# Create an enhanced summary table
enhanced_reference_table <- complete_reference_stats %>%
  mutate(
    # Add some calculated fields for better insights
    clone_percentage = round((number_of_clones / sum(number_of_clones)) * 100, 1),
    has_strong_relationships = first_degree + second_degree + number_of_clones,
    relationship_diversity = (first_degree > 0) + (second_degree > 0) + (third_degree > 0) + (number_of_clones > 0)
  ) %>%
  select(reference_variety, number_of_clones, clone_percentage, 
         overall_avg_kinship, has_strong_relationships, total_connections,
         first_degree, second_degree, third_degree, fourth_degree, 
         relationship_diversity) %>%
  arrange(desc(number_of_clones), desc(has_strong_relationships))

# Display enhanced table
print(enhanced_reference_table)

# Try to load kable packages for nice formatting
if(!require(knitr)) install.packages("knitr")
if(!require(kableExtra)) install.packages("kableExtra")
library(knitr)
library(kableExtra)

# Create formatted table
kable(head(enhanced_reference_table, 15), 
      caption = "Top 15 Cassava Reference Varieties by Farm Usage",
      col.names = c("Reference Variety", "Clones", "Clone %", "Avg Kinship", 
                    "Strong Rels", "Total Connections", "1st°", "2nd°", "3rd°", "4th°", "Diversity"))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. max_ref_relationships.csv -", nrow(max_ref_relationships), "samples\n")
cat("2. maximum_reference_relationships_plot.png\n")
cat("3. complete_reference_stats.csv - Reference variety analysis\n")
cat("4. Enhanced reference table displayed above\n")

# ==============================================================================
# UNCONNECTED SAMPLES RELATIONSHIPS (3rd, 4th, unrelated)
# ==============================================================================

# Step 1: Identify "unconnected" samples
# Find farm samples with NO clonal/1st/2nd degree relationships with ANY reference
unconnected_farm_samples <- combined_data %>%
  # Get all farm-reference relationships
  filter((IID1 %in% reference_varieties & !IID2 %in% reference_varieties) | 
           (!IID1 %in% reference_varieties & IID2 %in% reference_varieties)) %>%
  # Filter for strong relationships only
  filter(Relationship %in% c("Vegetative clones", "First degree", "Second degree")) %>%
  # Get the farm samples that DO have strong ref relationships
  mutate(
    farm_sample = case_when(
      IID1 %in% reference_varieties ~ IID2,
      IID2 %in% reference_varieties ~ IID1,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(farm_sample)) %>%
  pull(farm_sample) %>%
  unique()

# All farm samples that are NOT in the connected list
all_farm_samples <- combined_data %>%
  filter(!IID1 %in% reference_varieties | !IID2 %in% reference_varieties) %>%
  {unique(c(.[!.$IID1 %in% reference_varieties, "IID1"], 
            .[!.$IID2 %in% reference_varieties, "IID2"]))} %>%
  .[!is.na(.)]

unconnected_farms <- setdiff(all_farm_samples, unconnected_farm_samples)

cat("Total farm samples:", length(all_farm_samples), "\n")
cat("Farm samples connected to references:", length(unconnected_farm_samples), "\n") 
cat("Unconnected farm samples:", length(unconnected_farms), "\n")

# Find reference varieties with NO clonal/1st/2nd degree relationships with ANY other sample
connected_references <- combined_data %>%
  filter(IID1 %in% reference_varieties | IID2 %in% reference_varieties) %>%
  filter(Relationship %in% c("Vegetative clones", "First degree", "Second degree")) %>%
  {unique(c(.[.$IID1 %in% reference_varieties, "IID1"],
            .[.$IID2 %in% reference_varieties, "IID2"]))} %>%
  .[!is.na(.)]

unconnected_references <- setdiff(reference_varieties, connected_references)

cat("Total reference varieties:", length(reference_varieties), "\n")
cat("References with strong connections:", length(connected_references), "\n")
cat("Unconnected references:", length(unconnected_references), "\n")

# Combine all unconnected samples
all_unconnected_samples <- c(unconnected_farms, unconnected_references)
cat("Total unconnected samples:", length(all_unconnected_samples), "\n")

# Step 2: Create third CSV with ALL relationships among unconnected samples
unconnected_relationships <- combined_data %>%
  # Keep only relationships where BOTH samples are in the unconnected list
  filter(IID1 %in% all_unconnected_samples & IID2 %in% all_unconnected_samples) %>%
  # Remove self-comparisons
  filter(IID1 != IID2) %>%
  select(
    sample = IID1,
    partner = IID2,
    KINSHIP,
    IBS0,
    Relationship,
    NSNP
  ) %>%
  # Add relationship category
  mutate(
    sample_type = ifelse(sample %in% reference_varieties, "Reference", "Farm"),
    partner_type = ifelse(partner %in% reference_varieties, "Reference", "Farm"),
    relationship_category = case_when(
      sample_type == "Farm" & partner_type == "Farm" ~ "Farm-Farm",
      sample_type == "Reference" & partner_type == "Reference" ~ "Reference-Reference", 
      TRUE ~ "Farm-Reference"
    )
  )

# Save the unconnected relationships
write.csv(unconnected_relationships, "unconnected_relationships.csv", row.names = FALSE)

# Print summary
cat("\n=== UNCONNECTED RELATIONSHIPS SUMMARY ===\n")
cat("Total unconnected relationships:", nrow(unconnected_relationships), "\n")
cat("Relationship categories:\n")
print(table(unconnected_relationships$relationship_category))
cat("\nRelationship types:\n")
print(table(unconnected_relationships$Relationship))

# Show some examples
cat("\n=== SAMPLE UNCONNECTED RELATIONSHIPS ===\n")
print(head(unconnected_relationships %>% arrange(desc(KINSHIP)), 10))

# Print summary with min/max kinship
cat("\n=== UNCONNECTED RELATIONSHIPS SUMMARY ===\n")
cat("Total unconnected relationships:", nrow(unconnected_relationships), "\n")

# Overall kinship statistics
cat("\nOverall kinship statistics:\n")
cat("Min kinship:", round(min(unconnected_relationships$KINSHIP), 4), "\n")
cat("Max kinship:", round(max(unconnected_relationships$KINSHIP), 4), "\n")
cat("Mean kinship:", round(mean(unconnected_relationships$KINSHIP), 4), "\n")

# Kinship statistics by relationship category
cat("\nKinship statistics by relationship category:\n")
kinship_by_category <- unconnected_relationships %>%
  group_by(relationship_category) %>%
  summarise(
    count = n(),
    min_kinship = round(min(KINSHIP), 4),
    max_kinship = round(max(KINSHIP), 4),
    mean_kinship = round(mean(KINSHIP), 4),
    .groups = 'drop'
  )
print(kinship_by_category)

cat("\nRelationship categories:\n")
print(table(unconnected_relationships$relationship_category))

cat("\nRelationship types:\n")
print(table(unconnected_relationships$Relationship))

# Show strongest relationships among unconnected samples
cat("\n=== STRONGEST UNCONNECTED RELATIONSHIPS ===\n")
print(head(unconnected_relationships %>% arrange(desc(KINSHIP)), 10))

# ==============================================================================
# GGPLOT FOR MAX UNCONNECTED RELATIONSHIPS  
# ==============================================================================

library(ggplot2)

# Create MAX unconnected relationships (strongest per sample)
max_unconnected_relationships <- unconnected_relationships %>%
  group_by(sample) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup()

# Set up colors for relationship types
relationship_colors <- c(
  "Vegetative clones" = "#FF0000",
  "First degree" = "#0066CC", 
  "Second degree" = "#FF69B4",
  "Third degree" = "#8A2BE2",
  "Fourth degree and unrelated" = "#CCCCCC"
)

# Set up shapes for relationship categories
relationship_shapes <- c(
  "Farm-Farm" = 15,              # Square
  "Reference-Reference" = 16,    # Circle  
  "Farm-Reference" = 17          # Triangle
)

# Create the plot
p_max_unconnected <- ggplot(max_unconnected_relationships, aes(x = KINSHIP, y = IBS0, 
                                                               color = Relationship, 
                                                               shape = relationship_category)) +
  geom_point(size = 2, alpha = 0.7) +
  
  # KING threshold lines
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = 0.05, 
           label = c("3rd", "2nd", "1st", "Clone"), 
           color = "gray30", size = 3, angle = 90, vjust = -0.5) +
  
  scale_color_manual(values = relationship_colors, name = "Relationship Type") +
  scale_shape_manual(values = relationship_shapes, name = "Pairing Type") +
  
  labs(
    title = "Maximum relationships among unmatched farm samples",
    subtitle = paste("Strongest connections for", nrow(max_unconnected_relationships), "farm samples without reference matches"),
    x = "KING Kinship Coefficient (Φ)",
    y = "IBS0 Coefficient"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom"
  ) +
  coord_cartesian(xlim = c(0, 0.52), ylim = c(0, 0.1))

# Display and save plot
print(p_max_unconnected)
ggsave("maximum_unconnected_relationships_plot.png", p_max_unconnected, 
       width = 12, height = 8, dpi = 300)

# Save the data
write.csv(max_unconnected_relationships, "max_unconnected_relationships.csv", row.names = FALSE)

cat("\n=== UNCONNECTED ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. max_unconnected_relationships.csv -", nrow(max_unconnected_relationships), "samples\n")
cat("2. maximum_unconnected_relationships_plot.png\n")

# Print summary stats
cat("\nUnconnected relationship summary:\n")
print(table(max_unconnected_relationships$relationship_category))
cat("\nKinship range:", round(min(max_unconnected_relationships$KINSHIP), 3), 
    "to", round(max(max_unconnected_relationships$KINSHIP), 3), "\n")
