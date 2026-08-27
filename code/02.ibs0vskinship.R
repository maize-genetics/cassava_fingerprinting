# Analyze relationship data from PLINK and King for DArTseqLD data

# Load libraries ----
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gt)
library(webshot2)
# Set working directory ----
setwd("~/Documents/cassava_fingerprinting/")

# ==============================================================================
# 1. Read and process data
# ==============================================================================

# Read KING kinship data 
king_data <- read.table("output/plinkAndKing_geno0.2_mind0.2_maf0.01.kin0", 
                        header = TRUE, comment.char = "", stringsAsFactors = FALSE)
colnames(king_data) <- c("FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP")

# Read PLINK IBD data 
plink_data <- read.table("output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01.genome", 
                         header = TRUE, stringsAsFactors = FALSE)

# Rename X.FID1 column
colnames(king_data_check)[1] <- "FID1"

# Define the problematic patterns
fix_split_names <- function(data) {
  data %>%
    mutate(
      # Fix IID1 first
      IID1 = case_when(
        FID1 == "MASAKA" & IID1 == "LOCAL-2" ~ "MASAKA_LOCAL-2",
        FID1 == "NJULE" & IID1 == "RED" ~ "NJULE_RED",
        FID1 == "MASAKA" & IID1 == "LOCAL-1" ~ "MASAKA_LOCAL-1",
        FID1 == "OFUMBA" & IID1 == "CHAI" ~ "OFUMBA_CHAI",
        TRUE ~ IID1
      ),
      # Fix IID2
      IID2 = case_when(
        FID2 == "MASAKA" & IID2 == "LOCAL-2" ~ "MASAKA_LOCAL-2",
        FID2 == "NJULE" & IID2 == "RED" ~ "NJULE_RED",
        FID2 == "MASAKA" & IID2 == "LOCAL-1" ~ "MASAKA_LOCAL-1",
        FID2 == "OFUMBA" & IID2 == "CHAI" ~ "OFUMBA_CHAI",
        TRUE ~ IID2
      ),
      # Fix FID1 to match corrected IID1
      FID1 = case_when(
        IID1 == "MASAKA_LOCAL-2" ~ "MASAKA_LOCAL-2",
        IID1 == "NJULE_RED" ~ "NJULE_RED",
        IID1 == "MASAKA_LOCAL-1" ~ "MASAKA_LOCAL-1",
        IID1 == "OFUMBA_CHAI" ~ "OFUMBA_CHAI",
        TRUE ~ FID1
      ),
      # Fix FID2 to match corrected IID2
      FID2 = case_when(
        IID2 == "MASAKA_LOCAL-2" ~ "MASAKA_LOCAL-2",
        IID2 == "NJULE_RED" ~ "NJULE_RED",
        IID2 == "MASAKA_LOCAL-1" ~ "MASAKA_LOCAL-1",
        IID2 == "OFUMBA_CHAI" ~ "OFUMBA_CHAI",
        TRUE ~ FID2
      )
    )
}

# Apply fix to king_data
king_data <- fix_split_names(king_data)

# Apply fix to plink_data  
plink_data <- fix_split_names(plink_data)

# Verify the fix worked for both datasets
cat("=== VERIFICATION ===\n")
cat("KING data - checking all columns:\n")
king_verify <- king_data %>% 
  filter(IID1 %in% c("MASAKA_LOCAL-2", "NJULE_RED", "MASAKA_LOCAL-1", "OFUMBA_CHAI"))
print(head(king_verify, 4))

cat("\nPLINK data - checking all columns:\n")
plink_verify <- plink_data %>% 
  filter(IID1 %in% c("MASAKA_LOCAL-2", "NJULE_RED", "MASAKA_LOCAL-1", "OFUMBA_CHAI"))
print(head(plink_verify, 4))

# Define reference varieties
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

# Merge KING and PLINK data 
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
# 2. Classify relationships (can modify this)
# ==============================================================================

classify_relationships <- function(ibs0, kinship) {
  case_when(
    kinship >= 0.36 ~ "Highly related / clones",
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
# 3a. Get max relationship per reference x farm pair: MAX_REF_RELATIONSHIPS.CSV 
# ==============================================================================

# Get only farm samples and their strongest reference relationship
max_ref_relationships <- combined_data %>%
  filter(IID1 %in% reference_varieties | IID2 %in% reference_varieties) %>%
  # Filter out reference-reference pairs
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
# 3b. Get reference x reference data: REF_REF_RELATIONSHIPS.CSV
# ==============================================================================

# Get reference-reference relationships (all relationships, not just strongest or max) 
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

# Save reference-reference relationships
write.csv(ref_ref_relationships, "ref_ref_relationships.csv", row.names = FALSE)

cat("=== REF-REF RELATIONSHIPS SAVED ===\n")
cat("Total reference-reference relationships:", nrow(ref_ref_relationships), "\n")

# Summary of ref-ref relationships
cat("\nReference-Reference relationship types:\n")
print(table(ref_ref_relationships$Relationship))

cat("\nSample ref-ref relationships:\n")
print(head(ref_ref_relationships))

# ==============================================================================
# 4a. Create reference x farm max relationships plot (King coefficient x IBS0)
# ==============================================================================

# Set up colors
relationship_colors <- c(
  "Highly related / clones" = "#FF0000",
  "First degree" = "#0066CC", 
  "Second degree" = "#FF69B4",
  "Third degree" = "#8A2BE2",
  "Fourth degree and unrelated" = "#CCCCCC"
)

# Create plot
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
    legend.position = "right"
  ) +
  coord_cartesian(xlim = c(0, 0.52), ylim = c(0, 0.1))

# Display and save plot
print(p_max_ref)
ggsave("maximum_reference_relationships_plot.png", p_max_ref, 
       width = 12, height = 8, dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. max_ref_relationships.csv -", nrow(max_ref_relationships), "samples\n")
cat("2. maximum_reference_relationships_plot.png\n")

# ==============================================================================
# 4b. Create reference x reference figure (King coefficient x IBS0) 
# ==============================================================================

# Filter for strong relationships only in ref-ref relationships
strong_ref_ref_relationships <- ref_ref_relationships %>%
  filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree", "Third degree"))

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
  
  # Create the ref-ref plot
  p_strong_ref_ref <- ggplot(strong_ref_ref_relationships, aes(x = KINSHIP, y = IBS0, color = Relationship)) +
    geom_point(size = 3, alpha = 0.8) +
    
    # Use ggrepel for nicer non-overlapping labels
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
      legend.position = "right"
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
# 5. Make table for reference data
# ==============================================================================

# Create table for all max references with all relationship stats
complete_reference_stats <- data.frame(reference_variety = reference_varieties) %>%
  left_join(
    # Get vegetative clone counts
    max_ref_relationships %>%
      filter(Relationship == "Highly related / clones") %>%
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
        vegetative_clones = sum(Relationship == "Highly related / clones"),
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

# Create nicer summary table
enhanced_reference_table <- complete_reference_stats %>%
  mutate(
    # Add some extra stuff
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

# test nicer formatting?
if(!require(knitr)) install.packages("knitr")
if(!require(kableExtra)) install.packages("kableExtra")
library(knitr)
library(kableExtra)

# Create formatted table
kable(head(enhanced_reference_table, 15), 
      caption = "Top 15 Cassava Reference Varieties by Farm Usage",
      col.names = c("Reference Variety", "Highly related / clones", "Highly related / clones %", "Avg Kinship", 
                    "Strong Rels", "Total Connections", "1st°", "2nd°", "3rd°", "4th°", "Diversity"))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. max_ref_relationships.csv -", nrow(max_ref_relationships), "samples\n")
cat("2. maximum_reference_relationships_plot.png\n")
cat("3. complete_reference_stats.csv - Reference variety analysis\n")
cat("4. Enhanced reference table displayed above\n")

# ==============================================================================
# Look at unmatched farm samples (no clonal/1st/2nd degree relationships with any reference)
# Then find their TRUE strongest connection (to ANY sample, farm or reference)
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 1: Identify farm samples that DO have a strong reference match
# ------------------------------------------------------------------------------
connected_farms_via_reference <- combined_data %>%
  filter((IID1 %in% reference_varieties & !IID2 %in% reference_varieties) | 
           (!IID1 %in% reference_varieties & IID2 %in% reference_varieties)) %>%
  filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree")) %>%
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

# ------------------------------------------------------------------------------
# Step 2: Get ALL farm samples in the dataset
# ------------------------------------------------------------------------------
all_farm_samples <- combined_data %>%
  filter(!IID1 %in% reference_varieties | !IID2 %in% reference_varieties) %>%
  {unique(c(.[!.$IID1 %in% reference_varieties, "IID1"], 
            .[!.$IID2 %in% reference_varieties, "IID2"]))} %>%
  .[!is.na(.)]

# ------------------------------------------------------------------------------
# Step 3: Farms WITHOUT a strong reference match
# ------------------------------------------------------------------------------
farms_without_ref_match <- setdiff(all_farm_samples, connected_farms_via_reference)

cat("=== FARM SAMPLES: REFERENCE MATCH STATUS ===\n")
cat("Total farm samples:", length(all_farm_samples), "\n")
cat("Farms WITH strong reference match:", length(connected_farms_via_reference), "\n") 
cat("Farms WITHOUT strong reference match:", length(farms_without_ref_match), "\n\n")

# ------------------------------------------------------------------------------
# Step 4: For farms lacking a reference match, find TRUE strongest connection
# (search FULL combined_data - partner can be ANY sample, farm or reference)
# ------------------------------------------------------------------------------
true_max_for_unmatched_farms <- combined_data %>%
  filter(IID1 %in% farms_without_ref_match | IID2 %in% farms_without_ref_match) %>%
  mutate(
    farm_sample = case_when(
      IID1 %in% farms_without_ref_match ~ IID1,
      IID2 %in% farms_without_ref_match ~ IID2
    ),
    partner = case_when(
      IID1 %in% farms_without_ref_match ~ IID2,
      IID2 %in% farms_without_ref_match ~ IID1
    ),
    partner_type = ifelse(partner %in% reference_varieties, "Reference", "Farm")
  ) %>%
  group_by(farm_sample) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(farm_sample, partner, partner_type, KINSHIP, IBS0, Relationship, NSNP) %>%
  mutate(
    relationship_category = case_when(
      partner_type == "Reference" ~ "Farm-Reference",
      TRUE ~ "Farm-Farm"
    )
  ) %>%
  arrange(desc(KINSHIP))

# Save full results
write.csv(true_max_for_unmatched_farms, "true_max_unmatched_farms.csv", row.names = FALSE)

cat("=== TRUE STRONGEST CONNECTION FOR FARMS LACKING REFERENCE MATCH ===\n")
cat("Total farms analyzed:", nrow(true_max_for_unmatched_farms), "\n\n")

# ------------------------------------------------------------------------------
# Step 5: Split into two meaningful groups
# ------------------------------------------------------------------------------

# Group A: Farms with NO reference match, but STRONG match to another farm
strong_farm_farm_connections <- true_max_for_unmatched_farms %>%
  filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree"),
         partner_type == "Farm")

# Group B: Truly isolated farms - no strong match to reference OR farm
truly_isolated_farms <- true_max_for_unmatched_farms %>%
  filter(!Relationship %in% c("Highly related / clones", "First degree", "Second degree"))

cat("=== SUMMARY ===\n")
cat("Farms with NO reference match, but STRONG match to another farm:", 
    nrow(strong_farm_farm_connections), "\n")
cat("Farms truly isolated (no strong match to reference OR farm):", 
    nrow(truly_isolated_farms), "\n\n")

# Save both groups separately
write.csv(strong_farm_farm_connections, "farms_with_farm_farm_clones.csv", row.names = FALSE)
write.csv(truly_isolated_farms, "truly_isolated_farms.csv", row.names = FALSE)

cat("Relationship breakdown for farm-farm connections:\n")
print(table(strong_farm_farm_connections$Relationship))

cat("\nRelationship breakdown for truly isolated farms:\n")
print(table(truly_isolated_farms$Relationship))

# ==============================================================================
# Step 6: Also check reference varieties for the same issue (isolated references)
# ==============================================================================

connected_references <- combined_data %>%
  filter(IID1 %in% reference_varieties | IID2 %in% reference_varieties) %>%
  filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree")) %>%
  {unique(c(.[.$IID1 %in% reference_varieties, "IID1"],
            .[.$IID2 %in% reference_varieties, "IID2"]))} %>%
  .[!is.na(.)]

unconnected_references <- setdiff(reference_varieties, connected_references)

cat("\n=== REFERENCE VARIETIES: CONNECTION STATUS ===\n")
cat("Total reference varieties:", length(reference_varieties), "\n")
cat("References WITH strong connection:", length(connected_references), "\n")
cat("References WITHOUT strong connection:", length(unconnected_references), "\n")

# Find TRUE strongest match for isolated references (search full dataset)
true_max_isolated_refs <- combined_data %>%
  filter(IID1 %in% unconnected_references | IID2 %in% unconnected_references) %>%
  mutate(
    reference = case_when(
      IID1 %in% unconnected_references ~ IID1,
      IID2 %in% unconnected_references ~ IID2
    ),
    partner = case_when(
      IID1 %in% unconnected_references ~ IID2,
      IID2 %in% unconnected_references ~ IID1
    ),
    partner_type = ifelse(partner %in% reference_varieties, "Reference", "Farm")
  ) %>%
  group_by(reference) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(reference, partner, partner_type, KINSHIP, IBS0, Relationship, NSNP) %>%
  mutate(
    relationship_category = case_when(
      partner_type == "Reference" ~ "Reference-Reference",
      TRUE ~ "Reference-Farm"
    )
  ) %>%
  arrange(desc(KINSHIP))

write.csv(true_max_isolated_refs, "true_max_isolated_references.csv", row.names = FALSE)

cat("\n=== TRUE STRONGEST MATCH FOR ISOLATED REFERENCE VARIETIES ===\n")
print(true_max_isolated_refs, n = Inf)

# ==============================================================================
# Step 7: Plots
# ==============================================================================

relationship_colors <- c(
  "Highly related / clones" = "#FF0000",
  "First degree" = "#0066CC", 
  "Second degree" = "#FF69B4",
  "Third degree" = "#8A2BE2",
  "Fourth degree and unrelated" = "#CCCCCC"
)

relationship_shapes <- c(
  "Farm-Farm" = 15,              # Square
  "Reference-Reference" = 16,    # Circle  
  "Farm-Reference" = 17          # Triangle
)

# Combine both groups - relationship_category already exists from earlier code
# (Farm-Farm for strong_farm_farm_connections, Farm-Reference or Farm-Farm for truly_isolated_farms)
combined_farm_analysis <- bind_rows(
  strong_farm_farm_connections,
  truly_isolated_farms
)

p_combined_farms <- ggplot(combined_farm_analysis, 
                           aes(x = KINSHIP, y = IBS0, 
                               color = Relationship, 
                               shape = relationship_category)) +
  geom_point(size = 2, alpha = 0.7) +
  
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = 0.05, 
           label = c("3rd", "2nd", "1st", "Clone"), 
           color = "gray30", size = 3, angle = 90, vjust = -0.5) +
  
  scale_color_manual(values = relationship_colors, name = "Relationship Type") +
  scale_shape_manual(values = relationship_shapes, name = "Pairing Type") +
  
  labs(
    title = "Farm Samples Lacking Reference Match: True Strongest Connection",
    subtitle = paste("Strongest connection for", nrow(combined_farm_analysis), 
                     "farm samples without clonal/1st/2nd degree reference match"),
    x = "KING Kinship Coefficient (Φ)",
    y = "IBS0 Coefficient"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  ) +
  coord_cartesian(xlim = c(0, 0.52), ylim = c(0, 0.1))

print(p_combined_farms)

ggsave("combined_farm_connectivity_plot.png", p_combined_farms, 
       width = 12, height = 8, dpi = 300)

cat("Saved combined plot: combined_farm_connectivity_plot.png\n")

# ==============================================================================
# Step 8: Final Summary
# ==============================================================================

cat("\n=== FINAL SUMMARY: FARM SAMPLE CONNECTIVITY ===\n")
cat("Total farm samples:", length(all_farm_samples), "\n")
cat("  - Connected to reference (clone/1st/2nd degree):", length(connected_farms_via_reference), 
    " (", round(length(connected_farms_via_reference)/length(all_farm_samples)*100, 1), "%)\n")
cat("  - No reference match, but connected to another farm:", nrow(strong_farm_farm_connections),
    " (", round(nrow(strong_farm_farm_connections)/length(all_farm_samples)*100, 1), "%)\n")
cat("  - Truly isolated (no strong match to anything):", nrow(truly_isolated_farms),
    " (", round(nrow(truly_isolated_farms)/length(all_farm_samples)*100, 1), "%)\n")

cat("\n=== FINAL SUMMARY: REFERENCE VARIETY CONNECTIVITY ===\n")
cat("Total reference varieties:", length(reference_varieties), "\n")
cat("  - Connected (clone/1st/2nd degree to anything):", length(connected_references),
    " (", round(length(connected_references)/length(reference_varieties)*100, 1), "%)\n")
cat("  - Isolated (no strong connections):", length(unconnected_references),
    " (", round(length(unconnected_references)/length(reference_varieties)*100, 1), "%)\n")

cat("\n=== FILES SAVED ===\n")
cat("1. true_max_unmatched_farms.csv - all farms lacking reference match with their true best partner\n")
cat("2. farms_with_farm_farm_clones.csv -", nrow(strong_farm_farm_connections), "farms with farm-farm strong matches\n")
cat("3. truly_isolated_farms.csv -", nrow(truly_isolated_farms), "farms with no strong match anywhere\n")
cat("4. true_max_isolated_references.csv -", nrow(true_max_isolated_refs), "isolated reference varieties with true best match\n")
cat("5. farm_farm_clones_no_reference_match.png - plot\n")
cat("6. truly_isolated_farms.png - plot\n")

cat("\n=== ANALYSIS COMPLETE ===\n")


Copy code
# ==============================================================================
# Reference varieties without strong matches - TABLE ONLY
# ==============================================================================

# Step 1: Identify reference varieties WITH a strong connection to anything
connected_references <- combined_data %>%
  filter(IID1 %in% reference_varieties | IID2 %in% reference_varieties) %>%
  filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree")) %>%
  {unique(c(.[.$IID1 %in% reference_varieties, "IID1"],
            .[.$IID2 %in% reference_varieties, "IID2"]))} %>%
  .[!is.na(.)]

unconnected_references <- setdiff(reference_varieties, connected_references)

cat("Total reference varieties:", length(reference_varieties), "\n")
cat("References WITHOUT strong connection:", length(unconnected_references), "\n\n")

# Step 2: For isolated references, find TRUE strongest match (full dataset search)
true_max_isolated_refs <- combined_data %>%
  filter(IID1 %in% unconnected_references | IID2 %in% unconnected_references) %>%
  mutate(
    reference = case_when(
      IID1 %in% unconnected_references ~ IID1,
      IID2 %in% unconnected_references ~ IID2
    ),
    partner = case_when(
      IID1 %in% unconnected_references ~ IID2,
      IID2 %in% unconnected_references ~ IID1
    ),
    partner_type = ifelse(partner %in% reference_varieties, "Reference", "Farm")
  ) %>%
  group_by(reference) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    Reference_Variety = reference,
    Strongest_Match = partner,
    Match_Type = partner_type,
    Kinship = KINSHIP,
    IBS0,
    Relationship,
    NSNP
  ) %>%
  mutate(
    Kinship = round(Kinship, 4),
    IBS0 = round(IBS0, 4)
  ) %>%
  arrange(desc(Kinship))

cat("=== REFERENCE VARIETIES WITHOUT STRONG MATCH: TRUE STRONGEST CONNECTION ===\n")
print(true_max_isolated_refs, n = Inf)

write.csv(true_max_isolated_refs, "isolated_reference_varieties_true_max.csv", row.names = FALSE)

gt_table <- true_max_isolated_refs %>%
  gt() %>%
  tab_header(
    title = "Isolated Reference Varieties: Strongest Genetic Match",
    subtitle = paste(nrow(true_max_isolated_refs), 
                     "reference varieties lacking clone/1st/2nd degree relationships")
  ) %>%
  fmt_number(columns = c(Kinship, IBS0), decimals = 4) %>%
  data_color(
    columns = Kinship,
    colors = scales::col_numeric(
      palette = c("red", "white", "green"),
      domain = c(min(true_max_isolated_refs$Kinship), 0.088)
    )
  )

# Save as PNG - drag this into Google Slides
gtsave(gt_table, "isolated_reference_varieties_table.png")