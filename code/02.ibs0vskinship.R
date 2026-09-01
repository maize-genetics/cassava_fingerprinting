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
colnames(king_data)[1] <- "FID1"

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
write.csv(true_max_for_unmatched_farms, "max_unmatched_farms.csv", row.names = FALSE)

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


##### Look at inbreeding coefficient

het_data <- read.table("output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01.het", 
                       header = TRUE, stringsAsFactors = FALSE)

head(het_data)
colnames(het_data)

het_data <- het_data %>%
  mutate(
    IID = case_when(
      FID == "MASAKA" & IID == "LOCAL-2" ~ "MASAKA_LOCAL-2",
      FID == "NJULE" & IID == "RED" ~ "NJULE_RED",
      FID == "MASAKA" & IID == "LOCAL-1" ~ "MASAKA_LOCAL-1",
      FID == "OFUMBA" & IID == "CHAI" ~ "OFUMBA_CHAI",
      TRUE ~ IID
    ),
    FID = case_when(
      IID == "MASAKA_LOCAL-2" ~ "MASAKA_LOCAL-2",
      IID == "NJULE_RED" ~ "NJULE_RED",
      IID == "MASAKA_LOCAL-1" ~ "MASAKA_LOCAL-1",
      IID == "OFUMBA_CHAI" ~ "OFUMBA_CHAI",
      TRUE ~ FID
    )
  )

# Verify the fix worked
het_data %>% filter(IID %in% c("MASAKA_LOCAL-2", "NJULE_RED", "MASAKA_LOCAL-1", "OFUMBA_CHAI"))

# Plot kinship coefficient for max pair-wise farm x reference relationship by the inbreeding coefficient of the farm sample 

max_ref_relationships <- max_ref_relationships %>%
  left_join(
    het_data %>% select(IID, F),
    by = c("sample" = "IID")
  )

# Check the merge worked
head(max_ref_relationships)

# Check for any samples that didn't get an F value (missing from het_data)
cat("Samples missing F value:", sum(is.na(max_ref_relationships$F)), "\n")


p_max_ref_inbreeding <- ggplot(max_ref_relationships, aes(x = KINSHIP, y = F, 
                                                          color = Relationship, 
                                                          shape = relationship_category)) +
  geom_point(size = 2, alpha = 0.7) +
  
  # KING threshold lines
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = max(max_ref_relationships$F, na.rm = TRUE) * 0.95, 
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
    y = "Inbreeding Coefficient (F)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  ) +
  coord_cartesian(xlim = c(0, 0.52))

print(p_max_ref_inbreeding)

ggsave("maximum_reference_relationships_inbreeding_plot.png", p_max_ref_inbreeding, 
       width = 12, height = 8, dpi = 300)

# Check out doing the diffence in the inbreeding coefficient between the ref and the farm sample

max_ref_relationships <- max_ref_relationships %>%
  left_join(
    het_data %>% select(IID, F_ref = F),
    by = c("ref_partner" = "IID")
  )

# Rename the farm sample's F column for clarity (if not already done)
max_ref_relationships <- max_ref_relationships %>%
  rename(F_farm = F)

# Check the merge
head(max_ref_relationships)

# Check for missing values
cat("Missing F_ref values:", sum(is.na(max_ref_relationships$F_ref)), "\n")

# Calculate delta F
max_ref_relationships <- max_ref_relationships %>%
  mutate(delta_F = F_farm - F_ref)

# Filter to just First degree (blue) points
blue_only <- max_ref_relationships %>%
  filter(Relationship == "First degree")

# Correlation test: does lower kinship correlate with higher delta_F?
cor_test_deltaF <- cor.test(blue_only$KINSHIP, blue_only$delta_F, method = "pearson")
print(cor_test_deltaF)

# Plot
library(ggplot2)
p_blue_deltaF <- ggplot(blue_only, aes(x = KINSHIP, y = delta_F)) +
  geom_point(color = "#0066CC", size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    title = "First Degree Relationships: Kinship vs ΔF (Farm - Reference)",
    subtitle = paste0("Pearson r = ", round(cor_test_deltaF$estimate, 3), 
                      ", p = ", format.pval(cor_test_deltaF$p.value, digits = 3)),
    x = "KING Kinship Coefficient (Φ)",
    y = "ΔF (Farm Inbreeding − Reference Inbreeding)"
  ) +
  theme_minimal()

print(p_blue_deltaF)

p_all_deltaF <- ggplot(max_ref_relationships, aes(x = KINSHIP, y = delta_F, color = Relationship)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = max(max_ref_relationships$delta_F, na.rm = TRUE) * 0.95, 
           label = c("3rd", "2nd", "1st", "Clone"), 
           color = "gray30", size = 3, angle = 90, vjust = -0.5) +
  
  scale_color_manual(values = relationship_colors, name = "Relationship Type") +
  
  labs(
    title = "Kinship vs ΔF (Farm − Reference) Across All Relationship Types",
    subtitle = paste("All", nrow(max_ref_relationships), "farm-reference pairs"),
    x = "KING Kinship Coefficient (Φ)",
    y = "ΔF (Farm Inbreeding − Reference Inbreeding)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

print(p_all_deltaF)

ggsave("kinship_vs_deltaF_all_relationships.png", p_all_deltaF, width = 12, height = 8, dpi = 300)

# Boxplot F for ref versus farm

het_data_labeled <- het_data %>%
  mutate(sample_type = ifelse(IID %in% reference_varieties, "Reference", "Farm"))

# Summary stats
het_data_labeled %>%
  group_by(sample_type) %>%
  summarise(
    n = n(),
    mean_F = mean(F, na.rm = TRUE),
    median_F = median(F, na.rm = TRUE),
    sd_F = sd(F, na.rm = TRUE)
  )

# Statistical test
wilcox.test(F ~ sample_type, data = het_data_labeled)

# Visualize
library(ggplot2)
ggplot(het_data_labeled, aes(x = sample_type, y = F, fill = sample_type)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Inbreeding Coefficient: Reference Varieties vs Farm Samples",
       y = "Inbreeding Coefficient (F)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

# Same but by category

# ==============================================================================
# Check if negative F farm samples cluster within specific relationship categories
# ==============================================================================

# Summary: F_farm distribution by relationship category
F_by_relationship <- max_ref_relationships %>%
  group_by(Relationship) %>%
  summarise(
    n = n(),
    pct_negative_F = round(100 * mean(F_farm < 0, na.rm = TRUE), 1),
    mean_F_farm = round(mean(F_farm, na.rm = TRUE), 3),
    median_F_farm = round(median(F_farm, na.rm = TRUE), 3),
    min_F_farm = round(min(F_farm, na.rm = TRUE), 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(pct_negative_F))

print(F_by_relationship)

# Visualize: boxplot of F_farm split by relationship category
library(ggplot2)

# Reorder Relationship factor levels so red (clones) is furthest right
max_ref_relationships$Relationship <- factor(
  max_ref_relationships$Relationship,
  levels = c("Fourth degree and unrelated", "Third degree", "Second degree", 
             "First degree", "Highly related / clones")
)

# Recreate the plot with the new ordering
p_F_by_relationship <- ggplot(max_ref_relationships, aes(x = Relationship, y = F_farm, fill = Relationship)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = relationship_colors) +
  labs(
    title = "Farm Sample Inbreeding Coefficient (F) by Relationship to closest Reference",
    y = "Inbreeding Coefficient (F) - Farm Sample",
    x = "Relationship Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_F_by_relationship)

ggsave("F_farm_by_relationship_category.png", p_F_by_relationship, width = 10, height = 7, dpi = 300)

# ==============================================================================
# COMPREHENSIVE MULTI-PARENT / SPONTANEOUS CROSS CHECK
# Using ALL samples (farm + reference), NO nearest-neighbor pre-filtering
# ==============================================================================

library(dplyr)
library(igraph)

# ------------------------------------------------------------------------------
# Step 1: Build clone clusters using ALL samples (farm + reference)
# ------------------------------------------------------------------------------

all_clone_edges <- combined_data %>%
  filter(Relationship == "Highly related / clones") %>%
  filter(IID1 != IID2) %>%
  select(IID1, IID2)

clone_graph <- graph_from_data_frame(all_clone_edges, directed = FALSE)
clone_components <- components(clone_graph)

cluster_lookup <- data.frame(
  sample_name = names(clone_components$membership),
  cluster_id = clone_components$membership,
  stringsAsFactors = FALSE
)

# Add singleton samples (no clone relationships at all) as their own unique cluster
all_sample_names <- unique(c(combined_data$IID1, combined_data$IID2))
singleton_samples <- setdiff(all_sample_names, cluster_lookup$sample_name)

if (length(singleton_samples) > 0) {
  singleton_lookup <- data.frame(
    sample_name = singleton_samples,
    cluster_id = max(cluster_lookup$cluster_id) + seq_along(singleton_samples),
    stringsAsFactors = FALSE
  )
  cluster_lookup <- bind_rows(cluster_lookup, singleton_lookup)
}

# Create readable cluster labels, flagging whether each cluster contains reference(s)
cluster_labels <- cluster_lookup %>%
  mutate(is_ref = sample_name %in% reference_varieties) %>%
  group_by(cluster_id) %>%
  summarise(
    cluster_label = paste(sort(sample_name), collapse = " = "),
    n_members = n(),
    contains_reference = any(is_ref),
    .groups = 'drop'
  )

cluster_lookup <- cluster_lookup %>%
  left_join(cluster_labels, by = "cluster_id")

cat("=== CLONE CLUSTER SUMMARY (ALL SAMPLES) ===\n")
cat("Total samples:", length(all_sample_names), "\n")
cat("Total distinct clonal clusters:", n_distinct(cluster_lookup$cluster_id), "\n")
cat("Clusters with >1 member (true clone groups):", sum(cluster_labels$n_members > 1), "\n\n")

write.csv(cluster_lookup, "full_clone_cluster_lookup.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Step 2: For EVERY sample, pull ALL first-degree+ relationships (no pre-filtering)
# ------------------------------------------------------------------------------

all_strong_relationships <- combined_data %>%
  filter(IID1 != IID2) %>%
  filter(Relationship %in% c("Highly related / clones", "First degree")) %>%
  select(IID1, IID2, KINSHIP, IBS0, Relationship, NSNP)

# Make this symmetric (each sample sees all its partners, regardless of which column it was in)
symmetric_relationships <- bind_rows(
  all_strong_relationships %>% rename(sample = IID1, partner = IID2),
  all_strong_relationships %>% rename(sample = IID2, partner = IID1)
)

# Attach cluster info for the PARTNER (so we can count distinct clusters matched)
symmetric_relationships <- symmetric_relationships %>%
  left_join(cluster_lookup %>% select(sample_name, cluster_id, cluster_label, contains_reference),
            by = c("partner" = "sample_name"))

# ------------------------------------------------------------------------------
# Step 3: For every sample, count DISTINCT clusters it has a strong relationship to
#         EXCLUDING its own cluster (i.e., don't count relationships to its own clones)
# ------------------------------------------------------------------------------

# Get each sample's own cluster_id, so we can exclude self-cluster matches
own_cluster <- cluster_lookup %>% select(sample_name, own_cluster_id = cluster_id)

symmetric_relationships <- symmetric_relationships %>%
  left_join(own_cluster, by = c("sample" = "sample_name")) %>%
  filter(cluster_id != own_cluster_id)  # exclude matches to your own clone-cluster members

multi_parent_check_full <- symmetric_relationships %>%
  group_by(sample) %>%
  summarise(
    n_distinct_clusters = n_distinct(cluster_id),
    cluster_labels_matched = paste(unique(cluster_label), collapse = " ;; "),
    any_reference_involved = any(contains_reference),
    kinship_values = paste(round(KINSHIP, 3), collapse = ", "),
    relationship_types = paste(Relationship, collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_distinct_clusters))

cat("=== COMPREHENSIVE MULTI-PARENT SIGNATURE CHECK ===\n")
cat("(No nearest-neighbor filtering - checked ALL samples, ALL relationships)\n\n")

cat("Samples with strong relationships to 2+ distinct clonal lineages:\n")
multi_parent_candidates <- multi_parent_check_full %>% filter(n_distinct_clusters > 1)
print(as.data.frame(multi_parent_candidates))

cat("\nTotal candidates found:", nrow(multi_parent_candidates), "\n")
cat("Of these, how many involve at least one reference variety:", 
    sum(multi_parent_candidates$any_reference_involved), "\n")
cat("Of these, how many are PURELY farm-farm (no reference involved at all):", 
    sum(!multi_parent_candidates$any_reference_involved), "\n\n")

cat("Distribution of number of distinct clusters matched:\n")
cluster_count_dist <- multi_parent_check_full %>%
  count(n_distinct_clusters)
print(as.data.frame(cluster_count_dist))

write.csv(multi_parent_check_full, "comprehensive_multi_parent_check.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Step 4: For the candidates, check IBS0 to see if any show TRUE parent-offspring
#         signature (IBS0 near clone-floor ~0.002-0.003) vs sibling signature (elevated)
# ------------------------------------------------------------------------------

candidate_names <- multi_parent_candidates$sample
cat("\nNumber of candidate names:", length(candidate_names), "\n")

candidate_detail <- symmetric_relationships %>%
  filter(sample %in% candidate_names) %>%
  select(sample, partner, cluster_label, KINSHIP, IBS0, Relationship) %>%
  arrange(sample, desc(KINSHIP))

cat("Rows in candidate_detail:", nrow(candidate_detail), "\n\n")

cat("\n=== DETAILED IBS0 CHECK FOR ALL MULTI-PARENT CANDIDATES ===\n")
print(as.data.frame(candidate_detail))

# Flag which specific relationships look like TRUE parent-offspring (low IBS0)
# vs sibling-type (elevated IBS0), using clone-cluster floor as reference
clone_ibs0_floor <- combined_data %>%
  filter(Relationship == "Highly related / clones") %>%
  summarise(mean_floor = mean(IBS0), q90_floor = quantile(IBS0, 0.90))

cat("\nClone-level IBS0 floor (mean):", round(clone_ibs0_floor$mean_floor, 5), "\n")
cat("Clone-level IBS0 floor (90th percentile):", round(clone_ibs0_floor$q90_floor, 5), "\n\n")

candidate_detail <- candidate_detail %>%
  mutate(
    likely_direct_parent = IBS0 <= clone_ibs0_floor$q90_floor
  )

cat("Candidate relationships consistent with TRUE parent-offspring (low IBS0):\n")
print(as.data.frame(candidate_detail %>% filter(likely_direct_parent)))

cat("\nCandidate relationships consistent with SIBLING-type (elevated IBS0):\n")
print(as.data.frame(candidate_detail %>% filter(!likely_direct_parent)))

write.csv(candidate_detail, "multi_parent_candidates_ibs0_detail.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. full_clone_cluster_lookup.csv\n")
cat("2. comprehensive_multi_parent_check.csv\n")
cat("3. multi_parent_candidates_ibs0_detail.csv\n")

# ==============================================================================
# Break down the 644 multi-cluster candidates for interpretability
# ==============================================================================

# For each candidate, get their SINGLE strongest (max) relationship type
# so we can see if this is "new" cases or just the blue cases we already found
candidate_max_relationship <- symmetric_relationships %>%
  filter(sample %in% multi_parent_candidates$sample) %>%
  group_by(sample) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample, max_partner = partner, max_KINSHIP = KINSHIP, 
         max_IBS0 = IBS0, max_Relationship = Relationship)

# Join back to the candidate summary
multi_parent_candidates_labeled <- multi_parent_candidates %>%
  left_join(candidate_max_relationship, by = "sample")

# Breakdown by max relationship type (clone vs first degree)
cat("=== BREAKDOWN OF 644 CANDIDATES BY THEIR MAX RELATIONSHIP TYPE ===\n")
print(as.data.frame(multi_parent_candidates_labeled %>% count(max_Relationship)))

# Breakdown: is the sample itself a farm sample or reference?
multi_parent_candidates_labeled <- multi_parent_candidates_labeled %>%
  mutate(sample_type = ifelse(sample %in% reference_varieties, "Reference", "Farm"))

cat("\n=== BREAKDOWN BY SAMPLE TYPE (Farm vs Reference) ===\n")
print(as.data.frame(multi_parent_candidates_labeled %>% count(sample_type, max_Relationship)))

# Cross-tab: sample type x whether reference is involved in their multi-match
cat("\n=== CROSS-TAB: sample type x any_reference_involved ===\n")
print(as.data.frame(multi_parent_candidates_labeled %>% 
                      count(sample_type, any_reference_involved)))

# How many of the ORIGINAL 38 "First degree farm-reference max" samples 
# are actually IN this 644 list?
original_38 <- max_ref_relationships %>% 
  filter(Relationship == "First degree") %>% 
  pull(sample)

cat("\nOf the original 38 first-degree farm-reference samples, how many appear in the 644:\n")
cat(sum(original_38 %in% multi_parent_candidates$sample), "out of", length(original_38), "\n")

# IBS0 boxplot
ibs0_comparison_data <- combined_data %>%
  filter(IID1 != IID2) %>%
  filter(Relationship %in% c("Highly related / clones", "First degree"))

p_ibs0_comparison <- ggplot(ibs0_comparison_data, aes(x = Relationship, y = IBS0, fill = Relationship)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = clone_ibs0_floor$q90_floor, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Highly related / clones" = "red", "First degree" = "#0066CC")) +
  labs(
    title = "IBS0: Clone vs First-Degree Relationships",
    subtitle = "Dashed line = clone-level IBS0 90th percentile (parent-offspring threshold)",
    y = "IBS0 Coefficient",
    x = "Relationship Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_ibs0_comparison)
ggsave("ibs0_clone_vs_firstdegree_comparison.png", p_ibs0_comparison, width = 8, height = 6, dpi = 300)

# T test
t.test(IBS0 ~ Relationship, data = combined_data %>% 
         filter(Relationship %in% c("Highly related / clones", "First degree")),
       var.equal = FALSE)  # This confirms Welch's (default), as opposed to var.equal = TRUE for classic Student's t-test


# ==============================================================================
# ALL-PAIRWISE VERSION: Same analysis as max_ref_relationships, but using
# every farm-reference pair (not just each farm's single best match)
# ==============================================================================

# Build the all-pairwise farm x reference dataset (mirrors max_ref_relationships structure)
all_ref_relationships <- combined_data %>%
  filter(IID1 != IID2) %>%
  filter((IID1 %in% reference_varieties & !IID2 %in% reference_varieties) |
           (IID2 %in% reference_varieties & !IID1 %in% reference_varieties)) %>%
  mutate(
    sample = case_when(!IID1 %in% reference_varieties ~ IID1, 
                       !IID2 %in% reference_varieties ~ IID2),
    ref_partner = case_when(IID1 %in% reference_varieties ~ IID1, 
                            IID2 %in% reference_varieties ~ IID2),
    relationship_category = "Farm-Reference"
  ) %>%
  select(sample, ref_partner, KINSHIP, IBS0, Relationship, NSNP, relationship_category)

# Merge in F for farm sample
all_ref_relationships <- all_ref_relationships %>%
  left_join(het_data %>% select(IID, F), by = c("sample" = "IID"))

cat("Total all-pairwise farm-reference relationships:", nrow(all_ref_relationships), "\n")
cat("(compare to max-relationship version which had", nrow(max_ref_relationships), "rows)\n\n")

# ------------------------------------------------------------------------------
# Plot: Kinship vs F (all pairwise)
# ------------------------------------------------------------------------------

p_all_ref_inbreeding <- ggplot(all_ref_relationships, aes(x = KINSHIP, y = F, 
                                                          color = Relationship, 
                                                          shape = relationship_category)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = max(all_ref_relationships$F, na.rm = TRUE) * 0.95, 
           label = c("3rd", "2nd", "1st", "Clone"), 
           color = "gray30", size = 3, angle = 90, vjust = -0.5) +
  scale_color_manual(values = relationship_colors, name = "Relationship Type") +
  scale_shape_manual(values = c("Farm-Reference" = 17), name = "Pairing Type") +
  labs(
    title = "All Pairwise Reference x Farm Relationships",
    subtitle = paste("All", nrow(all_ref_relationships), "farm-reference pairs (not filtered to max)"),
    x = "KING Kinship Coefficient (Φ)",
    y = "Inbreeding Coefficient (F)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  ) +
  coord_cartesian(xlim = c(0, 0.52))

print(p_all_ref_inbreeding)
ggsave("all_pairwise_reference_relationships_inbreeding_plot.png", p_all_ref_inbreeding, 
       width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# ΔF: merge in F_ref, calculate delta_F
# ------------------------------------------------------------------------------

all_ref_relationships <- all_ref_relationships %>%
  left_join(het_data %>% select(IID, F_ref = F), by = c("ref_partner" = "IID")) %>%
  rename(F_farm = F) %>%
  mutate(delta_F = F_farm - F_ref)

cat("Missing F_ref values:", sum(is.na(all_ref_relationships$F_ref)), "\n")

# First Degree only - Kinship vs ΔF (NO correlation line, NO p-value, per your request)
blue_only_all <- all_ref_relationships %>% filter(Relationship == "First degree")

p_blue_deltaF_all <- ggplot(blue_only_all, aes(x = KINSHIP, y = delta_F)) +
  geom_point(color = "#0066CC", size = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "First Degree Relationships: Kinship vs ΔF (Farm - Reference)",
    subtitle = paste("All pairwise first-degree matches, n =", nrow(blue_only_all)),
    x = "KING Kinship Coefficient (Φ)",
    y = "ΔF (Farm Inbreeding − Reference Inbreeding)"
  ) +
  theme_minimal()

print(p_blue_deltaF_all)
ggsave("all_pairwise_first_degree_deltaF.png", p_blue_deltaF_all, width = 10, height = 7, dpi = 300)

# All relationship types - Kinship vs ΔF
p_all_deltaF_allpairs <- ggplot(all_ref_relationships, aes(x = KINSHIP, y = delta_F, color = Relationship)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(0.044, 0.088, 0.19, 0.36), 
             linetype = "dashed", alpha = 0.6, color = "gray30") +
  annotate("text", x = c(0.044, 0.088, 0.19, 0.36), y = max(all_ref_relationships$delta_F, na.rm = TRUE) * 0.95, 
           label = c("3rd", "2nd", "1st", "Clone"), 
           color = "gray30", size = 3, angle = 90, vjust = -0.5) +
  scale_color_manual(values = relationship_colors, name = "Relationship Type") +
  labs(
    title = "Kinship vs ΔF (Farm − Reference) Across All Relationship Types",
    subtitle = paste("All", nrow(all_ref_relationships), "farm-reference pairs (not filtered to max)"),
    x = "KING Kinship Coefficient (Φ)",
    y = "ΔF (Farm Inbreeding − Reference Inbreeding)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

print(p_all_deltaF_allpairs)
ggsave("all_pairwise_kinship_vs_deltaF_all_relationships.png", p_all_deltaF_allpairs, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# F_farm by Relationship category boxplot (all pairwise)
# ------------------------------------------------------------------------------

all_ref_relationships$Relationship <- factor(
  all_ref_relationships$Relationship,
  levels = c("Fourth degree and unrelated", "Third degree", "Second degree", 
             "First degree", "Highly related / clones")
)

F_by_relationship_allpairs <- all_ref_relationships %>%
  group_by(Relationship) %>%
  summarise(
    n = n(),
    pct_negative_F = round(100 * mean(F_farm < 0, na.rm = TRUE), 1),
    mean_F_farm = round(mean(F_farm, na.rm = TRUE), 3),
    median_F_farm = round(median(F_farm, na.rm = TRUE), 3),
    min_F_farm = round(min(F_farm, na.rm = TRUE), 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(pct_negative_F))

print(F_by_relationship_allpairs)

p_F_by_relationship_allpairs <- ggplot(all_ref_relationships, aes(x = Relationship, y = F_farm, fill = Relationship)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = relationship_colors) +
  labs(
    title = "Farm Sample Inbreeding Coefficient (F) by Relationship to Reference",
    subtitle = "All pairwise relationships (not filtered to max)",
    y = "Inbreeding Coefficient (F) - Farm Sample",
    x = "Relationship Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_F_by_relationship_allpairs)
ggsave("all_pairwise_F_farm_by_relationship_category.png", p_F_by_relationship_allpairs, width = 10, height = 7, dpi = 300)