# Load libraries ----
library(dplyr)
library(ggplot2)

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
# 3. CREATE MAX_REF_RELATIONSHIPS.CSV
# ==============================================================================

# Get all farm samples and their strongest reference relationship ----
farm_max_ref <- combined_data %>%
  filter(IID1 %in% reference_varieties | IID2 %in% reference_varieties) %>%
  mutate(
    farm_sample = case_when(
      IID1 %in% reference_varieties & IID2 %in% reference_varieties ~ NA_character_,
      IID1 %in% reference_varieties ~ IID2,
      IID2 %in% reference_varieties ~ IID1,
      TRUE ~ NA_character_
    ),
    ref_partner = case_when(
      IID1 %in% reference_varieties & IID2 %in% reference_varieties ~ NA_character_,
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

# Get reference-reference relationships ----
ref_max_ref <- combined_data %>%
  filter(IID1 %in% reference_varieties & IID2 %in% reference_varieties) %>%
  group_by(IID1) %>%
  slice_min(abs(KINSHIP - 0.5), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample = IID1, ref_partner = IID2, KINSHIP, IBS0, Relationship, NSNP) %>%
  mutate(relationship_category = "Reference-Reference")

# Combine and save ----
max_ref_relationships <- bind_rows(farm_max_ref, ref_max_ref)
write.csv(max_ref_relationships, "max_ref_relationships.csv", row.names = FALSE)

cat("=== MAX REFERENCE RELATIONSHIPS SAVED ===\n")
cat("Total samples:", nrow(max_ref_relationships), "\n")

# ==============================================================================
# 4. CREATE MAXIMUM REFERENCE RELATIONSHIPS PLOT
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
    title = "Maximum Reference Relationships per Sample",
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