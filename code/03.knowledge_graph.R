# Test Cytoscape figure and knowledge graph for cassava fingerprinting data
# For max relationships between farm x reference and reference x reference

setwd("~/Documents/cassava_fingerprinting/")
library(dplyr)
library(readr)

# Read both files and combine into one dataset
max_ref_temp <- read.csv("max_ref_relationships.csv")
ref_ref_temp <- read.csv("ref_ref_relationships.csv")

# Combine 
max_ref_relationships <- bind_rows(
  max_ref_temp,
  ref_ref_temp %>% 
    filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree"))
)

# Preview
head(max_ref_relationships)
cat("Total relationships:", nrow(max_ref_relationships), "\n")

# Create edge list for Cytoscape
edges <- max_ref_relationships %>%
  select(
    source = sample,           # Farm sample (source node)
    target = ref_partner,      # Reference variety (target node)  
    KINSHIP,                   # Edge weight/attribute
    IBS0,                      # Edge attribute
    Relationship,              # Edge attribute  
    NSNP,                      # Edge attribute
    relationship_category      # Edge attribute
  ) %>%
  # Add the interaction type as a new column
  mutate(interaction = "connects_to") %>%
  # Remove any rows with missing data
  filter(!is.na(source) & !is.na(target))

# Keep only strong relationships
strong_edges <- edges %>%
  filter(KINSHIP > 0.36 | Relationship %in% c("First degree", "Second degree"))

# Save edge list
write.csv(strong_edges, "cytoscape_edges.csv", row.names = FALSE)

# Get all unique nodes (samples + references)
all_nodes <- unique(c(max_ref_relationships$sample, max_ref_relationships$ref_partner))

# Define reference varieties (same as your original script)
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

# Create node attributes
nodes <- data.frame(
  name = all_nodes,  # Changed from 'id' to 'name'
  node_type = ifelse(all_nodes %in% reference_varieties, "Reference", "Farm"),
  stringsAsFactors = FALSE
) %>%
  # Add additional node attributes
  left_join(
    max_ref_relationships %>%
      group_by(sample) %>%
      summarise(
        max_kinship = max(KINSHIP, na.rm = TRUE),
        strongest_relationship = Relationship[which.max(KINSHIP)],
        .groups = 'drop'
      ),
    by = c("name" = "sample")  # Changed from 'id' to 'name'
  ) %>%
  # For reference nodes, get their info too
  left_join(
    max_ref_relationships %>%
      group_by(ref_partner) %>%
      summarise(
        connections_count = n(),
        .groups = 'drop'
      ),
    by = c("name" = "ref_partner")  # Changed from 'id' to 'name'
  ) %>%
  # Clean up missing values
  mutate(
    max_kinship = ifelse(is.na(max_kinship), 0, max_kinship),
    connections_count = ifelse(is.na(connections_count), 0, connections_count),
    strongest_relationship = ifelse(is.na(strongest_relationship), "None", strongest_relationship)
  )

# Save node list
write.csv(nodes, "cytoscape_nodes.csv", row.names = FALSE)
cat("Node file created:", nrow(nodes), "nodes\n")
cat("Reference nodes:", sum(nodes$node_type == "Reference"), "\n")
cat("Farm nodes:", sum(nodes$node_type == "Farm"), "\n")

# Network summary
cat("\n=== NETWORK SUMMARY FOR CYTOSCAPE ===\n")
cat("Total nodes:", nrow(nodes), "\n")
cat("Total edges:", nrow(edges), "\n")
cat("Reference varieties:", sum(nodes$node_type == "Reference"), "\n")
cat("Farm samples:", sum(nodes$node_type == "Farm"), "\n")

# Edge type distribution
cat("\nRelationship types:\n")
print(table(edges$Relationship))

# Preview files
cat("\n=== EDGE FILE PREVIEW ===\n")
print(head(edges))

cat("\n=== NODE FILE PREVIEW ===\n") 
print(head(nodes))

##### Make something in R instead #####

# ==============================================================================
# Knowledge graph: Farm x Reference and Reference x Reference
# Color = Relationship Type (consistent with other figures)
# Shape = Node Type (Reference vs Farm)
# Includes isolated/weak nodes (3rd/4th degree) as unconnected dots
# ==============================================================================

library(dplyr)
library(visNetwork)

# Combine max farm-reference relationships + strong ref-ref relationships
max_ref_relationships <- bind_rows(
  max_ref_temp,
  ref_ref_temp %>% 
    filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree"))
)

# Build edge list (full, unfiltered)
edges <- max_ref_relationships %>%
  select(source = sample, target = ref_partner, KINSHIP, IBS0, Relationship, NSNP) %>%
  filter(!is.na(source) & !is.na(target))

# Keep only strong relationships for EDGE visualization
vis_edges_ref <- edges %>%
  filter(KINSHIP > 0.088) %>%
  mutate(
    from_name = source,
    to_name = target,
    width = KINSHIP * 10,
    color = case_when(
      Relationship == "Highly related / clones" ~ "red",
      Relationship == "First degree" ~ "#0066CC",
      Relationship == "Second degree" ~ "#FF69B4",
      TRUE ~ "#CCCCCC"
    ),
    title = paste0(
      source, " ↔ ", target, "<br>",
      "Relationship: ", Relationship, "<br>",
      "Kinship: ", round(KINSHIP, 3), "<br>",
      "IBS0: ", round(IBS0, 4), "<br>",
      "NSNP: ", NSNP
    )
  )

# Build node list from FULL dataset (before edge filtering) - keeps isolated nodes
all_nodes_ref <- unique(c(max_ref_relationships$sample, max_ref_relationships$ref_partner))

# Node summary from FULL dataset, not just filtered edges
node_summary_ref <- bind_rows(
  max_ref_relationships %>% select(name = sample, KINSHIP, Relationship),
  max_ref_relationships %>% select(name = ref_partner, KINSHIP, Relationship)
) %>%
  group_by(name) %>%
  summarise(
    max_kinship = max(KINSHIP, na.rm = TRUE),
    strongest_relationship = Relationship[which.max(KINSHIP)],
    n_connections = n(),
    .groups = 'drop'
  )

nodes_ref <- data.frame(
  name = all_nodes_ref,
  node_type = ifelse(all_nodes_ref %in% reference_varieties, "Reference", "Farm"),
  stringsAsFactors = FALSE
) %>%
  left_join(node_summary_ref, by = "name") %>%
  mutate(
    max_kinship = ifelse(is.na(max_kinship), 0, max_kinship),
    strongest_relationship = ifelse(is.na(strongest_relationship), "None", strongest_relationship),
    n_connections = ifelse(is.na(n_connections), 0, n_connections)
  )

# Prepare nodes for visNetwork - shape by node type, color by relationship type
vis_nodes_ref <- nodes_ref %>%
  mutate(
    id = 1:nrow(nodes_ref),
    label = name,
    group = strongest_relationship,
    title = paste0(
      "<b>", name, "</b><br>",
      "Type: ", node_type, "<br>",
      "Max Kinship: ", round(max_kinship, 3), "<br>",
      "Strongest Relationship: ", strongest_relationship, "<br>",
      "Connections: ", n_connections
    ),
    # Shape by node type
    shape = ifelse(node_type == "Reference", "triangle", "dot"),
    # Size - slightly bigger for references so they stand out
    size = ifelse(node_type == "Reference", 25, 15),
    # Color by relationship type (matching your other figures)
    color = case_when(
      strongest_relationship == "Highly related / clones" ~ "red",
      strongest_relationship == "First degree" ~ "#0066CC",
      strongest_relationship == "Second degree" ~ "#FF69B4",
      TRUE ~ "#CCCCCC"
    )
  ) %>%
  select(id, label, group, title, shape, size, color)

# Prepare final edges with numeric IDs (using the SAME node list, so IDs match)
vis_edges_ref_final <- vis_edges_ref %>%
  mutate(
    from = match(from_name, nodes_ref$name),
    to = match(to_name, nodes_ref$name)
  ) %>%
  filter(!is.na(from) & !is.na(to)) %>%
  select(from, to, width, color, title)

cat("Total nodes:", nrow(vis_nodes_ref), "\n")
cat("Total edges:", nrow(vis_edges_ref_final), "\n")
cat("Isolated nodes (no edge drawn):", nrow(vis_nodes_ref) - length(unique(c(vis_edges_ref_final$from, vis_edges_ref_final$to))), "\n")

# Create the interactive knowledge graph
ref_knowledge_graph <- visNetwork(vis_nodes_ref, vis_edges_ref_final,
                                  main = "Cassava Farm x Reference Genetic Relationships",
                                  submain = "Shape = sample type | Color = relationship strength") %>%
  visLayout(randomSeed = 123) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50),
    stabilization = FALSE
  ) %>%
  visInteraction(
    navigationButtons = TRUE,
    hover = TRUE,
    selectConnectedEdges = TRUE,
    tooltipDelay = 100
  ) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    selectedBy = "group",
    nodesIdSelection = TRUE
  ) %>%
  visEdges(
    arrows = "none",
    smooth = list(enabled = TRUE, type = "continuous"),
    font = list(size = 10),
    scaling = list(min = 1, max = 10)
  ) %>%
  visNodes(
    borderWidth = 2,
    font = list(size = 16, face = "arial"),
    scaling = list(min = 10, max = 30)
  ) %>%
  visLegend(
    addNodes = list(
      list(label = "Reference Variety", shape = "triangle", color = "grey", size = 20),
      list(label = "Farm Sample", shape = "dot", color = "grey", size = 15)
    ),
    addEdges = data.frame(
      label = c("Highly related / clones", "First degree", "Second degree"),
      color = c("red", "#0066CC", "#FF69B4")
    ),
    useGroups = FALSE,
    position = "right"
  )

# Display
ref_knowledge_graph

# Save as HTML file
visSave(ref_knowledge_graph, "cassava_knowledge_graph_v2.html")

cat("\n=== KNOWLEDGE GRAPH SAVED ===\n")
cat("Nodes:", nrow(vis_nodes_ref), "\n")
cat("Edges:", nrow(vis_edges_ref_final), "\n")
cat("Reference nodes:", sum(nodes_ref$node_type == "Reference"), "\n")
cat("Farm nodes:", sum(nodes_ref$node_type == "Farm"), "\n")

# ==============================================================================
# Knowledge graph for farm samples lacking a strong reference match
# Using MAX relationship per sample (not all pairwise relationships)
# ==============================================================================

library(dplyr)
library(visNetwork)

# Step 1: Get ALL pairwise relationships among unmatched farms from combined_data
unmatched_all_relationships <- combined_data %>%
  filter(IID1 %in% farms_without_ref_match & IID2 %in% farms_without_ref_match) %>%
  filter(IID1 != IID2) %>%
  select(source = IID1, target = IID2, KINSHIP, IBS0, Relationship, NSNP)

# Step 2: For each sample, keep only their MAX relationship (strongest connection)
unmatched_max_relationships <- unmatched_all_relationships %>%
  group_by(source) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Total max relationships (edges):", nrow(unmatched_max_relationships), "\n")
print(table(unmatched_max_relationships$Relationship))

# Step 3: Build node list
all_nodes_unmatched <- unique(c(unmatched_max_relationships$source, unmatched_max_relationships$target))

node_summary_unmatched <- bind_rows(
  unmatched_max_relationships %>% select(name = source, KINSHIP, Relationship),
  unmatched_max_relationships %>% select(name = target, KINSHIP, Relationship)
) %>%
  group_by(name) %>%
  summarise(
    max_kinship = max(KINSHIP, na.rm = TRUE),
    strongest_relationship = Relationship[which.max(KINSHIP)],
    .groups = 'drop'
  )

nodes_unmatched <- data.frame(
  name = all_nodes_unmatched,
  node_type = "Farm",
  stringsAsFactors = FALSE
) %>%
  left_join(node_summary_unmatched, by = "name") %>%
  mutate(
    max_kinship = ifelse(is.na(max_kinship), 0, max_kinship),
    strongest_relationship = ifelse(is.na(strongest_relationship), "None", strongest_relationship)
  )


# Step 4: Prepare nodes for visNetwork
vis_nodes_unmatched <- nodes_unmatched %>%
  mutate(
    id = 1:nrow(nodes_unmatched),
    label = name,
    group = strongest_relationship,
    title = paste0(
      "<b>", name, "</b><br>",
      "Max Kinship: ", round(max_kinship, 3), "<br>",
      "Strongest Relationship: ", strongest_relationship
    ),
    size = 15,
    color = case_when(
      strongest_relationship == "Highly related / clones" ~ "red",
      strongest_relationship == "First degree" ~ "#0066CC",
      strongest_relationship == "Second degree" ~ "#FF69B4",
      strongest_relationship == "Third degree" ~ "#8A2BE2",
      TRUE ~ "#CCCCCC"
    )
  ) %>%
  select(id, label, group, title, size, color)

# Step 5: Prepare edges for visNetwork
vis_edges_unmatched <- unmatched_max_relationships %>%
  mutate(
    from = match(source, nodes_unmatched$name),
    to = match(target, nodes_unmatched$name),
    width = KINSHIP * 10,
    color = case_when(
      Relationship == "Highly related / clones" ~ "red",
      Relationship == "First degree" ~ "#0066CC",
      Relationship == "Second degree" ~ "#FF69B4",
      Relationship == "Third degree" ~ "#8A2BE2",
      TRUE ~ "#CCCCCC"
    ),
    title = paste0(
      source, " ↔ ", target, "<br>",
      "Relationship: ", Relationship, "<br>",
      "Kinship: ", round(KINSHIP, 3), "<br>",
      "IBS0: ", round(IBS0, 4), "<br>",
      "NSNP: ", NSNP
    )
  ) %>%
  filter(!is.na(from) & !is.na(to)) %>%
  select(from, to, width, color, title)

cat("\nTotal edges for visualization:", nrow(vis_edges_unmatched), "\n")

# Step 6: Create the interactive knowledge graph
unmatched_knowledge_graph <- visNetwork(vis_nodes_unmatched, vis_edges_unmatched) %>%
  visLayout(randomSeed = 123) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50),
    stabilization = FALSE
  ) %>%
  visInteraction(
    navigationButtons = TRUE,
    hover = TRUE,
    selectConnectedEdges = TRUE,
    tooltipDelay = 100
  ) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    selectedBy = "group",
    nodesIdSelection = TRUE
  ) %>%
  visEdges(
    arrows = "none",
    smooth = list(enabled = TRUE, type = "continuous"),
    font = list(size = 10),
    scaling = list(min = 1, max = 10)
  ) %>%
  visNodes(
    borderWidth = 2,
    font = list(size = 16, face = "arial"),
    scaling = list(min = 10, max = 30)
  ) %>%
  visLegend(
    addNodes = list(
      list(label = "Highly related / clones", color = "red", size = 15),
      list(label = "First degree", color = "#0066CC", size = 15),
      list(label = "Second degree", color = "#FF69B4", size = 15),
      list(label = "Third degree", color = "#8A2BE2", size = 15),
      list(label = "Fourth degree / unrelated", color = "#CCCCCC", size = 15)
    ),
    useGroups = FALSE,
    position = "right"
  )

# Display
unmatched_knowledge_graph

# Save as HTML file
visSave(unmatched_knowledge_graph, "unmatched_farms_knowledge_graph.html")

cat("\n=== KNOWLEDGE GRAPH SAVED ===\n")
cat("Nodes:", nrow(vis_nodes_unmatched), "\n")
cat("Edges:", nrow(vis_edges_unmatched), "\n")


# ==============================================================================
# Knowledge graph for farm samples lacking a strong reference match
# ONLY "Highly related / clones" (red) relationships
# Using MAX relationship per sample
# ==============================================================================

library(dplyr)
library(visNetwork)

# Step 1: Get ALL pairwise relationships among unmatched farms from combined_data
unmatched_all_relationships <- combined_data %>%
  filter(IID1 %in% farms_without_ref_match & IID2 %in% farms_without_ref_match) %>%
  filter(IID1 != IID2) %>%
  select(source = IID1, target = IID2, KINSHIP, IBS0, Relationship, NSNP)

# Step 2: For each sample, keep only their MAX relationship
unmatched_max_relationships <- unmatched_all_relationships %>%
  group_by(source) %>%
  slice_max(KINSHIP, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 3: Filter to ONLY "Highly related / clones"
red_only_relationships <- unmatched_max_relationships %>%
  filter(Relationship == "Highly related / clones")

cat("Total 'Highly related/clones' edges:", nrow(red_only_relationships), "\n")

# Step 4: Build node list
all_nodes_red <- unique(c(red_only_relationships$source, red_only_relationships$target))

node_summary_red <- bind_rows(
  red_only_relationships %>% select(name = source, KINSHIP, Relationship),
  red_only_relationships %>% select(name = target, KINSHIP, Relationship)
) %>%
  group_by(name) %>%
  summarise(
    max_kinship = max(KINSHIP, na.rm = TRUE),
    strongest_relationship = Relationship[which.max(KINSHIP)],
    .groups = 'drop'
  )

nodes_red <- data.frame(
  name = all_nodes_red,
  node_type = "Farm",
  stringsAsFactors = FALSE
) %>%
  left_join(node_summary_red, by = "name") %>%
  mutate(
    max_kinship = ifelse(is.na(max_kinship), 0, max_kinship),
    strongest_relationship = ifelse(is.na(strongest_relationship), "Highly related / clones", strongest_relationship)
  )

# Step 5: Prepare nodes for visNetwork
vis_nodes_red <- nodes_red %>%
  mutate(
    id = 1:nrow(nodes_red),
    label = name,
    group = strongest_relationship,
    title = paste0(
      "<b>", name, "</b><br>",
      "Max Kinship: ", round(max_kinship, 3)
    ),
    size = 15,
    color = "red"
  ) %>%
  select(id, label, group, title, size, color)

# Step 6: Prepare edges for visNetwork
vis_edges_red <- red_only_relationships %>%
  mutate(
    from = match(source, nodes_red$name),
    to = match(target, nodes_red$name),
    width = KINSHIP * 10,
    color = "red",
    title = paste0(
      source, " ↔ ", target, "<br>",
      "Relationship: ", Relationship, "<br>",
      "Kinship: ", round(KINSHIP, 3), "<br>",
      "IBS0: ", round(IBS0, 4), "<br>",
      "NSNP: ", NSNP
    )
  ) %>%
  filter(!is.na(from) & !is.na(to)) %>%
  select(from, to, width, color, title)

cat("Total edges for visualization:", nrow(vis_edges_red), "\n")

# Step 7: Create the interactive knowledge graph
red_knowledge_graph <- visNetwork(vis_nodes_red, vis_edges_red,
                                  main = "Farm x Farm Clonal Relationships",
                                  submain = "Farm samples lacking a reference match, connected to most related farm sample") %>%
  visLayout(randomSeed = 123) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50),
    stabilization = FALSE
  ) %>%
  visInteraction(
    navigationButtons = TRUE,
    hover = TRUE,
    selectConnectedEdges = TRUE,
    tooltipDelay = 100
  ) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    nodesIdSelection = TRUE
  ) %>%
  visEdges(
    arrows = "none",
    smooth = list(enabled = TRUE, type = "continuous"),
    font = list(size = 10),
    scaling = list(min = 1, max = 10)
  ) %>%
  visNodes(
    borderWidth = 2,
    font = list(size = 16, face = "arial"),
    scaling = list(min = 10, max = 30)
  ) %>%
  visLegend(
    addNodes = list(
      list(label = "Highly related / clones", color = "red", size = 15)
    ),
    useGroups = FALSE,
    position = "right"
  )

# Display
red_knowledge_graph

# Save as HTML file
visSave(red_knowledge_graph, "red_only_farms_knowledge_graph.html")

cat("\n=== KNOWLEDGE GRAPH SAVED ===\n")
cat("Nodes:", nrow(vis_nodes_red), "\n")
cat("Edges:", nrow(vis_edges_red), "\n")



# ==============================================================================
# Knowledge graph for Reference x Reference relationships
# "Highly related / clones" (red), "First degree" (blue), 
# "Second degree" (pink), "Third degree" (purple)
# ALL pairwise relationships (not just max)
# ==============================================================================

library(dplyr)
library(visNetwork)

# Step 1: Get ALL pairwise reference-reference relationships, filtered to red + blue + pink + purple
ref_ref_all_relationships <- combined_data %>%
  filter(IID1 %in% reference_varieties & IID2 %in% reference_varieties) %>%
  filter(IID1 != IID2) %>%
  filter(Relationship %in% c("Highly related / clones", "First degree", "Second degree", "Third degree")) %>%
  select(source = IID1, target = IID2, KINSHIP, IBS0, Relationship, NSNP)

cat("Total reference-reference relationships (clone + 1st + 2nd + 3rd degree):", nrow(ref_ref_all_relationships), "\n")
print(table(ref_ref_all_relationships$Relationship))

all_nodes_refref <- unique(c(ref_ref_all_relationships$source, ref_ref_all_relationships$target))

node_summary <- bind_rows(
  ref_ref_all_relationships %>% select(name = source, KINSHIP, Relationship),
  ref_ref_all_relationships %>% select(name = target, KINSHIP, Relationship)
) %>%
  group_by(name) %>%
  summarise(
    max_kinship = max(KINSHIP, na.rm = TRUE),
    strongest_relationship = Relationship[which.max(KINSHIP)],
    n_connections = n(),
    .groups = 'drop'
  )

nodes_refref <- data.frame(
  name = all_nodes_refref,
  node_type = "Reference",
  stringsAsFactors = FALSE
) %>%
  left_join(node_summary, by = "name") %>%
  mutate(
    max_kinship = ifelse(is.na(max_kinship), 0, max_kinship),
    strongest_relationship = ifelse(is.na(strongest_relationship), "None", strongest_relationship),
    n_connections = ifelse(is.na(n_connections), 0, n_connections)
  )

# Step 3: Prepare nodes for visNetwork
vis_nodes_refref <- nodes_refref %>%
  mutate(
    id = 1:nrow(nodes_refref),
    label = name,
    group = strongest_relationship,
    title = paste0(
      "<b>", name, "</b><br>",
      "Max Kinship: ", round(max_kinship, 3), "<br>",
      "Strongest Relationship: ", strongest_relationship, "<br>",
      "Connections: ", n_connections
    ),
    size = 20,
    color = case_when(
      strongest_relationship == "Highly related / clones" ~ "red",
      strongest_relationship == "First degree" ~ "#0066CC",
      strongest_relationship == "Second degree" ~ "#FF69B4",
      strongest_relationship == "Third degree" ~ "#8A2BE2",
      TRUE ~ "#CCCCCC"
    )
  ) %>%
  select(id, label, group, title, size, color)

# Step 4: Prepare edges for visNetwork
vis_edges_refref <- ref_ref_all_relationships %>%
  mutate(
    from = match(source, nodes_refref$name),
    to = match(target, nodes_refref$name),
    width = KINSHIP * 10,
    color = case_when(
      Relationship == "Highly related / clones" ~ "red",
      Relationship == "First degree" ~ "#0066CC",
      Relationship == "Second degree" ~ "#FF69B4",
      Relationship == "Third degree" ~ "#8A2BE2",
      TRUE ~ "#CCCCCC"
    ),
    title = paste0(
      source, " ↔ ", target, "<br>",
      "Relationship: ", Relationship, "<br>",
      "Kinship: ", round(KINSHIP, 3), "<br>",
      "IBS0: ", round(IBS0, 4), "<br>",
      "NSNP: ", NSNP
    )
  ) %>%
  filter(!is.na(from) & !is.na(to)) %>%
  select(from, to, width, color, title)

cat("\nTotal edges for visualization:", nrow(vis_edges_refref), "\n")

# Step 5: Create the interactive knowledge graph
refref_knowledge_graph <- visNetwork(vis_nodes_refref, vis_edges_refref,
                                     main = "Reference x Reference Relationships",
                                     submain = "Clonal, First, Second, and Third Degree Connections") %>%
  visLayout(randomSeed = 123) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50),
    stabilization = FALSE
  ) %>%
  visInteraction(
    navigationButtons = TRUE,
    hover = TRUE,
    selectConnectedEdges = TRUE,
    tooltipDelay = 100
  ) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    selectedBy = "group",
    nodesIdSelection = TRUE
  ) %>%
  visEdges(
    arrows = "none",
    smooth = list(enabled = TRUE, type = "continuous"),
    font = list(size = 10),
    scaling = list(min = 1, max = 10)
  ) %>%
  visNodes(
    borderWidth = 2,
    font = list(size = 16, face = "arial"),
    scaling = list(min = 10, max = 30)
  ) %>%
  visLegend(
    addNodes = list(
      list(label = "Highly related / clones", color = "red", size = 15),
      list(label = "First degree", color = "#0066CC", size = 15),
      list(label = "Second degree", color = "#FF69B4", size = 15),
      list(label = "Third degree", color = "#8A2BE2", size = 15)
    ),
    useGroups = FALSE,
    position = "right"
  )

# Display
refref_knowledge_graph

# Save as HTML file
visSave(refref_knowledge_graph, "ref_ref_clone_1st_2nd_3rd_degree_graph.html")

cat("\n=== KNOWLEDGE GRAPH SAVED ===\n")
cat("Nodes:", nrow(vis_nodes_refref), "\n")
cat("Edges:", nrow(vis_edges_refref), "\n")

