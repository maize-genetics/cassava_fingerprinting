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
    filter(Relationship %in% c("Vegetative clones", "First degree", "Second degree"))
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

library(visNetwork)
library(dplyr)

#  Prepare nodes for visNetwork (needs numeric IDs)
vis_nodes <- nodes %>%
  mutate(
    id = 1:nrow(nodes),  # visNetwork needs numeric IDs
    label = name,
    group = node_type,
    # Create hover tooltips - different for Reference vs Farm
    title = ifelse(
      node_type == "Reference",
      # Reference tooltip - just basic info
      paste0(
        "<b>", name, "</b><br>",
        "Type: ", node_type, "<br>",
        "Connections: ", connections_count
      ),
      # Farm sample tooltip - detailed kinship info
      paste0(
        "<b>", name, "</b><br>",
        "Type: ", node_type, "<br>",
        "Max Kinship: ", round(max_kinship, 3), "<br>",
        "Strongest Relationship: ", strongest_relationship, "<br>"
      )
    ),
    # Node sizes based on type
    size = ifelse(node_type == "Reference", 25, 15),
    # Node colors
    color = ifelse(node_type == "Reference", "#8B0000", "#2E8B57")
  ) %>%
  select(id, label, group, title, size, color)

# Prepare edges for visNetwork
vis_edges <- edges %>%
  # Filter for stronger relationships to avoid overwhelm
  filter(KINSHIP > 0.088) %>%
  mutate(
    from = match(source, nodes$name),  # Convert to numeric IDs
    to = match(target, nodes$name),
    # Edge width based on kinship strength
    width = KINSHIP * 10,
    # Edge colors based on relationship type
    color = case_when(
      Relationship == "First degree" ~ "#ff4444",
      Relationship == "Second degree" ~ "#ff8800", 
      Relationship == "Vegetative clones" ~ "#44ff44",
      TRUE ~ "#888888"
    ),
    # Hover tooltip for edges
    title = paste0(
      source, " ↔ ", target, "<br>",
      "Relationship: ", Relationship, "<br>",
      "Kinship: ", round(KINSHIP, 3), "<br>",
      "IBS0: ", round(IBS0, 4), "<br>",
      "NSNP: ", NSNP
    )
  ) %>%
  filter(!is.na(from) & !is.na(to)) %>%  # Remove any NA matches
  select(from, to, width, color, title)

# Create the interactive knowledge graph
knowledge_graph <- visNetwork(vis_nodes, vis_edges) %>%
  # Layout options
  visLayout(randomSeed = 123) %>%  # For reproducible layout
  
  # Physics settings for better clustering
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50),
    stabilization = FALSE
  ) %>%
  
  # Interactive options
  visInteraction(
    navigationButtons = TRUE,
    hover = TRUE,
    selectConnectedEdges = TRUE,
    tooltipDelay = 100
  ) %>%
  
  # Highlighting options
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    selectedBy = "group",
    nodesIdSelection = TRUE
    
  ) %>%
  
  # Edge styling
  visEdges(
    arrows = "none",
    smooth = list(enabled = TRUE, type = "continuous"),
    font = list(size = 10),
    scaling = list(min = 1, max = 10)
  ) %>%
  
  # Node styling  
  visNodes(
    borderWidth = 2,
    font = list(size = 16, face = "arial"),
    scaling = list(min = 10, max = 30)
  ) %>%
  
  # Legend
  visLegend(
    addNodes = list(
      list(label = "Reference Variety", color = "#8B0000", size = 25),
      list(label = "Farm Sample", color = "#2E8B57", size = 15)
    ),
    useGroups = FALSE,
    position = "right"
  )

# Display the knowledge graph
knowledge_graph

# Save as HTML file
visSave(knowledge_graph, "cassava_knowledge_graph.html")

# Print summary
cat("Knowledge Graph Created!\n")
cat("Nodes:", nrow(vis_nodes), "\n")
cat("Edges:", nrow(vis_edges), "\n") 
cat("Reference varieties:", sum(vis_nodes$group == "Reference"), "\n")
cat("Farm samples:", sum(vis_nodes$group == "Farm"), "\n")


