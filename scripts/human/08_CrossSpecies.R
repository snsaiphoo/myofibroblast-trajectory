############################################################
## Cross-species mechanotransduction comparison
## Mouse tendon repair vs Human tendon rupture
############################################################

library(Seurat)
library(dplyr)

############################################################
# Load data
############################################################

mesenchymal <- readRDS("../../data/mesenchymal_geneset_scored_refined.rds")
human <- readRDS("../../data/human_fib_ucell.rds")

############################################################
# UCell pathway names
############################################################

pathways <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

############################################################
# Mouse Repair Fibroblasts (I7D)
############################################################

mouse_rf <- subset(
  mesenchymal,
  subset =
    cell_type_refined == "Repair fibroblasts" &
    condition == "I7D"
)

mouse_rf_mean <- colMeans(
  mouse_rf@meta.data[, pathways],
  na.rm = TRUE
)

############################################################
# Mouse Repair-Activated Stromal Cells (I7D)
############################################################

mouse_ras <- subset(
  mesenchymal,
  subset =
    cell_type_refined == "Repair-activated stromal cells" &
    condition == "I7D"
)

mouse_ras_mean <- colMeans(
  mouse_ras@meta.data[, pathways],
  na.rm = TRUE
)

############################################################
# Human fibroblast populations
############################################################

human_sub <- subset(
  human,
  subset = author_cell_type %in% c(
    "ADAM12hi fibroblasts",
    "ABCA10hi fibroblasts",
    "FBLN1hi fibroblasts",
    "NR4A1hi fibroblasts"
  )
)

human_mean <- human_sub@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(
    across(
      all_of(pathways),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

############################################################
# Build comparison matrix
############################################################

# Start with the two mouse columns
comparison <- data.frame(
  Mouse_Repair_Fibroblasts = mouse_rf_mean,
  Mouse_Repair_Activated_Stromal = mouse_ras_mean
)

# Add one human population at a time
for (cell in human_mean$author_cell_type) {
  
  vals <- human_mean %>%
    filter(author_cell_type == cell) %>%
    select(-author_cell_type) %>%
    unlist(use.names = FALSE)
  
  comparison[[cell]] <- vals
}

# Use pathway names as row names
rownames(comparison) <- pathways

############################################################
# Rename pathways
############################################################

rownames(comparison) <- c(
  "Integrin-FAK",
  "RhoA-ROCK",
  "YAP/TEAD",
  "TGFβ",
  "SMAD",
  "Fibroblast Activation",
  "Myofibroblast",
  "Inflammation"
)

############################################################
# Reorder columns
############################################################

comparison <- comparison[, c(
  "Mouse_Repair_Fibroblasts",
  "Mouse_Repair_Activated_Stromal",
  "ADAM12hi fibroblasts",
  "ABCA10hi fibroblasts",
  "FBLN1hi fibroblasts",
  "NR4A1hi fibroblasts"
)]

############################################################
# Inspect results
############################################################

round(comparison, 3)

############################################################
# Save matrix
############################################################

write.csv(
  comparison,
  "../../results/07_snRNA/CrossSpecies_PathwayMatrix.csv",
  quote = FALSE
)

saveRDS(
  comparison,
  "../../data/CrossSpecies_PathwayMatrix.rds"
)
