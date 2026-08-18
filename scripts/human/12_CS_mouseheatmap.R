############################################################
## Mouse UCell Heatmap
## I7D mechanotransduction + TGFβ/SMAD
############################################################

library(Seurat)
library(dplyr)
library(pheatmap)
library(tibble)

############################################################
# Load mouse data
############################################################

mesenchymal <- readRDS(
  "../../data/mesenchymal_geneset_scored_refined.rds"
)

############################################################
# Pathways
############################################################

pathways <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell"
)

############################################################
# Pathway labels
############################################################

pathway_labels <- c(
  "Integrin–FAK",
  "RhoA–ROCK",
  "YAP/TEAD",
  "TGFβ",
  "SMAD"
)

############################################################
# Mouse populations
############################################################

mouse_populations <- c(
  "Repair fibroblasts",
  "Repair-activated stromal cells"
)

############################################################
# Calculate mean UCell scores
############################################################

mouse_scores <- mesenchymal@meta.data %>%
  
  filter(
    condition == "I7D",
    cell_type_refined %in% mouse_populations
  ) %>%
  
  group_by(
    cell_type_refined
  ) %>%
  
  summarise(
    across(
      all_of(pathways),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

############################################################
# Convert to matrix
############################################################

mouse_mat <- mouse_scores %>%
  
  column_to_rownames(
    "cell_type_refined"
  ) %>%
  
  t()

############################################################
# Rename rows
############################################################

rownames(mouse_mat) <- pathway_labels

############################################################
# Rename columns
############################################################

colnames(mouse_mat) <- c(
  "Repair Fibroblasts",
  "Repair-Activated Stromal"
)

############################################################
# Inspect actual UCell scores
############################################################

print(
  round(mouse_mat, 3)
)

############################################################
# Save matrix
############################################################

write.csv(
  mouse_mat,
  "../../results/07_snRNA/Mouse_I7D_Mechanotransduction_TGFb_Matrix.csv"
)

############################################################
# Create heatmap
############################################################

pheatmap(
  mouse_mat,
  
  # IMPORTANT: no scaling
  scale = "none",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  border_color = NA,
  
  fontsize_row = 11,
  fontsize_col = 11,
  
  angle_col = 45,
  
  # Show actual UCell scores
  display_numbers = TRUE,
  
  number_format = "%.2f",
  
  fontsize_number = 10,
  
  main =
    "Mouse I7D Mechanotransduction and TGFβ/SMAD Signatures",
  
  filename =
    "../../results/07_snRNA/Mouse_I7D_Mechanotransduction_TGFb_Heatmap.png",
  
  width = 7,
  height = 6
)
