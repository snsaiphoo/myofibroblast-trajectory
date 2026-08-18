############################################################
## Human UCell Heatmap
## Ruptured tendon fibroblast populations
## Donor-level averaging
############################################################

library(Seurat)
library(dplyr)
library(pheatmap)
library(tibble)

############################################################
# Load human data
############################################################

human <- readRDS(
  "../../data/human_fib_ucell.rds"
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
# Human fibroblast populations
############################################################

human_populations <- c(
  "ADAM12hi fibroblasts",
  "NR4A1hi fibroblasts",
  "ABCA10hi fibroblasts",
  "FBLN1hi fibroblasts"
)

############################################################
# Step 1:
# Calculate mean UCell score for each
# DONOR × CELL TYPE
############################################################

human_donor_scores <- human@meta.data %>%
  
  filter(
    !is.na(donor_id),
    condition == "Rupture",
    author_cell_type %in% human_populations
  ) %>%
  
  group_by(
    donor_id,
    author_cell_type
  ) %>%
  
  summarise(
    across(
      all_of(pathways),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  )

############################################################
# Inspect donor-level results
############################################################

print(human_donor_scores)

############################################################
# Step 2:
# Average donor-level scores
#
# Each donor contributes equally
############################################################

human_scores <- human_donor_scores %>%
  
  group_by(
    author_cell_type
  ) %>%
  
  summarise(
    across(
      all_of(pathways),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_donors = n(),
    .groups = "drop"
  )

############################################################
# Inspect final scores
############################################################

print(
  human_scores %>%
    mutate(
      across(all_of(pathways), ~ round(.x, 3))
    )
)

############################################################
# Convert to matrix
############################################################

human_mat <- human_scores %>%
  
  select(
    author_cell_type,
    all_of(pathways)
  ) %>%
  
  column_to_rownames(
    "author_cell_type"
  ) %>%
  
  t()

############################################################
# Rename rows
############################################################

rownames(human_mat) <- pathway_labels

############################################################
# Rename columns
############################################################

colnames(human_mat) <- c(
  "ADAM12hi",
  "NR4A1hi",
  "ABCA10hi",
  "FBLN1hi"
)

############################################################
# Inspect actual UCell scores
############################################################

print(
  round(human_mat, 3)
)

############################################################
# Save donor-level matrix
############################################################

write.csv(
  human_mat,
  "../../results/07_snRNA/Human_Rupture_DonorLevel_Mechanotransduction_TGFb_Matrix.csv"
)

############################################################
# Create heatmap
############################################################

pheatmap(
  human_mat,
  
  # IMPORTANT:
  # No scaling
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
    "Human Tendon Rupture Mechanotransduction and TGFβ/SMAD",
  
  filename =
    "../../results/07_snRNA/Human_Rupture_DonorLevel_Mechanotransduction_TGFb_Heatmap.png",
  
  width = 8,
  height = 6
)
