library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(stringr)

# ============================================================
# Load scored human fibroblast dataset
# ============================================================

fib <- readRDS("../../data/human_fib_ucell.rds")

# ============================================================
# Check UCell score columns
# ============================================================

score_cols <- grep("_UCell$", colnames(fib@meta.data), value = TRUE)
score_cols

# Summary statistics
summary_stats <- lapply(score_cols, function(x){
  summary(fib@meta.data[[x]])
})

names(summary_stats) <- score_cols

# ============================================================
# Average UCell scores by fibroblast subtype
# ============================================================

avg_scores <- fib@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(
    across(all_of(score_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# Save average scores
write.csv(
  avg_scores,
  "../../results/07_snRNA/Human_Fibroblast_UCell_MeanScores.csv",
  row.names = FALSE
)

# ============================================================
# Convert to matrix
# ============================================================

heatmap_mat <- avg_scores %>%
  column_to_rownames("author_cell_type") %>%
  as.matrix()

# ============================================================
# Row order
# ============================================================

celltype_order <- c(
  "ABCA10hi fibroblasts",
  "ADAM12hi fibroblasts",
  "FBLN1hi fibroblasts",
  "NR4A1hi fibroblasts"
)

heatmap_mat <- heatmap_mat[celltype_order, ]

# ============================================================
# Column order
# ============================================================

score_order <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell",
  "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

heatmap_mat <- heatmap_mat[, score_order]

# ============================================================
# Pretty labels
# ============================================================

colnames(heatmap_mat) <- c(
  "Integrin–FAK",
  "RhoA–ROCK",
  "YAP/TEAD",
  "TGFβ Signaling",
  "SMAD Signaling",
  "Fibroblast Activation",
  "Myofibroblast",
  "ECM Organization",
  "Inflammation"
)

# ============================================================
# Figure 1
# Mechanotransduction
# ============================================================

heatmap_mech <- heatmap_mat[, c(
  "Integrin–FAK",
  "RhoA–ROCK",
  "YAP/TEAD",
  "TGFβ Signaling",
  "SMAD Signaling"
)]

heatmap_mech <- scale(heatmap_mech)

pheatmap(
  heatmap_mech,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 35,
  cellheight = 25,
  border_color = NA,
  main = "Mechanotransduction and TGFβ/SMAD Signaling",
  filename = "../../results/07_snRNA/Figure1_Mechanotransduction_Heatmap_Human.png",
  width = 10,
  height = 5
)

# ============================================================
# Figure 2
# Fibroblast activation
# ============================================================

heatmap_fib <- heatmap_mat[, c(
  "Fibroblast Activation",
  "Myofibroblast"
)]

pheatmap(
  heatmap_fib,
  scale = "column",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 40,
  cellheight = 25,
  border_color = NA,
  main = "Fibroblast Activation Scores",
  filename = "../../results/07_snRNA/Figure2_FibroblastActivation_Heatmap_Human.png",
  width = 7,
  height = 5
)

# ============================================================
# Figure 3
# ECM organization & inflammation
# ============================================================

heatmap_ecm <- heatmap_mat[, c(
  "ECM Organization",
  "Inflammation"
)]

pheatmap(
  heatmap_ecm,
  scale = "column",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 40,
  cellheight = 25,
  border_color = NA,
  main = "ECM Organization and Inflammation",
  filename = "../../results/07_snRNA/Figure3_ECM_Inflammation_Heatmap_Human.png",
  width = 7,
  height = 5
)



