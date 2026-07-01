library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(scCustomize)
library(stringr)

# Load in data 
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")
gene_sets <- readRDS("../data/all_gene_sets_refined.rds")

# run ucell scoring
mesenchymal <- AddModuleScore_UCell(
  mesenchymal,
  features = gene_sets
)

# check score columns
score_cols <- grep("_UCell$", colnames(mesenchymal@meta.data), value = TRUE)
score_cols

# summary statistics 
summary_stats <- lapply(score_cols, function(x){
  summary(mesenchymal@meta.data[[x]])
})

names(summary_stats) <- score_cols

saveRDS(mesenchymal, "../data/mesenchymal_geneset_scored_refined.rds")

# average scores per gene set by cell type
avg_scores <- mesenchymal@meta.data %>%
  group_by(cell_type_refined) %>%
  summarise(
    across(all_of(score_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# convert to matrix
heatmap_mat <- avg_scores %>%
  column_to_rownames("cell_type_refined") %>%
  as.matrix()

# row order
celltype_order <- c(
  "Mature tenocytes",
  "Stromal progenitor-like cells",
  "Signaling stromal cells",
  "Fibrochondrocyte-like tenocytes",
  "Activated tenocytes",
  "Proinflammatory mesenchymal cells",
  "Proliferating mesenchymal cells",
  "Repair-activated stromal cells",
  "Repair fibroblasts",
  "ECM-remodelling tenocytes"
)

# column order using original UCell names
score_order <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell",
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

# apply row and column order
heatmap_mat <- heatmap_mat[celltype_order, score_order]

# clean labels for plotting
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

heatmap_mech <- heatmap_mat[, c(
  "Integrin–FAK",
  "RhoA–ROCK",
  "YAP/TEAD",
  "TGFβ Signaling",
  "SMAD Signaling"
)]

pheatmap(
  heatmap_mech,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 35,
  cellheight = 25,
  main = "Mechanotransduction and TGFβ/SMAD Signaling",
  border_color = NA,
  filename = "../figures2/Figure1_Mechanotransduction_Heatmap.png",
  width = 10,
  height = 6
)

heatmap_fib <- heatmap_mat[, c(
  "Fibroblast Activation",
  "Myofibroblast"
)]

pheatmap(
  heatmap_fib,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 40,
  cellheight = 25,
  main = "Fibroblast Activation Scores",
  border_color = NA,
  filename = "../figures2/Figure2_FibroblastActivation_Heatmap.png",
  width = 8,
  height = 6
)

heatmap_ecm <- heatmap_mat[, c(
  "ECM Organization",
  "Inflammation"
)]

pheatmap(
  heatmap_ecm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 40,
  cellheight = 25,
  main = "ECM Organization and Inflammation",
  border_color = NA,
  filename = "../figures2/Figure3_ECM_Inflammation_Heatmap.png",
  width = 9,
  height = 6
)
