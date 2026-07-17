# Cross species comparisons 
library(Seurat)
library(dplyr)
library(tibble)
library(pheatmap)

mesenchymal <- readRDS("../../data/mesenchymal_geneset_scored_refined.rds")
m_gene_sets <- readRDS("../../data/all_gene_sets_refined.rds")
human <- readRDS("../../data/human_fib_ucell.rds")
h_gene_sets <- readRDS("../../data/all_gene_sets_human_ensembl.rds")

# Mouse pathways (shared between species)
score_order <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

# Mouse
mouse_avg <-
  mesenchymal@meta.data %>%
  group_by(cell_type_refined) %>%  
  summarise(across(all_of(score_order), mean, na.rm = TRUE))

mouse_matrix <-
  mouse_avg %>%
  column_to_rownames("cell_type_refined") %>%
  as.matrix()

# Human

human_avg <-
  human@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(across(all_of(score_order), mean, na.rm = TRUE))

human_matrix <-
  human_avg %>%
  column_to_rownames("author_cell_type") %>%
  as.matrix()

mouse_scaled <- scale(mouse_matrix)
human_scaled <- scale(human_matrix)

# Figure 1: Cross-species correlation heatmap

cross_species_cor <- cor(
  t(mouse_scaled),
  t(human_scaled),
  method = "pearson"
)

pheatmap(
  cross_species_cor,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  border_color = "grey90",
  fontsize = 11,
  fontsize_number = 10,
  main = "Cross-species correlation of mechanotransduction signatures",
  filename = "../../results/07_snRNA/Figure1_CrossSpecies_CH.png",
  width = 9,
  height = 6
)


# meesenchymal and fibroblasts 
mouse_repair <- mouse_matrix[c(
  "Repair fibroblasts",
  "Repair-activated stromal cells",
  "Stromal progenitor-like cells",
  "Signaling stromal cells",
  "Proinflammatory mesenchymal cells",
  "Proliferating mesenchymal cells"
), ]

mouse_repair_scaled <- scale(mouse_repair)
human_scaled <- scale(human_matrix)

cross_species_cor <- cor(
  t(mouse_repair_scaled),
  t(human_scaled),
  method = "pearson"
)

pheatmap(
  cross_species_cor,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  border_color = "grey90",
  fontsize = 11,
  fontsize_number = 10,
  main = "Cross-species correlation of mechanotransduction signatures",
  filename = "../../results/07_snRNA/Figure1_CrossSpecies_CHR.png",
  width = 9,
  height = 6
)


