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
mesenchymal <- readRDS("../../data/mesenchymal_refined_annotated.rds")
gene_sets <- readRDS("../../data/all_gene_sets_refined.rds")

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

saveRDS(mesenchymal, "../../data/mesenchymal_geneset_scored_refined.rds")

# Checkpoint
mesenchymal <- readRDS("../../data/mesenchymal_geneset_scored_refined.rds")

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
  filename = "../../results/02_figures_refined/Figure1_Mechanotransduction_Heatmap.png",
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
  filename = "../../results/02_figures_refined/Figure2_FibroblastActivation_Heatmap.png",
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
  filename = "../../results/02_figures_refined/Figure3_ECM_Inflammation_Heatmap.png",
  width = 9,
  height = 6
)

# Condition-specific UCell heatmaps

condition_order <- c("WT", "I1D", "I7D", "I30D")

# Average UCell scores by condition and cell type
condition_celltype_scores <- mesenchymal@meta.data %>%
  group_by(condition, cell_type_refined) %>%
  summarise(
    across(all_of(score_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    cell_type_refined = factor(cell_type_refined, levels = celltype_order)
  ) %>%
  arrange(condition, cell_type_refined)


# Pretty column labels
score_label_map <- c(
  "Integrin_FAK_UCell" = "Integrin–FAK",
  "RhoA_ROCK_UCell" = "RhoA–ROCK",
  "YAP_TEAD_UCell" = "YAP/TEAD",
  "HALLMARK_TGF_BETA_SIGNALING_UCell" = "TGFβ Signaling",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell" = "SMAD Signaling",
  "GOBP_FIBROBLAST_ACTIVATION_UCell" = "Fibroblast Activation",
  "Myofibroblast_UCell" = "Myofibroblast",
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell" = "ECM Organization",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell" = "Inflammation"
)


# Figure 4: One heatmap per condition

score_order <- names(score_label_map)

mech_cols <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell"
)

fib_cols <- c(
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell"
)

ecm_cols <- c(
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

for (cond in condition_order) {
  
  cond_df <- condition_celltype_scores %>%
    filter(condition == cond) %>%
    select(cell_type_refined, all_of(score_order)) %>%
    column_to_rownames("cell_type_refined")
  
  cond_mat <- as.matrix(cond_df)
  
  row_order <- intersect(celltype_order, rownames(cond_mat))
  cond_mat <- cond_mat[row_order, , drop = FALSE]
  
  # Mechanotransduction / TGFβ / SMAD
  mat_mech <- cond_mat[, mech_cols, drop = FALSE]
  colnames(mat_mech) <- score_label_map[colnames(mat_mech)]
  
  pheatmap(
    mat_mech,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 9,
    angle_col = 45,
    main = paste("Mechanotransduction and TGFβ/SMAD -", cond),
    border_color = NA,
    filename = paste0("../../results/02_figures_refined/Figure4a_Mech_TGFB_SMAD_", cond, ".png"),
    width = 10,
    height = 6
  )
  
  # Fibroblast activation / Myofibroblast
  mat_fib <- cond_mat[, fib_cols, drop = FALSE]
  colnames(mat_fib) <- score_label_map[colnames(mat_fib)]
  
  pheatmap(
    mat_fib,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    main = paste("Fibroblast Activation -", cond),
    border_color = NA,
    filename = paste0("../../results/02_figures_refined/Figure4b_Fibroblast_", cond, ".png"),
    width = 10,
    height = 6
  )
  
  # ECM / Inflammation
  mat_ecm <- cond_mat[, ecm_cols, drop = FALSE]
  colnames(mat_ecm) <- score_label_map[colnames(mat_ecm)]
  
  pheatmap(
    mat_ecm,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    main = paste("ECM Organization and Inflammation -", cond),
    border_color = NA,
    filename = paste0("../../results/02_figures_refined/Figure4c_ECM_Inflammation_", cond, ".png"),
    width = 10,
    height = 6
  )
}

# Figure 5
mech_score_order <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell"
)

mech_label_map <- c(
  "Integrin_FAK_UCell" = "Integrin–FAK",
  "RhoA_ROCK_UCell" = "RhoA–ROCK",
  "YAP_TEAD_UCell" = "YAP/TEAD",
  "HALLMARK_TGF_BETA_SIGNALING_UCell" = "TGFβ Signaling",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell" = "SMAD Signaling"
)

mech_condition_scores <- mesenchymal@meta.data %>%
  group_by(condition) %>%
  summarise(
    across(all_of(mech_score_order), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(condition = factor(condition, levels = condition_order)) %>%
  arrange(condition)

mech_line_df <- mech_condition_scores %>%
  pivot_longer(
    cols = all_of(mech_score_order),
    names_to = "Signature",
    values_to = "Mean_Score"
  ) %>%
  mutate(
    Signature = recode(Signature, !!!mech_label_map),
    condition = factor(condition, levels = condition_order)
  )

ggplot(
  mech_line_df,
  aes(x = condition, y = Mean_Score, group = Signature, color = Signature)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "Mechanotransduction Signature Dynamics Across Healing Time",
    x = "Condition",
    y = "Mean UCell Score",
    color = "Signature"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "../../results/02_figures_refined/Figure5_Mechanotransduction_Dynamics_Lineplot.png",
  last_plot(),
  width = 10,
  height = 6,
  dpi = 300
)

VlnPlot(
  mesenchymal,
  features = "YAP_TEAD_UCell",
  group.by = "cell_type_refined",
  pt.size = 0
) +
  RotatedAxis()
