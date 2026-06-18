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

# load data
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")
gene_sets <- readRDS("../data/mechanotransduction_gene_sets.rds")
gene_set_groups <- readRDS("../data/mechanotransduction_gene_set_groups.rds")

# run ucell scoring
mesenchymal <- AddModuleScore_UCell(
  mesenchymal,
  features = gene_sets
)

# check score columns
score_cols <- grep("_UCell$", colnames(mesenchymal@meta.data), value = TRUE)
score_cols

# save scored object
saveRDS(mesenchymal, "../data/mesenchymal_gene_set_scored.rds")

# compute group scores by averaging ucell scores within each biological group
for (group_name in names(gene_set_groups)) {
  group_score_cols <- paste0(gene_set_groups[[group_name]], "_UCell")
  group_score_cols <- intersect(group_score_cols, colnames(mesenchymal@meta.data))
  mesenchymal[[paste0(group_name, "_Score")]] <- rowMeans(
    mesenchymal@meta.data[, group_score_cols, drop = FALSE],
    na.rm = TRUE
  )
}

# cell type order
celltype_order <- c(
  "Mature tenocytes",
  "Stromal progenitor-like cells",
  "Signaling stromal cells",
  "Fibrochondrocyte-like tenocytes",
  "Activated tenocytes",
  "Proinflammatory mesenchymal cells",
  "Proliferating mesenchymal cells",
  "Repair fibroblasts",
  "Repair-activated stromal cells",
  "ECM-remodelling tenocytes"
)

condition_order <- c("WT", "I1D", "I7D", "I30D")

# figure 1: heatmap of average ucell scores by cell type

avg_scores <- mesenchymal@meta.data %>%
  group_by(cell_type_refined) %>%
  dplyr::summarise(dplyr::across(ends_with("_UCell"), mean), .groups = "drop")

heatmap_mat <- avg_scores %>%
  column_to_rownames("cell_type_refined") %>%
  as.matrix()

heatmap_mat <- heatmap_mat[celltype_order, ]

heatmap_scaled <- scale(heatmap_mat)

colnames(heatmap_scaled) <- colnames(heatmap_scaled) %>%
  gsub("_UCell", "", .) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("REACTOME_", "", .) %>%
  gsub("GOBP_", "", .) %>%
  gsub("_", " ", .)

pheatmap(
  heatmap_scaled,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 7,
  angle_col = 90,
  main = "Average UCell Scores by Cell Type",
  border_color = NA,
  filename = "../figures/heatmap_ucell_by_celltype.png",
  width = 12,
  height = 6
)

# figure 2: heatmap of biological program activity by cell type x condition

# heatmap of biological program activity by condition
condition_scores <- mesenchymal@meta.data %>%
  group_by(condition) %>%
  dplyr::summarise(
    Mechanotransduction = mean(Mechanotransduction_Score),
    TGFB_SMAD = mean(TGFB_SMAD_Score),
    ECM_Remodeling = mean(ECM_Remodeling_Score),
    Inflammation  = mean(Inflammation_Score),
    Repair_Tendon_State = mean(Repair_Tendon_State_Score),
    .groups = "drop"
  ) %>%
  mutate(condition = factor(condition, levels = condition_order)) %>%
  arrange(condition)

condition_mat <- condition_scores %>%
  column_to_rownames("condition") %>%
  as.matrix()

pheatmap(
  condition_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 12,
  cellheight = 40,
  cellwidth = 80,
  angle_col = 45,
  main = "Biological Program Activity by Condition",
  border_color = NA,
  filename = "../figures/heatmap_programs_by_condition.png",
  width = 8,
  height = 5
)

# figure 3: line plot of biological program dynamics across healing time 

# compute group scores by cell type and condition
group_scores_condition <- mesenchymal@meta.data %>%
  group_by(cell_type_refined, condition) %>%
  dplyr::summarise(
    Mechanotransduction = mean(Mechanotransduction_Score),
    TGFB_SMAD  = mean(TGFB_SMAD_Score),
    ECM_Remodeling = mean(ECM_Remodeling_Score),
    Inflammation = mean(Inflammation_Score),
    Repair_Tendon_State = mean(Repair_Tendon_State_Score),
    .groups = "drop"
  ) %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    cell_type_refined = factor(cell_type_refined, levels = celltype_order)
  ) %>%
  arrange(condition, cell_type_refined)

line_df <- group_scores_condition %>%
  pivot_longer(
    cols = c(Mechanotransduction, TGFB_SMAD, ECM_Remodeling,
             Inflammation, Repair_Tendon_State),
    names_to = "Biological_Program",
    values_to = "Mean_Score"
  )

ggplot(
  line_df,
  aes(
    x = condition,
    y = Mean_Score,
    group = Biological_Program,
    color = Biological_Program
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ cell_type_refined, scales = "free_y", ncol = 3) +
  theme_classic() +
  labs(
    title  = "Biological Program Dynamics Across Healing Time",
    x = "Condition",
    y = "Mean UCell Group Score",
    color = "Biological Program"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "../figures/lineplot_program_dynamics.png",
  last_plot(),
  width = 14, height = 10, dpi = 300
)

# figure 4: heatmap per condition 