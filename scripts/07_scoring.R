library(Seurat)
library(UCell)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load mesenchymal object
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")

# Load curated gene sets
gene_sets <- readRDS("../data/mechanotransduction_gene_sets.rds")
gene_set_groups <- readRDS("../data/mechanotransduction_gene_set_groups.rds")

# Run UCell scoring
mesenchymal <- AddModuleScore_UCell(
  mesenchymal,
  features = gene_sets
)

# Check new score columns
score_cols <- grep("_UCell$", colnames(mesenchymal@meta.data), value = TRUE)
score_cols

# Save scored object
saveRDS(
  mesenchymal,
  "../data/mesenchymal_gene_set_scored.rds"
)

# Supporting visualizations
FeaturePlot(
  mesenchymal,
  features = c(
    "YAP_TAZ_CURATED_UCell",
    "HALLMARK_TGF_BETA_SIGNALING_UCell",
    "MYOFIBROBLAST_CURATED_UCell"
  ),
  ncol = 3
)

VlnPlot(
  mesenchymal,
  features = c(
    "YAP_TAZ_CURATED_UCell",
    "HALLMARK_TGF_BETA_SIGNALING_UCell",
    "MYOFIBROBLAST_CURATED_UCell"
  ),
  group.by = "cell_type_manual",
  pt.size = 0
)

library(dplyr)

mesenchymal@meta.data %>%
  group_by(cell_type_manual) %>%
  summarise(
    YAP = mean(YAP_TAZ_CURATED_UCell),
    TGFB = mean(HALLMARK_TGF_BETA_SIGNALING_UCell),
    MYO = mean(MYOFIBROBLAST_CURATED_UCell)
  ) %>%
  arrange(desc(YAP))

# YAP/TAZ figure by condition
library(scCustomize)
library(stringr)

plot_feature_by_condition <- function(seurat_object, feature, title = NULL,
                                      conditions = c("WT", "I1D", "I7D", "I30D"),
                                      ncol = 2, pt.size = 0.5) {
  
  clean_title <- if (!is.null(title)) {
    title
  } else {
    feature %>%
      str_remove("_UCell$") %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  }
  
  plots <- lapply(conditions, function(cond) {
    FeaturePlot_scCustom(
      seurat_object = subset(seurat_object, subset = condition == cond),
      features = feature,
      pt.size = pt.size
    ) +
      ggtitle(cond) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_blank(),
        strip.background = element_blank()
      )
  })
  
  wrap_plots(plots, ncol = ncol) +
    plot_annotation(title = clean_title) &
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
}

# Usage
plot_feature_by_condition(mesenchymal, "HALLMARK_TGF_BETA_SIGNALING_UCell",
                          title = "TGF-β Signaling UCell Score by Condition")

plot_feature_by_condition(mesenchymal, "MYOFIBROBLAST_CURATED_UCell",
                          title = "Myofibroblast UCell Score by Condition")

plot_feature_by_condition(mesenchymal, "YAP_TAZ_CURATED_UCell",
                          title = "YAP/TAZ UCell Score by Condition")

scores <- mesenchymal@meta.data %>%
       select(
             YAP_TAZ_CURATED_UCell,
             HALLMARK_TGF_BETA_SIGNALING_UCell,
             MYOFIBROBLAST_CURATED_UCell
         )

cor(scores)


avg_scores <- mesenchymal@meta.data %>%
       group_by(cell_type_manual) %>%
       summarise(
           across(ends_with("_UCell"), mean)
       )
library(pheatmap)
library(tibble)
library(dplyr)

# Convert avg_scores to matrix
heatmap_mat <- avg_scores %>%
  column_to_rownames("cell_type_manual") %>%
  as.matrix()

# Scale each gene set column
heatmap_scaled <- scale(heatmap_mat)

# Shorten column names
colnames(heatmap_scaled) <- colnames(heatmap_scaled) %>%
  gsub("_UCell", "", .) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("REACTOME_", "", .) %>%
  gsub("GOBP_", "", .) %>%
  gsub("_", " ", .)


pheatmap(
  heatmap_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 7,
  angle_col = 90,
  main = "Average UCell Scores by Cell Type",
  border_color = NA
)

# Based on 5 biological groups
for (group_name in names(gene_set_groups)) {
  
  group_score_cols <- paste0(
    gene_set_groups[[group_name]],
    "_UCell"
  )
  
  group_score_cols <- intersect(
    group_score_cols,
    colnames(mesenchymal@meta.data)
  )
  
  mesenchymal[[paste0(group_name, "_Score")]] <- rowMeans(
    mesenchymal@meta.data[, group_score_cols, drop = FALSE],
    na.rm = TRUE
  )
}

library(dplyr)

group_scores <- mesenchymal@meta.data %>%
  group_by(cell_type_manual) %>%
  summarise(
    Mechanotransduction = mean(Mechanotransduction_Score),
    TGFB_SMAD = mean(TGFB_SMAD_Score),
    ECM_Fibrosis = mean(ECM_Fibrosis_Score),
    Inflammation = mean(Inflammation_Score),
    Repair_Tendon_State = mean(Repair_Tendon_State_Score)
  )

group_scores

library(pheatmap)
library(tibble)

group_mat <- group_scores %>%
  column_to_rownames("cell_type_manual") %>%
  as.matrix()

group_scaled <- scale(group_mat)

desired_order <- c(
  "Homeostatic fibroblasts",
  "Mature tenocytes",
  "Proinflammatory tenocytes",
  "Activated tenocytes",
  "Repair fibroblasts",
  "Proliferating tenocytes",
  "ECM-remodelling tenocytes",
  "Fibrochondrocyte-like tenocytes",
  "Signaling tenocytes"
)

group_mat <- group_mat[desired_order, ]

pheatmap(
  scale(group_mat),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 11,
  fontsize_col = 11,
  main = "Biological Program Activity by Cell Type",
  border_color = NA
)


# By Condition/Timepoint
library(dplyr)
library(tibble)
library(pheatmap)

condition_group_scores <- mesenchymal@meta.data %>%
  group_by(condition) %>%
  summarise(
    Mechanotransduction = mean(Mechanotransduction_Score),
    TGFB_SMAD = mean(TGFB_SMAD_Score),
    ECM_Fibrosis = mean(ECM_Fibrosis_Score),
    Inflammation = mean(Inflammation_Score),
    Repair_Tendon_State = mean(Repair_Tendon_State_Score),
    .groups = "drop"
  )

condition_order <- c("WT", "I1D", "I7D", "I30D")

condition_group_scores <- condition_group_scores %>%
  mutate(condition = factor(condition, levels = condition_order)) %>%
  arrange(condition)

condition_mat <- condition_group_scores %>%
  column_to_rownames("condition") %>%
  as.matrix()

pheatmap(
  condition_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Biological Program Activity by Condition",
  border_color = NA
)

DotPlot(
  mesenchymal,
  features = c(
    "YAP_TAZ_CURATED_UCell",
    "TGFB_SMAD_Score",
    "ECM_Fibrosis_Score",
    "Inflammation_Score"
  ),
  group.by = "cell_type_manual"
)


library(dplyr)
library(tibble)
library(pheatmap)

condition_order <- c("WT", "I1D", "I7D", "I30D")

celltype_order <- c(
  "Homeostatic fibroblasts",
  "Mature tenocytes",
  "Proinflammatory tenocytes",
  "Activated tenocytes",
  "Repair fibroblasts",
  "Proliferating tenocytes",
  "ECM-remodelling tenocytes",
  "Fibrochondrocyte-like tenocytes",
  "Signaling tenocytes"
)

group_scores_condition_ordered <- group_scores_condition %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    cell_type_manual = factor(cell_type_manual, levels = celltype_order),
    cell_condition = paste(condition, cell_type_manual, sep = " | ")
  ) %>%
  arrange(condition, cell_type_manual)

group_condition_mat <- group_scores_condition_ordered %>%
  select(
    cell_condition,
    Mechanotransduction,
    TGFB_SMAD,
    ECM_Fibrosis,
    Inflammation,
    Repair_Tendon_State
  ) %>%
  column_to_rownames("cell_condition") %>%
  as.matrix()

group_condition_scaled <- scale(group_condition_mat)

pheatmap(
  group_condition_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 11,
  main = "Biological Program Activity by Cell Type and Condition",
  border_color = NA
)

library(dplyr)
library(tidyr)
library(ggplot2)

line_df <- group_scores_condition %>%
  mutate(
    condition = factor(condition, levels = c("WT", "I1D", "I7D", "I30D"))
  ) %>%
  pivot_longer(
    cols = c(
      Mechanotransduction,
      TGFB_SMAD,
      ECM_Fibrosis,
      Inflammation,
      Repair_Tendon_State
    ),
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
  facet_wrap(~ cell_type_manual, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Biological Program Dynamics Across Healing Time",
    x = "Condition",
    y = "Mean UCell Group Score",
    color = "Biological Program"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
