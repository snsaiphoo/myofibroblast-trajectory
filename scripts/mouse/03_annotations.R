library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

set.seed(123)

data_dir <- "../../data"
fig_dir <- "../results/01_figures_all"

# Load processed Seurat object
combined <- readRDS(file.path(data_dir, "combined_processed_resolutions.rds"))

# Set chosen clustering resolution
Idents(combined) <- "RNA_snn_res.0.5"

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(combined)

# Load mouse reference
ref <- MouseRNAseqData()

# Run SingleR
singleR_results <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

# Add SingleR labels to Seurat metadata
combined$SingleR_label <- singleR_results$labels
combined$SingleR_pruned <- singleR_results$pruned.labels

# Plot SingleR annotations
singleR_umap <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "SingleR_label",
  label = TRUE,
  repel = TRUE
)

singleR_pruned_umap <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "SingleR_pruned",
  label = TRUE,
  repel = TRUE
)

singleR_umap
singleR_pruned_umap

# Compare SingleR labels with Seurat clusters
singleR_cluster_table <- table(
  combined$RNA_snn_res.0.5,
  combined$SingleR_label
)

singleR_cluster_table

# Save SingleR results
saveRDS(
  singleR_results,
  file.path(data_dir, "singleR_results_res0.5.rds")
)

write.csv(
  as.data.frame(singleR_cluster_table),
  file.path(data_dir, "singleR_cluster_table_res0.5.csv"),
  row.names = FALSE
)

ggsave(
  file.path(fig_dir, "singleR_labels_umap.png"),
  singleR_umap,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "singleR_pruned_umap.png"),
  singleR_pruned_umap,
  width = 10,
  height = 8,
  dpi = 300
)

# Save Seurat object with SingleR labels
saveRDS(
  combined,
  file.path(data_dir, "combined_singler_annotated.rds")
)