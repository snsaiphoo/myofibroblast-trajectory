library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(clustree)

set.seed(123)

combined <- readRDS("../data/combined_singlets_counts_only.rds")

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))

DimPlot(combined, reduction = "pca")
ElbowPlot(combined, ndims = 50)

pcs_to_use <- 1:20

combined <- FindNeighbors(
  combined,
  dims = pcs_to_use
)

# Leiden algorithm is better suited for similar cell types
combined <- FindClusters(
  combined,
  resolution = seq(0.1, 0.8, 0.1),
  algorithm = 4
)

saveRDS(combined, "../data/checkpoint_03_clusters.rds")

sapply(seq(0.1, 0.8, 0.1), function(res) {
  length(unique(combined[[paste0("RNA_snn_res.", res)]][,1]))
})

saveRDS(combined, "../data/combined_processed_resolutions_REPRODUCED.rds")

combined <- RunUMAP(
  combined,
  dims = pcs_to_use
)

# Compare UMAPs to choose resolution
resolutions <- seq(0.1, 0.8, 0.1)

umap_list <- lapply(resolutions, function(res) {
  DimPlot(
    combined,
    group.by = paste0("RNA_snn_res.", res),
    label = TRUE
  ) +
    ggtitle(paste("Resolution", res))
})

combined_umap <- wrap_plots(umap_list, ncol = 3)

combined_umap

ggsave(
  "../figures/umap_resolutions.png",
  combined_umap,
  width = 18,
  height = 12,
  dpi = 300
)

# Clustertree to determine resolution
clust <- clustree(
  combined,
  prefix = "RNA_snn_res."
)

clust

ggsave(
  "../figures/cluster_tree.png",
  clust,
  width = 18,
  height = 12,
  dpi = 300
)

# Choose a resolution
Idents(combined) <- "RNA_snn_res.0.5"

umap_full <- DimPlot(
  combined,
  reduction = "umap",
  label = TRUE
)

umap_condition <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "condition"
)

umap_full
umap_condition

saveRDS(combined, "../data/combined_processed_resolutions_REPRODUCED_umap.rds")

ggsave("../figures/umap_full.png", umap_full, width = 8, height = 6, dpi = 300)
ggsave("../figures/umap_condition.png", umap_condition, width = 8, height = 6, dpi = 300)

saveRDS(combined, "../data/combined_processed_resolutions.rds")

# Find all markers 
Idents(combined) <- "RNA_snn_res.0.5"

markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

saveRDS(
  markers,
  "../data/markers_res0.5.rds"
)

write.csv(
  markers,
  "../data/markers_res0.5.csv",
  row.names = FALSE
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

View(top_markers)

write.csv(markers, "../figures/cluster_markers_res0.5.csv", row.names = FALSE)
write.csv(top_markers, "../figures/top10_markers_res0.5.csv", row.names = FALSE)

# Tenocyte/Fibroblast
FeaturePlot(
  combined,
  features = c("Col1a1", "Dcn", "Tnmd", "Scx")
)

FeaturePlot(
  combined,
  features = c("Thbs4", "Comp", "Fmod", "Mkx")
)

DotPlot(
  combined,
  features = c(
    "Scx",
    "Tnmd",
    "Thbs4",
    "Comp"
  )
) + RotatedAxis()

