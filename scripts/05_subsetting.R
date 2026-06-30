library(Seurat)
library(dplyr)
library(plyr)
library(patchwork)
library(ggplot2)
library(clustree)

# Subsetting the mesenchymal populations
combined <- readRDS("../data/combined_manual_annotated.rds")

mesenchymal <- subset(
  combined,
  idents = c(
    "ECM-remodelling tenocytes",
    "Repair fibroblasts",
    "Homeostatic fibroblasts",
    "Proliferating tenocytes",
    "Fibrochondrocyte-like tenocytes",
    "Signaling tenocytes",
    "Mature tenocytes",
    "Proinflammatory tenocytes",
    "Activated tenocytes",
    "Stromal cells"
  )
)

saveRDS(mesenchymal, "../data/mesenchymal_subset.rds")

# reprocess

mesenchymal <- NormalizeData(mesenchymal)
mesenchymal <- FindVariableFeatures(mesenchymal, selection.method = "vst", nfeatures = 2000)
mesenchymal <- ScaleData(mesenchymal)
mesenchymal <- RunPCA(mesenchymal, features = VariableFeatures(mesenchymal))

# Select number of PCs
ElbowPlot(mesenchymal, ndims = 50)
pcs_to_use <- 1:20

mesenchymal <- FindNeighbors(mesenchymal, dims = pcs_to_use)
mesenchymal <- FindClusters(
  mesenchymal,
  resolution = seq(0.1, 0.8, 0.1),
  algorithm = 4,
  random.seed = 123
)
mesenchymal <- RunUMAP(mesenchymal, dims = pcs_to_use, seed.use = 123)

# Check cluster counts per resolution
sapply(seq(0.1, 0.8, 0.1), function(res) {
  length(unique(mesenchymal[[paste0("RNA_snn_res.", res)]][, 1]))
})

# select new resolution

# Plot all resolutions
umap_list <- lapply(seq(0.1, 0.8, 0.1), function(res) {
  DimPlot(mesenchymal, group.by = paste0("RNA_snn_res.", res), label = TRUE) +
    ggtitle(paste("Resolution", res)) +
    theme(legend.position = "none")
})

wrap_plots(umap_list, ncol = 3)

ggsave(
  "../figures/mesenchymal_umap_resolutions.png",
  last_plot(),
  width = 18, height = 12, dpi = 300
)

# Clustree for resolution stability
clustree(mesenchymal, prefix = "RNA_snn_res.")

ggsave(
  "../figures/mesenchymal_cluster_tree.png",
  last_plot(),
  width = 18, height = 12, dpi = 300
)

# Set chosen resolution
Idents(mesenchymal) <- "RNA_snn_res.0.4"
mesenchymal$mesenchymal_cluster <- Idents(mesenchymal)

DimPlot(mesenchymal, reduction = "umap", group.by = "mesenchymal_cluster", label = TRUE)

ggsave(
  "../figures/mesenchymal_umap_res0.4.png",
  last_plot(),
  width = 8, height = 6, dpi = 300
)

saveRDS(mesenchymal, "../data/mesenchymal_processed_res0.4.rds")

# marker gene identification
mesenchymal <- readRDS("../data/mesenchymal_processed_res0.4.rds")

set.seed(123)

markers_mesenchymal <- FindAllMarkers(
  mesenchymal,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.50
)

top15_markers <- markers_mesenchymal %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15)

write.csv(markers_mesenchymal, "../data/markers_mesenchymal_res0.4.csv", row.names = FALSE)
write.csv(top15_markers, "../data/top15_markers_mesenchymal_res0.4.csv", row.names = FALSE)
saveRDS(markers_mesenchymal, "../data/markers_mesenchymal_res0.4.rds")

mesenchymal <- readRDS("../data/mesenchymal_processed_res0.4.rds")

# Tenocyte identity markers across all clusters
VlnPlot(
  mesenchymal,
  features = c("Scx", "Tnmd", "Mkx", "Egr1"),
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
)

# Progenitor markers for Cluster 3
VlnPlot(
  mesenchymal,
  features = c("Cd55", "Procr", "Pdgfra", "Ly6a", "Tppp3", "Cd248"),
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
)

# refined annotation

Idents(mesenchymal) <- "RNA_snn_res.0.4"

# Remove Cluster 11 (nerve-associated stromal cells)
mesenchymal_clean <- subset(mesenchymal, idents = "11", invert = TRUE)
Idents(mesenchymal_clean) <- "RNA_snn_res.0.4"

# Assign refined labels
mesenchymal_clean$cell_type_refined <- plyr::mapvalues(
  x = as.character(Idents(mesenchymal_clean)),
  from = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  to = c(
    "Repair fibroblasts",
    "ECM-remodelling tenocytes",
    "Stromal progenitor-like cells",
    "Fibrochondrocyte-like tenocytes",
    "Repair-activated stromal cells",
    "Signaling stromal cells",
    "Mature tenocytes",
    "Proinflammatory mesenchymal cells",
    "Activated tenocytes",
    "Proliferating mesenchymal cells"
  )
)

Idents(mesenchymal_clean) <- "cell_type_refined"

# Check cell counts
table(mesenchymal_clean$cell_type_refined)

# final plots
# UMAP with refined labels
DimPlot(
  mesenchymal_clean,
  reduction = "umap",
  group.by = "cell_type_refined",
  label = TRUE,
)

ggsave(
  "../figures/mesenchymal_umap_refined_labels.png",
  last_plot(),
  width = 12, height = 8, dpi = 300
)

# UMAP split by condition
DimPlot(
  mesenchymal_clean,
  reduction = "umap",
  group.by = "cell_type_refined",  
  split.by = "condition",
  ncol = 2
)

ggsave(
  "../figures/mesenchymal_umap_refined_by_condition.png",
  last_plot(),
  width = 16, height = 8, dpi = 300
)

# Annotation validation DotPlot
DotPlot(
  mesenchymal_clean,
  features = c(
    "Tnmd", "Fmod", "Comp", "Chad",
    "Scx", "Mkx",
    "Acta2", "Postn", "Col1a1",
    "Sparcl1", "Cilp", "Aspn",
    "Acan", "Col11a1",
    "Cd55", "Procr",
    "Sfrp4", "Gdf10",
    "S100a8", "Cxcl5", "Il1b",
    "Prg4", "Timp1",
    "Ccna2", "Plk1"
  ),
  group.by = "cell_type_refined"
) + RotatedAxis()

ggsave(
  "../figures/mesenchymal_annotation_validation_dotplot.png",
  last_plot(),
  width = 18, height = 8, dpi = 300
)

saveRDS(mesenchymal_clean, "../data/mesenchymal_refined_annotated.rds")
