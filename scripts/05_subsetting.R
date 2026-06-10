library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(clustree)

combined <- readRDS("../data/combined_manual_annotated.rds")
levels(combined)

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

# Rerunning analysis on mesenchymal subset
mesenchymal <- NormalizeData(mesenchymal)
mesenchymal <- FindVariableFeatures(mesenchymal, selection.method = "vst", nfeatures = 2000)
mesenchymal <- ScaleData(mesenchymal)
mesenchymal <- RunPCA(mesenchymal, features = VariableFeatures(mesenchymal))


DimPlot(mesenchymal, reduction = "pca")
ElbowPlot(mesenchymal, ndims = 50)
pcs_to_use <- 1:20

pcs_to_use <- 1:20

mesenchymal <- FindNeighbors(mesenchymal, dims = pcs_to_use)

mesenchymal <- FindClusters(
  mesenchymal,
  resolution = seq(0.1, 0.8, 0.1),
  algorithm = 4,
  random.seed = 123
)

sapply(seq(0.1, 0.8, 0.1), function(res) {
  length(unique(mesenchymal[[paste0("RNA_snn_res.", res)]][,1]))
})

mesenchymal <- RunUMAP(
  mesenchymal,
  dims = pcs_to_use,
  seed.use = 123
)

resolutions <- seq(0.1, 0.8, 0.1)

umap_list <- lapply(resolutions, function(res) {
  DimPlot(
    mesenchymal,
    group.by = paste0("RNA_snn_res.", res),
    label = TRUE
  ) +
    ggtitle(paste("Resolution", res))
})

mesenchymal_umap <- wrap_plots(umap_list, ncol = 3)

mesenchymal_umap

ggsave(
  "../figures/mesenchymal_umap_resolutions.png",
  mesenchymal_umap,
  width = 18,
  height = 12,
  dpi = 300
)

clust <- clustree(
  mesenchymal,
  prefix = "RNA_snn_res."
)

clust

ggsave(
  "../figures/mesenchymal_cluster_tree.png",
  clust,
  width = 18,
  height = 12,
  dpi = 300
)

# 0.4 resolution 
Idents(mesenchymal) <- "RNA_snn_res.0.4"

# Store cluster assignments in metadata
mesenchymal$mesenchymal_cluster <- Idents(mesenchymal)

# Plot chosen clustering
umap_mesenchymal_res04 <- DimPlot(
  mesenchymal,
  reduction = "umap",
  group.by = "mesenchymal_cluster",
  label = TRUE
)

umap_mesenchymal_res04

# Save figure
ggsave(
  "../figures/mesenchymal_umap_res0.4.png",
  umap_mesenchymal_res04,
  width = 8,
  height = 6,
  dpi = 300
)

# Save object
saveRDS(
  mesenchymal,
  "../data/mesenchymal_processed_res0.4.rds"
)

set.seed(123)
# Find the marker genes for this subset
markers_mesenchymal <- FindAllMarkers(
  mesenchymal,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


saveRDS(
  markers_mesenchymal,
  "../data/markers_mesenchymal_res0.4.rds"
)

write.csv(
  markers_mesenchymal,
  "../data/markers_mesenchymal_res0.4.csv",
  row.names = FALSE
)

# top markers
top_markers_mesenchymal <- markers_mesenchymal %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15)

View(top_markers_mesenchymal)

write.csv(
  top_markers_mesenchymal,
  "../data/top15_markers_mesenchymal_res0.4.csv",
  row.names = FALSE
)

# umap by condition
m_umap_condition <- DimPlot(
  mesenchymal,
  group.by = "RNA_snn_res.0.4",
  split.by = "condition",
  ncol = 2,
  label = TRUE
)

ggsave(
  "../figures/mesenchymal_umap_condition.png",
  m_umap_condition,
  width = 12,
  height = 8,
  dpi = 300
)

DimPlot(
  mesenchymal,
  group.by = "cell_type_manual",
  label = TRUE
)

DimPlot(
  mesenchymal,
  group.by = "RNA_snn_res.0.4",
  label = TRUE
)

DimPlot(
  mesenchymal,
  group.by = "RNA_snn_res.0.4",
  split.by = "condition",
  ncol = 2,
  label = TRUE
)

# Checking for more associations
FindMarkers(
  mesenchymal,
  ident.1 = "5",
  ident.2 = "10"
)

DotPlot(
  mesenchymal,
  features = c(
    "Saa3", "Ereg", "Hp", "Serpina3g", "Tnn", "Adamts14",
    "Ccna2", "Plk1", "Kif11", "Kif2c", "Dlgap5", "Hmmr"
  ),
  group.by = "mesenchymal_cluster"
) + RotatedAxis()


FeaturePlot(
  mesenchymal,
  features = c(
    "Lyz2",
    "Tyrobp",
    "Csf1r",
    "Adgre1",
    "Cd68",
    "Fcgr3"
  )
)

VlnPlot(
  mesenchymal,
  features = c(
    "Lyz2",
    "Tyrobp",
    "Csf1r"
  ),
  group.by = "RNA_snn_res.0.4"
)
# Adding in final annotations 
library(Seurat)
library(dplyr)
library(plyr)

# Make sure identities are set to the chosen mesenchymal clustering
Idents(mesenchymal) <- "RNA_snn_res.0.4"

# Remove cluster 11
mesenchymal_clean <- subset(
  mesenchymal,
  idents = "11",
  invert = TRUE
)

# Reset identities after subsetting
Idents(mesenchymal_clean) <- "RNA_snn_res.0.4"

# Add refined labels
mesenchymal_clean$cell_type_refined <- plyr::mapvalues(
  x = as.character(Idents(mesenchymal_clean)),
  from = c("1","2","3","4","5","6","7","8","9","10"),
  to = c(
    "Repair fibroblasts",
    "ECM-remodelling tenocytes",
    "Homeostatic fibroblasts",
    "Fibrochondrocyte-like tenocytes",
    "Repair tenocytes",
    "Regulatory/signaling tenocytes",
    "Mature tenocytes",
    "Proinflammatory tenocytes",
    "Activated tenocytes",
    "Proliferating tenocytes"
  )
)

# Check labels
table(mesenchymal_clean$cell_type_refined)

# Plot refined labels
umap_refined <- DimPlot(
  mesenchymal_clean,
  reduction = "umap",
  group.by = "cell_type_refined",
  label = TRUE
)

umap_refined

# Save figure
ggsave(
  "../figures/mesenchymal_umap_refined_labels.png",
  umap_refined,
  width = 10,
  height = 7,
  dpi = 300
)

# Save final refined object
saveRDS(
  mesenchymal_clean,
  "../data/mesenchymal_refined_annotated.rds"
)

# Rename "Regulatory/signaling tenocytes" to "Signaling tenocytes"

mesenchymal_clean$cell_type_refined <- as.character(mesenchymal_clean$cell_type_refined)

mesenchymal_clean$cell_type_refined[
  mesenchymal_clean$cell_type_refined == "Regulatory/signaling tenocytes"
] <- "Signaling tenocytes"

mesenchymal_clean$cell_type_refined <- factor(mesenchymal_clean$cell_type_refined)

# Check
levels(mesenchymal_clean$cell_type_refined)
table(mesenchymal_clean$cell_type_refined)

# Replot
replot <- DimPlot(
  mesenchymal_clean,
  reduction = "umap",
  group.by = "cell_type_refined",
  label = TRUE
)

# Save figure
ggsave(
  "../figures/mesenchymal_umap_refined_labels.png",
  replot,
  width = 12,
  height = 8,
  dpi = 300
)

# Save updated object
saveRDS(
  mesenchymal_clean,
  "../data/mesenchymal_refined_annotated.rds"
)


DimPlot(mesenchymal_clean, group.by = "condition")
DimPlot(mesenchymal_clean, split.by = "condition", ncol = 2)
