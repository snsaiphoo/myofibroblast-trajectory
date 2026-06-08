library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

combined <- readRDS("../data/combined_processed_resolutions.rds")
markers <- readRDS("../data/markers_res0.5.rds")

Idents(combined) <- "RNA_snn_res.0.5"

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15)

View(top_markers) 

write.csv(
  top_markers,
  "../data/top15_markers_res0.5_for_manual_annotation.csv",
  row.names = FALSE
)

DimPlot(
  combined,
  reduction = "umap",
  split.by = "condition",
  label = TRUE,
  ncol = 2
)


FeaturePlot(combined, features = c(
  "Pdgfra",
  "Pi16",
  "Cd34",
  "Dpp4",
  "Ly6a",
  "Postn",
  "Acta2",
  "Tagln"
))

FeaturePlot(
  combined,
  features = c(
    "Ctsk",
    "Bmp2",
    "Sox9",
    "Runx2"
  )
)

# manual annotation
# Rename clusters
combined <- RenameIdents(
  combined,
  "1" = "Tendon-resident macrophages",
  "2" = "Late-stage ECM-remodelling tenocytes",
  "3" = "FAP-like repair fibroblasts",
  "4" = "Homeostatic fibroblasts",
  "5" = "Proliferating tenocytes",
  "6" = "Fibrochondrocyte-like tenocytes",
  "7" = "Activated signaling tenocytes",
  "8" = "Inflammatory myeloid cells",
  "9" = "Mature tenocytes",
  "10" = "Proinflammatory tenocytes",
  "11" = "Vascular endothelial cells",
  "12" = "Early injury-activated tenocytes",
  "13" = "Muscle-associated myogenic cells",
  "14" = "Dendritic cells",
  "15" = "Pericytes / vascular smooth muscle cells",
  "16" = "Nerve-associated cells",
  "17" = "Neutrophils",
  "18" = "Mast cells"
)

# Save annotations to metadata
combined$cell_type_manual <- Idents(combined)

# UMAP with legend only
p_manual <- DimPlot(
  combined,
  group.by = "cell_type_manual"
)

ggsave(
  filename = "../figures/manual_annotation_umap.png",
  plot = p_manual,
  width = 10,
  height = 8,
  dpi = 300
)

# UMAP split by condition, 2x2 grid
p_manual_split <- DimPlot(
  combined,
  group.by = "cell_type_manual",
  split.by = "condition",
  ncol = 2
)

ggsave(
  filename = "../figures/manual_annotation_umap_by_condition.png",
  plot = p_manual_split,
  width = 14,
  height = 10,
  dpi = 300
)
