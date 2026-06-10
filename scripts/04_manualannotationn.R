library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(123)
combined <- readRDS("../data/combined_processed_resolutions_REPRODUCED_umap.rds")
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
  "2" = "ECM-remodelling tenocytes",
  "3" = "Repair fibroblasts",
  "4" = "Homeostatic fibroblasts",
  "5" = "Proliferating tenocytes",
  "6" = "Fibrochondrocyte-like tenocytes",
  "7" = "Signaling tenocytes",
  "8" = "Inflammatory myeloid cells",
  "9" = "Mature tenocytes",
  "10" = "Proinflammatory tenocytes",
  "11" = "Vascular endothelial cells",
  "12" = "Activated tenocytes",
  "13" = "Muscle-associated myogenic cells",
  "14" = "Dendritic cells",
  "15" = "Mural cells",
  "16" = "Stromal cells",
  "17" = "Neutrophils",
  "18" = "Mast cells"
)

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

FeaturePlot(
  combined,
  features = c("Acta2", "Tagln", "Myl9", "Tpm2", "Cnn1", "Col1a1", "Col1a2"),
  ncol = 3
)

VlnPlot(
  combined,
  features = c("Acta2", "Tagln", "Myl9", "Tpm2", "Cnn1"),
  group.by = "cell_type_manual",
  pt.size = 0
)

celltype_markers <- c(
  "C1qa",     # macrophages
  "Dlk1",     # FAPs
  "Procr",    # stromal fibroblasts
  "Ccnb1",    # proliferating
  "Acan",     # fibrochondrocytes
  "Gdf10",    # activated tenocytes
  "Fmod",     # mature tenocytes
  "Cxcl5",    # inflammatory tenocytes
  "Emcn",     # endothelial
  "Flt3",     # dendritic
  "Myh11",    # mural/pericytes
  "A2m",      # stromal cells
  "S100a8",   # neutrophils
  "Tpsab1"    # mast cells
)

DotPlot(
  combined,
  features = celltype_markers,
  group.by = "cell_type_manual"
) +
  RotatedAxis()


Idents(combined) <- "cell_type_manual"

markers_10_vs_12 <- FindMarkers(
  combined,
  ident.1 = "Proinflammatory tenocytes",
  ident.2 = "Activated tenocytes"
)

head(markers_10_vs_12)

# cluster 10 and 12 should remain separated

saveRDS(
  combined,
  "../data/combined_manual_annotated.rds"
)
