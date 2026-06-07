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
