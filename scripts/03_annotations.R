library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(BiocParallel)
library(ggplot2)

set.seed(123)

# Load processed Seurat object
combined <- readRDS("../data/combined_processed_resolutions.rds")

# Set clustering resolution
Idents(combined) <- "RNA_snn_res.0.5"

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(combined)

# References
ref_mouse <- MouseRNAseqData()
ref_immune <- ImmGenData()

# Parallelization
param <- SnowParam(
  workers = 4,
  type = "SOCK"
)

# Run SingleR using both references
pred <- SingleR(
  test = sce,
  ref = list(
    MouseRNAseq = ref_mouse,
    ImmGen = ref_immune
  ),
  labels = list(
    ref_mouse$label.main,
    ref_immune$label.main
  ),
  BPPARAM = param
)

# Add labels to Seurat object
combined$cell_type_singleR <- pred$labels

# UMAP of SingleR labels
p_singleR <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "cell_type_singleR",
  label = TRUE,
  repel = TRUE,
  label.size = 4.5
) +
  ggtitle("SingleR Annotation") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

p_singleR

ggsave(
  "../figures/umap_singleR.png",
  p_singleR,
  width = 12,
  height = 9,
  dpi = 300
)

saveRDS(
  combined,
  "../data/combined_singler_multi_annotated.rds"
)
