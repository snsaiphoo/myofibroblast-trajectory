library(Seurat)

combined <- readRDS("../data/combined_processed_resolutions.rds")

Idents(combined) <- "RNA_snn_res.0.5"
DimPlot(combined, reduction = "umap", label = TRUE)

FeaturePlot(
  combined,
  features = c("Mest", "Dlk1", "H19", "Ptn", "Postn", "Plagl1"),
  cols = c("red")
)

DotPlot(
  combined,
  features = c(
    "Mest", "Dlk1", "H19", "Ptn", "Plagl1",
    "Postn", "Cthrc1", "Lrrc15", "Fn1", "Tnc",
    "Scx", "Tnmd", "Thbs4", "Comp"
  )
) + RotatedAxis

FeaturePlot(
  combined,
  features = c("Olfml2a", "Procr", "C7", "Cd55", "Gsn", "Pcolce2")
)

DotPlot(
  combined,
  features = c(
    "Col1a1", "Dcn", "Lum", "Fmod",
    "Scx", "Tnmd", "Thbs4", "Comp",
    "Postn", "Cthrc1", "Fn1", "Acta2", "Tagln"
  )
) + RotatedAxis()
