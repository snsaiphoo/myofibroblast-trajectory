# Initial violin plots
plot_qc <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  print(
  VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  ))
  
  return(seurat_obj)

}

# Filtering based off of the Violin Plots
filter_qc <- function(seurat_obj,
                      min_features,
                      max_features,
                      max_percent_mt) {
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^mt-"
  )
  
  seurat_obj <- subset(
    seurat_obj,
    subset =
      nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )
  
  return(seurat_obj)
}

# preprocess the pca 
preprocess_pca <- function(seurat_obj) {
  
  seurat_obj <- NormalizeData(seurat_obj)
  
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000
  )
  
  plot1 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )
  
  plot2 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
  
  print(plot1 + plot2)
  
  all.genes <- rownames(seurat_obj)
  
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(object = seurat_obj)
  )
  
  print(VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca"))
  
  print(DimPlot(seurat_obj, reduction = "pca") + NoLegend())
  
  return(seurat_obj)
}

# 1. Show elbow plot so you can choose PCs
plot_elbow <- function(seurat_obj, ndims = 50) {
  
  print(ElbowPlot(seurat_obj, ndims = ndims))
}

# 2. Run pK sweep and view best pK table
find_best_pk <- function(seurat_obj, pcs) {
  
  sweep.res <- paramSweep(seurat_obj, PCs = pcs)
  sweep.stats <- summarizeSweep(sweep.res)
  bcmvn <- find.pK(sweep.stats)
  
  View(bcmvn)
  
  return(bcmvn)
}

# 3. Run DoubletFinder
run_doubletfinder <- function(seurat_obj,
                              pcs,
                              pK,
                              pN = 0.25,
                              doublet_rate = 0.008) {
  
  n_cells <- ncol(seurat_obj)
  
  nExp <- round(
    doublet_rate * (n_cells / 1000) * n_cells
  )
  
  cat("Using pK =", pK, "\n")
  cat("Expected doublets =", nExp, "\n")
  
  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = pcs,
    pN = pN,
    pK = pK,
    nExp = nExp
  )
  
  return(seurat_obj)
}