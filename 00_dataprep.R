library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)

# Extract data
# unzip("data/GSM8767956_WT_matrix.zip", exdir = "data/")
# unzip("data/GSM8767957_I1D_matrix.zip", exdir = "data/")
# unzip("data/GSM8767958_I7D_matrix.zip", exdir = "data/")
# unzip("data/GSM8767959_I30D_matrix.zip", exdir = "data/")

# Read in the data 
wt_counts <- Read10X(data.dir = "data/WT_matrix/")
i1d_counts <- Read10X(data.dir = "data/I1D_matrix")
i7d_counts <- Read10X(data.dir = "data/I7D_matrix")
i30d_counts <- Read10X(data.dir = "data/I30D_matrix")

# Convert into Seurat objects
wt <- CreateSeuratObject(counts = wt_counts, project = "WT")
i1d <- CreateSeuratObject(counts = i1d_counts, project = "I1D")
i7d <- CreateSeuratObject(counts = i7d_counts, project = "I7D")
i30d <- CreateSeuratObject(counts = i30d_counts, project = "I30D")

set.seed(123)

# Add metadata for tracking
wt$condition <- "WT"
i1d$condition <- "I1D"
i7d$condition <- "I7D"
i30d$condition <- "I30D"

# preprocessing for WT 
wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^mt-")

VlnPlot(wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

wt <- subset(
  wt,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5000 &
    percent.mt < 10
)

wt <- NormalizeData(wt)

wt <- FindVariableFeatures(wt, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(wt), 10)

# plot variable features with and without labels
plot1 <- FeatureScatter(wt,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")

plot2 <- FeatureScatter(wt,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
plot1 + plot2

all.genes <- rownames(wt)
wt <- ScaleData(wt, features = all.genes)

wt <- RunPCA(wt, features = VariableFeatures(object = wt))

VizDimLoadings(wt, dims = 1:2, reduction = "pca")
DimPlot(wt, reduction = "pca") + NoLegend()

# Doublet Finder
# 1. Find best pK
sweep.res <- paramSweep(wt, PCs = 1:20)
sweep.stats <- summarizeSweep(sweep.res)
bcmvn <- find.pK(sweep.stats)

View(bcmvn)

n_cells <- ncol(wt)
nExp <- round(0.008 * (n_cells / 1000) * n_cells)
nExp

wt <- doubletFinder(
  wt,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.15,
  nExp = nExp
)

wt_clean <- subset(
  wt,
  subset = DF.classifications_0.25_0.15_185 == "Singlet"
)

ncol(wt)
ncol(wt_clean)

wt <- RunUMAP(wt, dims = 1:30)

DimPlot(
  wt,
  group.by = "DF.classifications_0.25_0.15_185"
)
