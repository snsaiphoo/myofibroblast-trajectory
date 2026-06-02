source("functions/seuratprocess.R")

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
wt <- plot_qc(wt)
wt <- filter_qc(wt, 500, 5000, 10)
wt <- preprocess_pca(wt)

# Doublet Finder for WT
plot_elbow(wt, ndims = 50)
pcs_to_use <- 1:25
bcmvn <- find_best_pk(wt, pcs = pcs_to_use)
wt <- run_doubletfinder(wt, pcs = 1:25, pK = 0.15)

df_col <- grep("DF.classifications", colnames(wt@meta.data), value = TRUE)

DimPlot(wt, group.by = df_col)

singlet_cells <- rownames(wt@meta.data)[wt@meta.data[[df_col]] == "Singlet"]

wt_clean <- subset(wt, cells = singlet_cells)

ncol(wt)
ncol(wt_clean)

cat("Before doublet removal:", ncol(wt), "cells\n")
cat("After doublet removal:", ncol(wt_clean), "cells\n")
cat("Removed:", ncol(wt) - ncol(wt_clean), "doublets\n")
