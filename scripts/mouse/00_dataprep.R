source("../../functions/seuratprocess.R")

library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)

# # Extract data
# unzip("data/GSM8767956_WT_matrix.zip", exdir = "data/")
# unzip("data/GSM8767957_I1D_matrix.zip", exdir = "data/")
# unzip("data/GSM8767958_I7D_matrix.zip", exdir = "data/")
# unzip("data/GSM8767959_I30D_matrix.zip", exdir = "data/")
# 
# # Read in the data 
wt_counts <- Read10X(data.dir = "../../data/WT_matrix/")
i1d_counts <- Read10X(data.dir = "../../data/I1D_matrix")
i7d_counts <- Read10X(data.dir = "../../data/I7D_matrix")
i30d_counts <- Read10X(data.dir = "../../data/I30D_matrix")

# # Convert into Seurat objects
wt <- CreateSeuratObject(counts = wt_counts, min.cells = 3, min.features = 200, project = "WT")
i1d <- CreateSeuratObject(counts = i1d_counts, min.cells = 3, min.features = 200, project = "I1D")
i7d <- CreateSeuratObject(counts = i7d_counts, min.cells = 3, min.features = 200, project = "I7D")
i30d <- CreateSeuratObject(counts = i30d_counts, min.cells = 3, min.features = 200, project = "I30D")
# 
set.seed(123)

# Add metadata for tracking
wt$condition <- "WT"
i1d$condition <- "I1D"
i7d$condition <- "I7D"
i30d$condition <- "I30D"

saveRDS(wt, file = "../../data/wt.rds")
saveRDS(i1d, file = "../../data/i1d.rds")
saveRDS(i7d, file = "../../data/i7d.rds")
saveRDS(i30d, file = "../../data/i30d.rds")

wt <- readRDS("../../data/wt.rds")
i1d <- readRDS("../../data/i1d.rds")
i7d <- readRDS("../../data/i7d.rds")
i30d <- readRDS("../../data/i30d.rds")

# Preprocessing for WT 
wt <- plot_qc(wt)
wt <- filter_qc(wt, 500, 5000, 10)
wt <- preprocess_pca(wt)

# Doublet Finder for WT
plot_elbow(wt, ndims = 50)
pcs_to_use <- 1:25
bcmvn <- find_best_pk(wt, pcs = pcs_to_use)
wt <- run_doubletfinder(wt, pcs = 1:25, pK = 0.28)

df_col <- grep("DF.classifications", colnames(wt@meta.data), value = TRUE)

DimPlot(wt, group.by = df_col)

singlet_cells <- rownames(wt@meta.data)[wt@meta.data[[df_col]] == "Singlet"]

wt_clean <- subset(wt, cells = singlet_cells)

ncol(wt)
ncol(wt_clean)

cat("Before doublet removal:", ncol(wt), "cells\n")
cat("After doublet removal:", ncol(wt_clean), "cells\n")
cat("Removed:", ncol(wt) - ncol(wt_clean), "doublets\n")

saveRDS(wt_clean, "../../data/wt_singlets.rds")

# Preprocessing for I1D
i1d <- plot_qc(i1d)
i1d <- filter_qc(i1d, 1000, 6000, 8)
i1d <- preprocess_pca(i1d)

# Doublet Finder for I1D
plot_elbow(i1d, ndims = 50)
pcs_to_use <- 1:25
bcmvn <- find_best_pk(i1d, pcs = pcs_to_use)
i1d <- run_doubletfinder(i1d, pcs = 1:25, pK = 0.23)

df_col <- grep("DF.classifications", colnames(i1d@meta.data), value = TRUE)

DimPlot(i1d, group.by = df_col)

singlet_cells <- rownames(i1d@meta.data)[i1d@meta.data[[df_col]] == "Singlet"]

i1d_clean <- subset(i1d, cells = singlet_cells)

ncol(i1d)
ncol(i1d_clean)

cat("Before doublet removal:", ncol(i1d), "cells\n")
cat("After doublet removal:", ncol(i1d_clean), "cells\n")
cat("Removed:", ncol(i1d) - ncol(i1d_clean), "doublets\n")

saveRDS(i1d_clean, "../../data/i1d_singlets.rds")

# Preprocessing for I7D
i7d <- plot_qc(i7d)
i7d <- filter_qc(i7d, 1500, 6500, 8)
i7d <- preprocess_pca(i7d)

# Doublet Finder for I7D
plot_elbow(i7d, ndims = 50)
pcs_to_use <- 1:25
bcmvn <- find_best_pk(i7d, pcs = pcs_to_use)
i7d <- run_doubletfinder(i7d, pcs = 1:25, pK = 0.05)

df_col <- grep("DF.classifications", colnames(i7d@meta.data), value = TRUE)

DimPlot(i7d, group.by = df_col)

singlet_cells <- rownames(i7d@meta.data)[i7d@meta.data[[df_col]] == "Singlet"]

i7d_clean <- subset(i7d, cells = singlet_cells)

ncol(i7d)
ncol(i7d_clean)

cat("Before doublet removal:", ncol(i7d), "cells\n")
cat("After doublet removal:", ncol(i7d_clean), "cells\n")
cat("Removed:", ncol(i7d) - ncol(i7d_clean), "doublets\n")

saveRDS(i7d_clean, "../../data/i7d_singlets.rds")

# Preprocessing for I30D
i30d <- plot_qc(i30d)
i30d <- filter_qc(i30d, 1000, 6000, 10)
i30d <- preprocess_pca(i30d)

# Doublet Finder for I3D
plot_elbow(i30d, ndims = 50)
pcs_to_use <- 1:25
bcmvn <- find_best_pk(i30d, pcs = pcs_to_use)
i30d <- run_doubletfinder(i30d, pcs = 1:25, pK = 0.21)

df_col <- grep("DF.classifications", colnames(i30d@meta.data), value = TRUE)

DimPlot(i30d, group.by = df_col)

singlet_cells <- rownames(i30d@meta.data)[i30d@meta.data[[df_col]] == "Singlet"]

i30d_clean <- subset(i30d, cells = singlet_cells)

ncol(i30d)
ncol(i30d_clean)

cat("Before doublet removal:", ncol(i30d), "cells\n")
cat("After doublet removal:", ncol(i30d_clean), "cells\n")
cat("Removed:", ncol(i30d) - ncol(i30d_clean), "doublets\n")

saveRDS(i30d_clean, "../../data/i30d_singlets.rds")

