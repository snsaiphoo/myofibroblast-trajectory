library(Seurat)
library(dplyr)
library(patchwork)

# Extract data
unzip("data/GSM8767956_WT_matrix.zip", exdir = "data/")
unzip("data/GSM8767957_I1D_matrix.zip", exdir = "data/")
unzip("data/GSM8767958_I7D_matrix.zip", exdir = "data/")
unzip("data/GSM8767959_I30D_matrix.zip", exdir = "data/")

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

# Add metadata for tracking
wt$condition <- "WT"
i1d$condition <- "I1D"
i7d$condition <- "I7D"
i30d$condition <- "I30D"

#


