library(Seurat)
library(dplyr)
library(patchwork)

set.seed(123)

# Read cleaned singlet objects
wt_clean <- readRDS("../data/wt_singlets.rds")
i1d_clean <- readRDS("../data/i1d_singlets.rds")
i7d_clean <- readRDS("../data/i7d_singlets.rds")
i30d_clean <- readRDS("../data/i30d_singlets.rds")

# Make sure condition metadata is present
wt_clean$condition <- "WT"
i1d_clean$condition <- "I1D"
i7d_clean$condition <- "I7D"
i30d_clean$condition <- "I30D"

# Reduce object size before merging
wt_clean <- DietSeurat(wt_clean, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
i1d_clean <- DietSeurat(i1d_clean, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
i7d_clean <- DietSeurat(i7d_clean, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
i30d_clean <- DietSeurat(i30d_clean, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)

# Merge cleaned count-only objects
combined <- merge(
  wt_clean,
  y = c(i1d_clean, i7d_clean, i30d_clean),
  add.cell.ids = c("WT", "I1D", "I7D", "I30D"),
  project = "TendonHealing"
)

# Join layers after merge
combined <- JoinLayers(combined)

# Keep only counts in merged object
combined <- DietSeurat(
  combined,
  assays = "RNA",
  layers = "counts",
  dimreducs = NULL,
  graphs = NULL
)

# Check merged object
combined
table(combined$condition)
Layers(combined[["RNA"]])
format(object.size(combined), units = "GB")

# Save merged count-only object
saveRDS(combined, "../data/combined_singlets_counts_only.rds")

# Clean memory
rm(wt_clean, i1d_clean, i7d_clean, i30d_clean)
gc()
