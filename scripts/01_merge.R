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

# Merge cleaned objects
combined <- merge(
  wt_clean,
  y = c(i1d_clean, i7d_clean, i30d_clean),
  add.cell.ids = c("WT", "I1D", "I7D", "I30D"),
  project = "TendonHealing"
)

# Check merged object
combined
table(combined$condition)

# Save merged object
saveRDS(combined, "../data/combined_singlets_raw.rds")

combined <- readRDS("../data/combined_singlets_raw.rds")
combined
dim(combined)

ncol(combined)   # cells
nrow(combined)   # genes

Assays(combined)
DefaultAssay(combined)

object.size(combined)
format(object.size(combined), units = "GB")

# object has old data layers, pre-scaling, following code will clean this up 
# Join count layers into one counts layer
combined <- JoinLayers(combined)

Layers(combined[["RNA"]])

combined@reductions <- list()
combined@graphs <- list()

combined_light <- DietSeurat(
  combined,
  assays = "RNA",
  layers = "counts",
  dimreducs = NULL,
  graphs = NULL
)

format(object.size(combined_light), units = "GB")
Layers(combined_light[["RNA"]])

saveRDS(combined_light, "../data/combined_singlets_counts_only.rds")

rm(combined)
rm(combined_light)

gc()

