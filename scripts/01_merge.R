library(Seurat)
library(dplyr)
library(patchwork)

set.seed(123)

# Read cleaned singlet objects
wt_clean <- readRDS("data/wt_singlets.rds")
i1d_clean <- readRDS("data/i1d_singlets.rds")
i7d_clean <- readRDS("data/i7d_singlets.rds")
i30d_clean <- readRDS("data/i30d_singlets.rds")

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
saveRDS(combined, "data/combined_singlets_raw.rds")
