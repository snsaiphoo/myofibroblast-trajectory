library(Seurat)
library(UCell)

# Load fibroblast Seurat object
fib <- readRDS("../../data/human_fib_sub.rds")

# Load Ensembl gene sets
gene_sets <- readRDS("../../data/all_gene_sets_human_ensembl.rds")

# Check gene sets
names(gene_sets)
sapply(gene_sets, length)

# Run UCell scoring
fib <- AddModuleScore_UCell(
  fib,
  features = gene_sets,
)

# ============================================================
# Donor-level UCell scores
# Each donor contributes equally to the analysis
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)

# Check UCell score columns
score_cols <- grep(
  "_UCell$",
  colnames(fib@meta.data),
  value = TRUE
)

score_cols

# ============================================================
# Calculate mean UCell score per donor and cell type
# ============================================================

fib_donor_scores <- fib@meta.data %>%
  filter(
    !is.na(donor_id),
    !is.na(condition),
    !is.na(author_cell_type)
  ) %>%
  group_by(
    donor_id,
    condition,
    author_cell_type
  ) %>%
  summarise(
    across(
      all_of(score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  )

head(fib_donor_scores)

# ============================================================
# Condition-level summary
# Each donor has equal weight
# ============================================================

fib_condition_scores <- fib_donor_scores %>%
  group_by(
    condition,
    author_cell_type
  ) %>%
  summarise(
    across(
      all_of(score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

head(fib_condition_scores)

# ============================================================
# Donor-level mean and SD
# ============================================================

fib_condition_summary <- fib_donor_scores %>%
  group_by(
    condition,
    author_cell_type
  ) %>%
  summarise(
    across(
      all_of(score_cols),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n_donors = n(),
    .groups = "drop"
  )

head(fib_condition_summary)

# Check that scores were added
grep("_UCell$", colnames(fib@meta.data), value = TRUE)

# Save scored object
saveRDS(fib, "../../data/human_fib_ucell.rds")
