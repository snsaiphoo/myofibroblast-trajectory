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

# Check that scores were added
grep("_UCell$", colnames(fib@meta.data), value = TRUE)

# Save scored object
saveRDS(fib, "../../data/human_fib_ucell.rds")
