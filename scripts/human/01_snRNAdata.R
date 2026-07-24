library(Seurat)
library(reticulate)
library(anndata)
library(dplyr)
library(ggplot2)
library(patchwork)

# data <- read_h5ad("../../data/human_rupture.h5ad")
# 
# data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 200, min.cells = 3)
# 
# saveRDS(data,"../../data/human_rupture.rds")

human <- readRDS("../../data/human_rupture.rds")

# Data Inspection
table(
  human$author_cell_type,
  human$disease_status
)

human

dim(human)

Reductions(human)
Assays(human)
table(human$disease_status)
table(human$author_cell_type)

# Create the condition variable 
human$condition <- ifelse(
  human$disease_status == "Quadriceps rupture",
  "Rupture",
  "Healthy"
)

table(human$condition)

# Data Quality Check
VlnPlot(
  human,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "subsets_mito_percent"
  ),
  group.by = "condition",
  pt.size = 0
)

set.seed(123)
# preprocess the data
human <- NormalizeData(human)

human <- FindVariableFeatures(
  human,
  selection.method = "vst",
  nfeatures = 2000
)

human <- ScaleData(human)

human <- RunPCA(human)

ElbowPlot(human, ndims=50)

human <- RunUMAP(
  human,
  dims = 1:20
)

saveRDS(human,"../../data/human_rupture_processed.rds")

human <- readRDS("../../data/human_rupture_processed.rds")

DimPlot(
  human,
  reduction = "umap",
  group.by = "disease"
)

# subset the fibroblasts
fib <- subset(
  human,
  subset = author_cell_type %in% c(
    "ABCA10hi fibroblasts",
    "ADAM12hi fibroblasts",
    "FBLN1hi fibroblasts",
    "NR4A1hi fibroblasts"
  )
)

saveRDS(fib,"../../data/human_fib_sub.rds")

DimPlot(
  fib,
  group.by = "author_cell_type",
  label = TRUE,
  repel = TRUE
)

prop.table(
  table(
    fib$author_cell_type,
    fib$condition
  ),
  margin = 2
)

library(dplyr)
library(ggplot2)

# Fibroblast metadata
fib_meta <- human@meta.data %>%
  filter(author_cell_type %in% c(
    "ABCA10hi fibroblasts",
    "ADAM12hi fibroblasts",
    "FBLN1hi fibroblasts",
    "NR4A1hi fibroblasts"
  ))

# Count cells
fib_counts <- fib_meta %>%
  count(condition, author_cell_type) %>%
  group_by(condition) %>%
  mutate(percent = n / sum(n) * 100)

fib_counts

p <- ggplot(fib_counts,
       aes(x = condition,
           y = percent,
           fill = author_cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    x = "",
    y = "Fibroblast composition (%)",
    fill = "Fibroblast subtype"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_brewer(palette = "Set2")

ggsave(
  filename = "../../results/07_snRNA/Human_Fibroblast_Composition.png",
  plot = p,
  width = 5.5,
  height = 4.5,
  units = "in",
  dpi = 600,
)

p <- p +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

ggsave(
  filename = "../../results/07_snRNA/Human_Fibroblast_Composition.png",
  plot = p,
  width = 5.5,
  height = 4.5,
  units = "in",
  dpi = 600,
)
