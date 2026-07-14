# clean monocle3

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)

mesenchymal <- readRDS("../../data/mesenchymal_refined_annotated.rds")

p <- DimPlot(mesenchymal, group.by = "cell_type_refined", label = TRUE)

ggsave(
  filename = "../results/03_monocle3_v2/monocle3_step3_refined_celltypes.png",
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

cds <- as.cell_data_set(mesenchymal)

colnames(colData(cds))

table(colData(cds)$cell_type_refined)

cds <- preprocess_cds(cds, num_dim = 30)

cds <- reduce_dimension(cds, reduction_method = "UMAP")

cds <- cluster_cells(cds, reduction_method = "UMAP")

table(partitions(cds))

#    1     2     3     4     5
#10102  6035  3897    76    38

cds <- learn_graph(cds, use_partition = TRUE)

p1 <- plot_cells(
  cds,
  color_cells_by = "cell_type_refined",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3
)

p2 <- plot_cells(
  cds,
  color_cells_by = "condition",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3
)

ggsave(
  "../results/03_monocle3_v2/monocle3_trajectory_cell_type_refined.png",
  p1,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "../results/03_monocle3_v2/monocle3_trajectory_condition_refined.png",
  p2,
  width = 10,
  height = 8,
  dpi = 300
)

# partition 1
cds_p1 <- cds[, partitions(cds) == 1]

table(colData(cds_p1)$cell_type_refined)

table(
  colData(cds)$cell_type_refined,
  colData(cds)$condition
)

#                                     I1D I30D  I7D   WT
#  Activated tenocytes               1067   19   35   51
#  ECM-remodelling tenocytes            0 2968    7   13
#  Fibrochondrocyte-like tenocytes      0 2569   84    6
#  Mature tenocytes                    16  190    3 1192
#  Proinflammatory mesenchymal cells 1278    6    0    0
#  Proliferating mesenchymal cells     11   31  748    0
#  Repair fibroblasts                   0   59 3139    0
#  Repair-activated stromal cells       1   33 2089    0
#  Signaling stromal cells              0 1473   25  264
#  Stromal progenitor-like cells        9  774   30 1958

root_cells_p1 <- colnames(cds_p1)[
  colData(cds_p1)$cell_type_refined == "Stromal progenitor-like cells" &
    colData(cds_p1)$condition == "WT"
]

cds_p1 <- order_cells(cds_p1, root_cells = root_cells_p1)

library(ggplot2)
library(viridis)

# Pseudotime trajectory
p_pseudotime <- plot_cells(
  cds_p1,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_principal_points = FALSE
) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Pseudotime"
  ) +
  coord_cartesian(
    xlim = c(-1.0, 6.5),
    ylim = c(-9.5, 3)
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Make trajectory graph thicker
for (i in seq_along(p_pseudotime$layers)) {
  if (inherits(p_pseudotime$layers[[i]]$geom, "GeomPath")) {
    p_pseudotime$layers[[i]]$aes_params$colour <- "black"
    p_pseudotime$layers[[i]]$aes_params$linewidth <- 1.5
  }
}

# Save figure
ggsave(
  filename = "../results/03_monocle3_v2/monocle3_p1_pseudotime_final.png",
  plot = p_pseudotime,
  width = 7,
  height = 5,
  dpi = 300
)

# partition 2
cds_p2 <- cds[, partitions(cds) == 2]

table(colData(cds_p2)$cell_type_refined)

table(
  colData(cds_p2)$cell_type_refined,
  colData(cds_p2)$condition
)

root_cells_p2 <- colnames(cds_p2)[
  colData(cds_p2)$cell_type_refined == "Repair fibroblasts"
]

cds_p2 <- order_cells(
  cds_p2,
  root_cells = root_cells_p2
)

# Pseudotime trajectory
p_pseudotime <- plot_cells(
  cds_p2,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_principal_points = FALSE
) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Pseudotime"
  ) +
  coord_cartesian(
    xlim = c(-8, 4),
    ylim = c(5, 11)
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Make trajectory graph thicker
for (i in seq_along(p_pseudotime$layers)) {
  if (inherits(p_pseudotime$layers[[i]]$geom, "GeomPath")) {
    p_pseudotime$layers[[i]]$aes_params$colour <- "red"
    p_pseudotime$layers[[i]]$aes_params$linewidth <- 1.5
  }
}

# Save figure
ggsave(
  filename = "../results/03_monocle3_v2/monocle3_p2_pseudotime_final.png",
  plot = p_pseudotime,
  width = 7,
  height = 5,
  dpi = 300
)

# partition 3

# partition 2
cds_p3 <- cds[, partitions(cds) == 3]

table(colData(cds_p3)$cell_type_refined)

table(
  colData(cds_p3)$cell_type_refined,
  colData(cds_p3)$condition
)

root_cells_p3 <- colnames(cds_p3)[
  colData(cds_p3)$cell_type_refined == "Mature tenocytes"
]

cds_p3 <- order_cells(
  cds_p3,
  root_cells = root_cells_p3
)

# Pseudotime trajectory
p_pseudotime <- plot_cells(
  cds_p3,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_principal_points = FALSE
) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Pseudotime"
  ) +
  coord_cartesian(
    xlim = c(-12, -1),
    ylim = c(-7.0, 2)
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Make trajectory graph thicker
for (i in seq_along(p_pseudotime$layers)) {
  if (inherits(p_pseudotime$layers[[i]]$geom, "GeomPath")) {
    p_pseudotime$layers[[i]]$aes_params$colour <- "red"
    p_pseudotime$layers[[i]]$aes_params$linewidth <- 1.5
  }
}

# Save figure
ggsave(
  filename = "../results/03_monocle3_v2/monocle3_p3_pseudotime_final.png",
  plot = p_pseudotime,
  width = 7,
  height = 5,
  dpi = 300
)
