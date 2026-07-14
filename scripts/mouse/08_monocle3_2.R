# Run 08_monocle3.R first to establish environment for script

root_partition <- function(cds_part, outdir, root_condition = NULL, root_celltype = NULL, root_cells = NULL) {
  if (is.null(root_cells)) {
    keep <- rep(TRUE, ncol(cds_part))

    if (!is.null(root_condition)) {
      keep <- keep & colData(cds_part)$condition == root_condition
    }

    if (!is.null(root_celltype)) {
      keep <- keep & colData(cds_part)$cell_type_refined == root_celltype
    }

    root_cells <- colnames(cds_part)[keep]
  }

  if (length(root_cells) == 0) {
    stop("No root cells found. Check the condition/cell type labels.")
  }

  cat("Number of root cells:", length(root_cells), "\n")

  cds_part <- order_cells(
    cds_part,
    root_cells = root_cells
  )

  saveRDS(
    cds_part,
    file.path(outdir, "Partition_rooted.rds")
  )

  # Cell Types

  p <- plot_cells(
    cds_part,
    color_cells_by = "cell_type_refined",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 4,
    cell_size = 0.8
  ) +
    scale_color_manual(values = celltype_cols) +
    theme_monocle()

  print(p)

  ggsave(
    file.path(outdir, "CellTypes.png"),
    p,
    width = 9,
    height = 8,
    dpi = 600,
    bg = "white"
  )

  # Timepoints

  p <- plot_cells(
    cds_part,
    color_cells_by = "condition",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 4,
    cell_size = 0.8
  ) +
    scale_color_manual(values = condition_cols) +
    theme_monocle()

  print(p)

  ggsave(
    file.path(outdir, "Timepoints.png"),
    p,
    width = 9,
    height = 8,
    dpi = 600,
    bg = "white"
  )

  # Pseudotime

  p <- plot_cells(
    cds_part,
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 4,
    cell_size = 0.8
  ) +
    scale_color_viridis_c(option = "magma") +
    theme_monocle()

  print(p)

  ggsave(
    file.path(outdir, "Pseudotime.png"),
    p,
    width = 9,
    height = 8,
    dpi = 600,
    bg = "white"
  )

  return(cds_part)
}

# Partition 1 (Homeostasis / Remodeling)

cds_p1 <- root_partition(
  cds_part = cds_p1,
  outdir = "../results/03_monocle3_v2/partition1",
  root_condition = "WT",
  root_celltype = "Stromal progenitor-like cells"
)

# Partition 2 (Peak Repair)

cds_p2 <- root_partition(
  cds_part = cds_p2,
  outdir = "../results/03_monocle3_v2/partition2",
  root_celltype = "Repair fibroblasts"
)

# Partition 3 (Acute Injury)

cds_p3 <- root_partition(
  cds_part = cds_p3,
  outdir = "../results/03_monocle3_v2/partition3",
  root_condition = "WT",
  root_celltype = "Mature tenocytes"
)
