# monocle3 ucell overlay
# continued environment from 08_monocle3_2.R
# output directories already exist from previous script

mesenchymal_scored <- readRDS("../data/mesenchymal_geneset_scored_refined.rds")

# transferring scores to monocle3 cds

# UCell columns
ucell_cols <- grep(
  "UCell",
  colnames(mesenchymal_scored@meta.data),
  value = TRUE
)

# Add them to the Monocle object
colData(cds)[, ucell_cols] <-
  mesenchymal_scored@meta.data[
    colnames(cds),
    ucell_cols
  ]

# recreated after adding the metadata

cds_p1 <- cds[, partitions(cds) == 1]
cds_p2 <- cds[, partitions(cds) == 2]
cds_p3 <- cds[, partitions(cds) == 3]

cds_p1 <- order_cells(cds_p1, root_cells = root_cells_p1)
cds_p2 <- order_cells(cds_p2, root_cells = root_cells_p2)
cds_p3 <- order_cells(cds_p3, root_cells = root_cells_p3)

# all 9 gene sets 
gene_sets <- c(
  "YAP_TEAD_UCell",
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "Myofibroblast_UCell",
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell"
)

legend_names <- c(
  YAP_TEAD_UCell = "YAP/TEAD",
  Integrin_FAK_UCell = "Integrin-FAK",
  RhoA_ROCK_UCell = "RhoA-ROCK",
  Myofibroblast_UCell = "Myofibroblast",
  GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell = "ECM Organization",
  GOBP_FIBROBLAST_ACTIVATION_UCell = "Fibroblast Activation",
  HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammatory Response",
  HALLMARK_TGF_BETA_SIGNALING_UCell = "TGF-β Signaling",
  REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell = "SMAD Signaling"
)

# reusable plotting function
library(ggplot2)
library(viridis)

plot_ucell_overlay <- function(cds_obj,
                               score,
                               xlim = NULL,
                               ylim = NULL){

  p <- plot_cells(
    cds_obj,
    color_cells_by = score,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE,
    label_principal_points = FALSE
  ) +
    scale_color_viridis_c(
      option = "magma",
      name = legend_names[[score]]
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )

  if (!is.null(xlim) & !is.null(ylim)) {
    p <- p +
      coord_cartesian(
        xlim = xlim,
        ylim = ylim
      )
  }

  ## Darker trajectory graph
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPath")) {
      p$layers[[i]]$aes_params$colour <- "black"
      p$layers[[i]]$aes_params$linewidth <- 1.5
    }
  }

  return(p)
}

# partition 1
for(score in gene_sets){

  p <- plot_ucell_overlay(
    cds_p1,
    score,
    xlim = c(-1, 6.5),
    ylim = c(-9.5, 3)
  )

  ggsave(
    filename = paste0(
      "../results/monocle3_v2/partition1/",
      score,
      ".png"
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# partition 2
for(score in gene_sets){

  p <- plot_ucell_overlay(cds_p2, score, xlim = c(-8, 4), ylim = c(5, 11))

  ggsave(
    filename = paste0(
      "../results/monocle3_v2/partition2/",
      score,
      ".png"
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# partition 3 
for(score in gene_sets){

  p <- plot_ucell_overlay(cds_p3, score, xlim = c(-12, -1), ylim = c(-7.0, 2))

  ggsave(
    filename = paste0(
      "../results/monocle3_v2/partition3/",
      score,
      ".png"
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}
