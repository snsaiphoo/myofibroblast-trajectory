# slingshot code
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)

seu <- readRDS("../../data/mesenchymal_geneset_scored_refined.rds")

sce <- as.SingleCellExperiment(seu)

# Check conversion
sce

# Confirm metadata transferred
colnames(colData(sce))

# Confirm cell type labels transferred
table(sce$cell_type_refined)

# Confirm condition transferred
table(sce$condition)

sce <- slingshot(
  sce,
  clusterLabels = "cell_type_refined",
  reducedDim = "UMAP",
  start.clus = "Mature tenocytes"
)

# Check results
sce
SlingshotDataSet(sce)

slingLineages(sce)


# For Slingshot, the starting cluster was specified as Mature tenocytes. Although Slingshot initializes trajectories at the cluster level rather than individual cells, this cluster was composed predominantly of WT cells (85.1%), making it closely comparable to the WT Mature tenocyte root used for Monocle3.


# Basic Slingshot trajectory plot on UMAP
png(
  filename = "../results/05_slingshot/Figure1_Slingshot_Trajectory_CellTypes.png",
  width = 3000,
  height = 2600,
  res = 300
)

umap <- reducedDims(sce)$UMAP

cell_cols <- as.factor(sce$cell_type_refined)

plot(
  umap,
  col = as.numeric(cell_cols),
  pch = 16,
  asp = 1,
  cex = 0.35,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "Slingshot trajectories by refined cell type"
)

lines(
  SlingshotDataSet(sce),
  lwd = 3,
  col = "black"
)

legend(
  "topright",
  legend = levels(cell_cols),
  col = seq_along(levels(cell_cols)),
  pch = 16,
  cex = 0.8,
  bty = "n"
)

dev.off()


# CLEANED FIGURE 1 per lineages

## UMAP coordinates
umap <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap) <- c("UMAP_1", "UMAP_2")

plot_df <- cbind(
  umap,
  cell_type_refined = sce$cell_type_refined,
  condition = sce$condition
)

## Lineage membership weights
weights <- slingCurveWeights(sce)

## Curves
curves <- slingCurves(sce)

for(i in seq_along(curves)){
  
  df <- plot_df
  
  ## Highlight cells assigned to this lineage
  df$highlight <- weights[,i] > 0.05
  
  curve_df <- data.frame(
    UMAP_1 = curves[[i]]$s[,1],
    UMAP_2 = curves[[i]]$s[,2]
  )
  
  p <- ggplot() +
    
    ## Background cells
    geom_point(
      data = df[!df$highlight,],
      aes(UMAP_1, UMAP_2),
      color = "grey85",
      size = 0.22
    ) +
    
    ## Lineage cells
    geom_point(
      data = df[df$highlight,],
      aes(UMAP_1, UMAP_2,
          color = cell_type_refined),
      size = 0.28,
      alpha = 0.9
    ) +
    
    ## Slingshot curve
    geom_path(
      data = curve_df,
      aes(UMAP_1, UMAP_2),
      linewidth = 1.2,
      color = "black"
    ) +
    
    coord_equal() +
    theme_classic(base_size = 15) +
    
    theme(
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 14)
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(size = 2, alpha = 1)
      )
    ) +
    
    labs(
      title = paste("Slingshot Lineage", i),
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Cell type"
    )
  
  ggsave(
    paste0("../results/05_slingshot/Figure1", LETTERS[i],
           "_LineageHighlight.png"),
    p,
    width = 8,
    height = 6.5,
    dpi = 300
  )
}

# Slingshot Pseudotime

## UMAP coordinates
umap <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap) <- c("UMAP_1", "UMAP_2")

plot_df <- cbind(
  umap,
  cell_type_refined = sce$cell_type_refined,
  condition = sce$condition
)

## Slingshot pseudotime and weights
pt <- slingPseudotime(sce)
weights <- slingCurveWeights(sce)
curves <- slingCurves(sce)

for(i in seq_along(curves)){
  
  df <- plot_df
  df$pseudotime <- pt[, i]
  df$highlight <- weights[, i] > 0.05 & !is.na(df$pseudotime)
  
  curve_df <- data.frame(
    UMAP_1 = curves[[i]]$s[, 1],
    UMAP_2 = curves[[i]]$s[, 2]
  )
  
  p <- ggplot() +
    
    geom_point(
      data = df[!df$highlight, ],
      aes(UMAP_1, UMAP_2),
      color = "grey85",
      size = 0.22
    ) +
    
    geom_point(
      data = df[df$highlight, ],
      aes(UMAP_1, UMAP_2, color = pseudotime),
      size = 0.35,
      alpha = 0.9
    ) +
    
    geom_path(
      data = curve_df,
      aes(UMAP_1, UMAP_2),
      color = "black",
      linewidth = 1.1
    ) +
    
    coord_equal() +
    scale_color_viridis_c(option = "plasma", na.value = "grey85") +
    theme_classic(base_size = 15) +
    
    theme(
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 14)
    ) +
    
    labs(
      title = paste("Slingshot Lineage", i, "pseudotime"),
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Pseudotime"
    )
  
  ggsave(
    paste0("../results/05_slingshot/Figure2", LETTERS[i],
           "_Slingshot_Pseudotime.png"),
    p,
    width = 8,
    height = 6.5,
    dpi = 300
  )
}

#Despite differences in trajectory construction, Slingshot recovered three biologically meaningful lineages originating from mature tenocytes and showed consistent progression toward repair-associated and inflammatory cell states, supporting the robustness of the Monocle3-derived trajectory.

# UCELL FIGURES

pt <- slingPseudotime(sce)
weights <- slingCurveWeights(sce)

meta_df <- seu@meta.data

for(i in 1:ncol(pt)) {
  
  lineage_name <- paste0("Lineage", i)
  
  df <- meta_df %>%
    mutate(
      pseudotime = pt[, i],
      weight = weights[, i]
    ) %>%
    filter(!is.na(pseudotime), weight > 0.05) %>%
    select(cell_type_refined, condition, pseudotime, all_of(ucell_cols)) %>%
    pivot_longer(
      cols = all_of(ucell_cols),
      names_to = "program",
      values_to = "UCell_score"
    ) %>%
    mutate(
      program = recode(program, !!!program_labels)
    )
  
  p <- ggplot(df, aes(x = pseudotime, y = UCell_score)) +
    geom_point(alpha = 0.18, size = 0.25) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
    facet_wrap(~ program, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 14) +
    theme(
      strip.text = element_text(size = 13, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 11)
    ) +
    labs(
      title = paste0("UCell programs along Slingshot ", lineage_name),
      x = "Slingshot pseudotime",
      y = "UCell score"
    )
  
  ggsave(
    filename = paste0("../results/05_slingshot/Figure3_", lineage_name,
                      "_UCell_vs_Pseudotime.png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Slingshot validated the Monocle3-derived repair landscape by recovering mature tenocytes as the common origin and resolving the main repair-associated partition into distinct remodeling, stromal/progenitor, and inflammatory lineages. Mechanotransduction-associated programs were most consistently enriched along the remodeling lineage, supporting a state-specific rather than globally linear activation model.

# ordering the pseudotime 

pt <- slingPseudotime(sce)
weights <- slingCurveWeights(sce)

lineage_summary <- data.frame(
  cell_type_refined = sce$cell_type_refined,
  condition = sce$condition,
  Lineage1_pseudotime = pt[, "Lineage1"],
  Lineage2_pseudotime = pt[, "Lineage2"],
  Lineage3_pseudotime = pt[, "Lineage3"],
  Lineage1_weight = weights[, "Lineage1"],
  Lineage2_weight = weights[, "Lineage2"],
  Lineage3_weight = weights[, "Lineage3"]
)

## Median pseudotime by cell type for each lineage
for(i in 1:3){
  
  pt_col <- paste0("Lineage", i, "_pseudotime")
  wt_col <- paste0("Lineage", i, "_weight")
  
  summary_i <- lineage_summary %>%
    filter(!is.na(.data[[pt_col]]), .data[[wt_col]] > 0.05) %>%
    group_by(cell_type_refined) %>%
    summarise(
      n_cells = n(),
      median_pseudotime = median(.data[[pt_col]], na.rm = TRUE),
      mean_pseudotime = mean(.data[[pt_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_pseudotime)
  
  print(paste("Lineage", i))
  print(summary_i)
  
  write.csv(
    summary_i,
    paste0("../results/05_slingshot/Lineage", i, "_median_pseudotime_by_celltype.csv"),
    row.names = FALSE
  )
}

# Pseudotime Heatmaps

pt <- slingPseudotime(sce)
weights <- slingCurveWeights(sce)
meta_df <- seu@meta.data

ucell_cols <- c(
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

program_labels <- c(
  YAP_TEAD_UCell = "YAP/TEAD",
  Integrin_FAK_UCell = "Integrin/FAK",
  RhoA_ROCK_UCell = "RhoA/ROCK",
  Myofibroblast_UCell = "Myofibroblast",
  GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell = "ECM organization",
  GOBP_FIBROBLAST_ACTIVATION_UCell = "Fibroblast activation",
  HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammation",
  HALLMARK_TGF_BETA_SIGNALING_UCell = "TGF-beta",
  REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell = "TGF-beta/SMAD"
)

make_lineage_heatmap <- function(lineage_number, n_bins = 50) {
  
  lineage_name <- paste0("Lineage", lineage_number)
  
  df <- meta_df %>%
    mutate(
      pseudotime = pt[, lineage_name],
      weight = weights[, lineage_name]
    ) %>%
    filter(!is.na(pseudotime), weight > 0.05) %>%
    select(pseudotime, all_of(ucell_cols)) %>%
    mutate(
      pseudotime_bin = ntile(pseudotime, n_bins)
    ) %>%
    group_by(pseudotime_bin) %>%
    summarise(
      across(all_of(ucell_cols), mean, na.rm = TRUE),
      mean_pseudotime = mean(pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = all_of(ucell_cols),
      names_to = "program",
      values_to = "mean_UCell"
    ) %>%
    group_by(program) %>%
    mutate(
      scaled_score = as.numeric(scale(mean_UCell))
    ) %>%
    ungroup() %>%
    mutate(
      program = recode(program, !!!program_labels),
      program = factor(
        program,
        levels = c(
          "YAP/TEAD",
          "Integrin/FAK",
          "RhoA/ROCK",
          "TGF-beta",
          "TGF-beta/SMAD",
          "Myofibroblast",
          "ECM organization",
          "Fibroblast activation",
          "Inflammation"
        )
      )
    )
  
  p <- ggplot(df, aes(x = pseudotime_bin, y = program, fill = scaled_score)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "navy",
      mid = "white",
      high = "firebrick",
      midpoint = 0,
      name = "Scaled\nUCell"
    ) +
    theme_classic(base_size = 15) +
    theme(
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 11),
      axis.title = element_text(size = 15),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    ) +
    labs(
      title = paste0("UCell program dynamics along Slingshot ", lineage_name),
      x = "Slingshot pseudotime bin",
      y = "UCell program"
    )
  
  ggsave(
    filename = paste0("../results/05_slingshot/Figure4_", lineage_name, "_UCell_Pseudotime_Heatmap.png"),
    plot = p,
    width = 9,
    height = 5.5,
    dpi = 300
  )
  
  write.csv(
    df,
    file = paste0("../results/05_slingshot/Figure4_", lineage_name, "_UCell_Pseudotime_Heatmap_Data.csv"),
    row.names = FALSE
  )
}

make_lineage_heatmap(1)
make_lineage_heatmap(2)
make_lineage_heatmap(3)

# Pseudotime UMAP UCELL

umap <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap) <- c("UMAP_1", "UMAP_2")

plot_df <- cbind(
  umap,
  seu@meta.data
)

curves <- slingCurves(sce)

curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  data.frame(
    UMAP_1 = curves[[i]]$s[, 1],
    UMAP_2 = curves[[i]]$s[, 2],
    lineage = paste0("Lineage ", i)
  )
}))

ucell_cols <- c(
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

program_labels <- c(
  YAP_TEAD_UCell = "YAP/TEAD",
  Integrin_FAK_UCell = "Integrin/FAK",
  RhoA_ROCK_UCell = "RhoA/ROCK",
  Myofibroblast_UCell = "Myofibroblast",
  GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell = "ECM organization",
  GOBP_FIBROBLAST_ACTIVATION_UCell = "Fibroblast activation",
  HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammation",
  HALLMARK_TGF_BETA_SIGNALING_UCell = "TGF-beta",
  REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell = "TGF-beta/SMAD"
)

for(score in ucell_cols){
  
  label <- program_labels[[score]]
  file_label <- gsub("[^A-Za-z0-9]+", "_", label)
  
  p <- ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
    geom_point(
      aes(color = .data[[score]]),
      size = 0.35,
      alpha = 0.9
    ) +
    geom_path(
      data = curve_df,
      aes(x = UMAP_1, y = UMAP_2, group = lineage),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.9
    ) +
    coord_equal() +
    scale_color_viridis_c(option = "plasma", name = paste0(label, "\nUCell")) +
    theme_classic(base_size = 15) +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 20)
    ) +
    labs(
      title = paste0(label, " activity along Slingshot trajectories"),
      x = "UMAP 1",
      y = "UMAP 2"
    )
  
  ggsave(
    filename = paste0("../results/05_slingshot/Figure5_", file_label, "_UCell_Slingshot_UMAP.png"),
    plot = p,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )
}