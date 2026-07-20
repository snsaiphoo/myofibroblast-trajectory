library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

comparison <- readRDS(
  "../../data/CrossSpecies_PathwayMatrix.rds"
)

comparison <- comparison[, c(
  "Mouse_Repair_Activated_Stromal",
  "Mouse_Repair_Fibroblasts",
  "ADAM12hi fibroblasts",
  "FBLN1hi fibroblasts",
  "NR4A1hi fibroblasts",
  "ABCA10hi fibroblasts"
)]

colnames(comparison) <- c(
  "Activated\nStromal",
  "Repair",
  "ADAM12hi",
  "FBLN1hi",
  "NR4A1hi",
  "ABCA10hi"
)

comparison_scaled <- t(
  scale(
    t(as.matrix(comparison))
  )
)

species <- c(
  "Mouse",
  "Mouse",
  "Human",
  "Human",
  "Human",
  "Human"
)

population <- c(
  "Repair",
  "Repair",
  "Repair",
  "Homeostatic",
  "Homeostatic",
  "Homeostatic"
)

ha <- HeatmapAnnotation(
  Species = species,
  Population = population,
  col = list(
    Species = c(
      Mouse = "#4F81BD",
      Human = "#E67E22"
    ),
    Population = c(
      Repair = "#C0392B",
      Homeostatic = "#27AE60"
    )
  ),
  annotation_name_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  )
)


ht <- Heatmap(
  comparison_scaled,
  name = "Z-score",
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_names_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  row_names_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  row_names_side = "left",
  column_title =
    "Cross-species conservation of mechanotransduction pathways",
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c(
      "#3B4CC0",
      "white",
      "#B40426"
    )
  )
)


png(
  "../../results/07_snRNA/Figure6_CrossSpeciesHeatmap.png",
  width = 2500,
  height = 1800,
  res = 300
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

####

# comparison <- readRDS(
#   "../../data/CrossSpecies_PathwayMatrix.rds"
# )

# colnames(comparison) <- c(
#   "Repair",
#   "Activated Stromal",
#   "ADAM12hi",
#   "NR4A1hi",
#   "ABCA10hi",
#   "FBLN1hi"
# )

mat <- as.matrix(comparison)

library(pheatmap)

annotation <- data.frame(
  
  Species = c(
    "Mouse",
    "Mouse",
    "Human",
    "Human",
    "Human",
    "Human"
  ),
  
  Population = c(
    "Repair",
    "Repair",
    "Repair",
    "Homeostatic",
    "Homeostatic",
    "Homeostatic"
  )
  
)

rownames(annotation) <- colnames(mat)


ann_colors <- list(
  Species = c(
    Mouse = "#4F81BD",
    Human = "#E67E22"
  ),
  Population = c(
    Repair = "#C0392B",
    Intermediate = "#E67E22",
    Homeostatic = "#27AE60"
  )
)

png(
  "../../results/07_snRNA/FigureS1_UnscaledHeatmap.png",
  width = 3200,
  height = 1800,
  res = 300
)

pheatmap(
  mat,
  name = "Average\nUCell score",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(c(
    "white",
    "#F7DADA",
    "#E88989",
    "#D14B4B",
    "#B40426"
  ))(100),
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 16,
  border_color = "grey90",
  main = "Average mechanotransduction pathway activity"
)

dev.off()



