## Cell Cycle Analysis (Mouse)

library(Seurat)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)

sce <- readRDS("../data/sce_slingshot_cellcycle.rds")

outdir <- "../results/cellcycle"

# Convert Seurat's built-in human gene lists to mouse

s.genes.mouse <- str_to_title(tolower(cc.genes.updated.2019$s.genes))
g2m.genes.mouse <- str_to_title(tolower(cc.genes.updated.2019$g2m.genes))

# Keep only genes present in the dataset
s.genes.mouse <- intersect(s.genes.mouse, rownames(seu))
g2m.genes.mouse <- intersect(g2m.genes.mouse, rownames(seu))

cat("Number of S phase genes found:", length(s.genes.mouse), "\n")
cat("Number of G2M phase genes found:", length(g2m.genes.mouse), "\n")

# Score cell cycle

seu <- CellCycleScoring(
  object = seu,
  s.features = s.genes.mouse,
  g2m.features = g2m.genes.mouse,
  set.ident = FALSE
)

saveRDS(seu, "../data/mesenchymal_cellcycle_scored.rds")

# Figures

library(ggplot2)

p_s <- FeaturePlot(
  seu,
  features = "S.Score",
  reduction = "umap",
  order = TRUE,
  pt.size = 0.4,
  cols = c("grey95", "#FDB863", "#D73027")
) +
  ggtitle("S Phase Score") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


ggsave(
  filename = file.path(outdir, "Figure1A_UMAP_SScore.png"),
  plot = p_s,
  width = 7,
  height = 6,
  dpi = 300
)

p_g2m <- FeaturePlot(
  seu,
  features = "G2M.Score",
  reduction = "umap",
  order = TRUE,
  pt.size = 0.4,
  cols = c("grey95", "#74ADD1", "#313695")
) +
  ggtitle("G2/M Phase Score") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_g2m

ggsave(
  filename = file.path(outdir, "Figure1B_UMAP_G2MScore.png"),
  plot = p_g2m,
  width = 7,
  height = 6,
  dpi = 300
)

p_phase <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "Phase",
  pt.size = 0.4
) +
  ggtitle("Cell Cycle Phase") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_phase

ggsave(
  filename = file.path(outdir, "Figure1C_UMAP_Phase.png"),
  plot = p_phase,
  width = 7,
  height = 6,
  dpi = 300
)

# Figure 2

# marker_dir <- file.path("cellcycle", "marker_plots")
# 
# if (!dir.exists(marker_dir)) {
#   dir.create(marker_dir)
# }

markers <- c(
  "Mki67",
  "Top2a",
  "Pcna",
  "Cdk1",
  "Mcm5",
  "Tyms"
)

pt <- slingPseudotime(sce)

df <- data.frame(
  cell = colnames(sce),
  Lineage1 = pt[,1],
  Lineage2 = pt[,2],
  Lineage3 = pt[,3],
  S.Score = colData(sce)$S.Score,
  G2M.Score = colData(sce)$G2M.Score,
  Phase = colData(sce)$Phase,
  CellType = colData(sce)$cell_type_refined
)

df1 <- df %>%
  filter(!is.na(Lineage1))

p1 <- ggplot(df1,
             aes(Lineage1, S.Score)) +
  geom_point(
    aes(color = CellType),
    alpha = 0.45,
    size = 0.7
  ) +
  geom_smooth(
    color = "black",
    linewidth = 1.2,
    se = TRUE
  ) +
  labs(
    title = "Lineage 1: S Phase Score",
    x = "Slingshot pseudotime",
    y = "S.Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p1

ggsave(
  file.path(
    pseudotime_dir,
    "Lineage1_SScore.png"
  ),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

# Lineage 1
p2 <- ggplot(df1,
             aes(Lineage1, G2M.Score)) +
  geom_point(
    aes(color = CellType),
    alpha = 0.45,
    size = 0.7
  ) +
  geom_smooth(
    color = "black",
    linewidth = 1.2,
    se = TRUE
  ) +
  labs(
    title = "Lineage 1: G2/M Phase Score",
    x = "Slingshot pseudotime",
    y = "G2M.Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p2

ggsave(
  file.path(
    outdir,
    "Lineage1_G2MScore.png"
  ),
  p2,
  width = 8,
  height = 5,
  dpi = 300
)

# bin pseudotime

library(dplyr)

df_phase <- df %>%
  filter(!is.na(Lineage1)) %>%
  mutate(
    PT_bin = cut(
      Lineage1,
      breaks = 20,
      include.lowest = TRUE
    )
  )

phase_summary <- df_phase %>%
  group_by(PT_bin, Phase) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(PT_bin) %>%
  mutate(
    proportion = n / sum(n)
  )

phase_summary$PT_bin <- factor(
  phase_summary$PT_bin,
  levels = unique(df_phase$PT_bin)
)

library(ggplot2)

p_phase_pt <- ggplot(
  phase_summary,
  aes(
    x = PT_bin,
    y = proportion,
    fill = Phase
  )
) +
  geom_col(width = 0.95) +
  scale_fill_manual(
    values = c(
      G1 = "#F8766D",
      S = "#619CFF",
      G2M = "#00BA38"
    )
  ) +
  labs(
    title = "Lineage 1: Cell-cycle phase along pseudotime",
    x = "Slingshot pseudotime",
    y = "Cell proportion"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_phase_pt

ggsave(
  file.path(
    outdir,
    "Lineage1_CellCyclePhase.png"
  ),
  p_phase_pt,
  width = 9,
  height = 5,
  dpi = 300
)

# Dotplot
## Figure 8 - Canonical proliferation markers DotPlot

library(Seurat)
library(ggplot2)
library(viridis)

# Create output folder if needed
if (!dir.exists("cellcycle")) {
  dir.create("cellcycle")
}

# Proliferation markers
markers <- c(
  "Mki67",
  "Top2a",
  "Cdk1",
  "Pcna",
  "Mcm5",
  "Tyms"
)

# Keep only genes present in the dataset
markers_present <- markers[markers %in% rownames(seu)]

cat("Markers found:\n")
print(markers_present)

# Generate DotPlot
p_dot <- DotPlot(
  seu,
  features = markers_present,
  group.by = "cell_type_refined"
) +
  RotatedAxis() +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Canonical Proliferation Markers Across Cell Types",
    x = "",
    y = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

# Display plot
p_dot

# Save figure
ggsave(
  filename = file.path(
    outdir,
    "Figure8_Proliferation_DotPlot.png"
  ),
  plot = p_dot,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Slingshot trajectory analysis indicates that cells within the terminal repair-associated branch exhibit progressively increasing cell-cycle activity. Combined with the enrichment of canonical proliferation markers (Mki67, Top2a, Cdk1, Tyms, and Mcm5), these findings suggest that the proliferating mesenchymal cluster represents a cycling subset of late-stage repair cells rather than a distinct differentiation endpoint."