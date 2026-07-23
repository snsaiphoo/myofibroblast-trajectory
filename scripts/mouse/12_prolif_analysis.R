# Assess Cell Type Composition Along the Repair Trajectory

library(dplyr)
library(ggplot2)

# cds_p2 was saved in workspace 
# Extract metadata and pseudotime
df <- as.data.frame(colData(cds_p2))
df$pseudotime <- pseudotime(cds_p2)

# Check pseudotime range
range(df$pseudotime)

# Divide pseudotime into 0.2-unit bins
df$pt_bin <- cut(
  df$pseudotime,
  breaks = seq(0, 1.8, by = 0.2),
  include.lowest = TRUE
)

# Cell type composition within each pseudotime bin

composition <- df %>%
  group_by(pt_bin, cell_type_refined) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pt_bin) %>%
  mutate(percent = 100 * n / sum(n))

print(composition, n = 100)

# Total number of cells in each pseudotime bin

cell_counts <- df %>%
  count(pt_bin)

print(cell_counts)

# Simplify cell types for visualization

composition_plot <- composition %>%
  mutate(
    CellType = case_when(
      cell_type_refined == "Repair fibroblasts" ~ "Repair fibroblasts",
      cell_type_refined == "Repair-activated stromal cells" ~ "Repair-activated stromal",
      cell_type_refined == "Proliferating mesenchymal cells" ~ "Proliferating",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(pt_bin, CellType) %>%
  summarise(percent = sum(percent), .groups = "drop")


# Position for cell count labels


cell_counts$label_y <- 103


# Plot
p <- ggplot(
  composition_plot,
  aes(
    x = pt_bin,
    y = percent,
    fill = CellType
  )
) +
  geom_col(width = 0.8) +
  
  geom_text(
    data = cell_counts,
    aes(
      x = pt_bin,
      y = label_y,
      label = paste0("n=",n)
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  
  scale_y_continuous(
    limits = c(0, 110),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Cell Composition Across Repair Pseudotime",
    subtitle = "n = total number of cells within each pseudotime bin",
    x = "Pseudotime Bin",
    y = "Cell Composition (%)",
    fill = "Cell Type"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size=10),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p)

#resave
ggsave(
  "../../results/02_figures_refined/Cell_Composition_Across_Pseudotime.png",
  width = 8,
  height = 5,
  dpi = 300
)

# LOESS curve per gene set
# Restricted to pseudotime <= 0.8

library(dplyr)
library(tidyr)
library(ggplot2)


# Gene sets

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

# Legend names

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

# Build dataframe

df <- as.data.frame(colData(cds_p2))

df$pseudotime <- pseudotime(cds_p2)

# Restrict to the portion of the trajectory before proliferating cells dominate
df <- df %>%
  filter(pseudotime <= 0.8)


# Prepare data for plotting

plot_df <- df %>%
  select(pseudotime, all_of(gene_sets)) %>%
  pivot_longer(
    cols = all_of(gene_sets),
    names_to = "GeneSet",
    values_to = "Score"
  ) %>%
  group_by(GeneSet) %>%
  mutate(
    Score_Z = as.numeric(scale(Score))
  ) %>%
  ungroup()

plot_df$GeneSet <- factor(
  plot_df$GeneSet,
  levels = names(legend_names),
  labels = legend_names
)

# LOESS plot

p <- ggplot(
  plot_df,
  aes(
    x = pseudotime,
    y = Score_Z,
    colour = GeneSet
  )
) +
  geom_smooth(
    method = "loess",
    span = 0.30,
    se = FALSE,
    linewidth = 1.3
  ) +
  theme_classic() +
  labs(
    title = "Mechanotransduction Pathway Dynamics Along Early Repair Pseudotime",
    subtitle = "Pseudotime ≤ 0.8",
    x = "Pseudotime",
    y = "Relative UCell Activity (Z-score)",
    colour = "Gene Set"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Save
#re save
ggsave(
  "../../results/03_monocle3_v2/partition3/Partition3_LOESS_Pseudotime_0.8.png",
  p,
  width = 9,
  height = 6,
  dpi = 600
)