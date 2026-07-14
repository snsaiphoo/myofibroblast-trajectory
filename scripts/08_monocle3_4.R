# loess curve per gene set
# script was adjusted for the main 3 partitions manually 

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

# cleaned up names

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

# build the dataframe for plotting
df <- as.data.frame(colData(cds_p3))

df$pseudotime <- pseudotime(cds_p3)

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
    x = "Pseudotime",
    y = "Relative UCell Activity (Z-score)",
    colour = "Gene Set"
  ) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

# change directory for each partition
ggsave(
  "../results/monocle3_v2/partition3/Partition3_LOESS_AllGeneSets.png",
  p,
  width = 9,
  height = 6,
  dpi = 300
)