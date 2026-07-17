# Cross species comparisons 
library(Seurat)
library(dplyr)
library(tibble)
library(pheatmap)
library(forcats)
library(tidyr)
library(ggplot2)

mesenchymal <- readRDS("../../data/mesenchymal_geneset_scored_refined.rds")
m_gene_sets <- readRDS("../../data/all_gene_sets_refined.rds")
human <- readRDS("../../data/human_fib_ucell.rds")
h_gene_sets <- readRDS("../../data/all_gene_sets_human_ensembl.rds")

# Mouse pathways (shared between species)
score_order <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

# Mouse
mouse_avg <-
  mesenchymal@meta.data %>%
  group_by(cell_type_refined) %>%  
  summarise(across(all_of(score_order), mean, na.rm = TRUE))

mouse_matrix <-
  mouse_avg %>%
  column_to_rownames("cell_type_refined") %>%
  as.matrix()

# Human

human_avg <-
  human@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(across(all_of(score_order), mean, na.rm = TRUE))

human_matrix <-
  human_avg %>%
  column_to_rownames("author_cell_type") %>%
  as.matrix()

mouse_scaled <- scale(mouse_matrix)
human_scaled <- scale(human_matrix)

# Figure 1: 
# Figure per pathway

# Mouse averages
mouse_avg <- mesenchymal@meta.data %>%
  group_by(cell_type_refined) %>%
  summarise(across(all_of(score_order), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-cell_type_refined, names_to = "Pathway", values_to = "UCell") %>%
  mutate(Species = "Mouse",
         CellType = cell_type_refined) %>%
  select(Species, CellType, Pathway, UCell)

# Human averages
human_avg <- fib@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(across(all_of(score_order), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-author_cell_type, names_to = "Pathway", values_to = "UCell") %>%
  mutate(Species = "Human",
         CellType = author_cell_type) %>%
  select(Species, CellType, Pathway, UCell)

combined <- bind_rows(mouse_avg, human_avg)

mouse_levels <- c(
  "Activated tenocytes",
  "ECM-remodelling tenocytes",
  "Fibrochondrocyte-like tenocytes",
  "Mature tenocytes",
  "Repair fibroblasts",
  "Repair-activated stromal cells",
  "Proliferating mesenchymal cells",
  "Stromal progenitor-like cells",
  "Signaling stromal cells",
  "Proinflammatory mesenchymal cells"
)

human_levels <- c(
  "ABCA10hi fibroblasts",
  "ADAM12hi fibroblasts",
  "FBLN1hi fibroblasts",
  "NR4A1hi fibroblasts"
)

combined <- combined %>%
  mutate(
    Display = case_when(
      Species == "Mouse" ~ paste0("1_", CellType),
      TRUE ~ paste0("2_", CellType)
    )
  )

display_levels <- c(
  paste0("1_", mouse_levels),
  paste0("2_", human_levels)
)

combined$Display <- factor(combined$Display,
                           levels = rev(display_levels))

pathway_names <- c(
  Integrin_FAK_UCell = "Integrin–FAK",
  RhoA_ROCK_UCell = "RhoA–ROCK",
  YAP_TEAD_UCell = "YAP/TEAD",
  HALLMARK_TGF_BETA_SIGNALING_UCell = "TGF-β",
  REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell = "SMAD",
  GOBP_FIBROBLAST_ACTIVATION_UCell = "Fibroblast Activation",
  Myofibroblast_UCell = "Myofibroblast",
  HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammatory Response"
)

dir.create("../../results/07_snRNA/CrossSpecies_PathwayHeatmaps",
           recursive = TRUE,
           showWarnings = FALSE)

for(path in score_order){
  
  df <- combined %>% filter(Pathway == path)
  lims <- range(df$UCell, na.rm = TRUE)
  
  p <- ggplot(df,
              aes(x = Species,
                  y = Display,
                  fill = UCell)) +
    
    geom_tile(width = 0.9,
              height = 0.9,
              colour = "white",
              linewidth = 0.5) +
    
    scale_fill_viridis_c(
      option = "plasma",
      limits = lims,
      name = "Mean\nUCell"
    ) +
    
    scale_y_discrete(
      labels = function(x) sub("^[12]_", "", x)
    ) +
    
    labs(
      title = pathway_names[path],
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(base_size = 14) +
    
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 10),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold",
                                hjust = 0.5,
                                size = 18),
      legend.position = "right"
    ) +
    
    annotate(
      "segment",
      x = 0.5,
      xend = 2.5,
      y = 4.5,
      yend = 4.5,
      linewidth = 0.8
    )
  
  print(p)
  
  ggsave(
    filename = file.path(
      "../../results/07_snRNA/CrossSpecies_PathwayHeatmaps",
      paste0(gsub("[/ ]", "_", pathway_names[path]), ".png")
    ),
    plot = p,
    width = 5.5,
    height = 8,
    dpi = 300
  )
}
