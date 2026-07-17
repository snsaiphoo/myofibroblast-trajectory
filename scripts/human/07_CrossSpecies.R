# Cross-species overlap coefficient analysis
# Compares the distribution of UCell pathway scores between
# mouse mesenchymal populations and human fibroblast populations.

library(Seurat)
library(dplyr)
library(pheatmap)

# Load data


mesenchymal <- readRDS("../../data/mesenchymal_geneset_scored_refined.rds")
human <- readRDS("../../data/human_fib_ucell.rds")

human <- subset(
  human,
  subset = condition == "Rupture"
)

###############################################################
# Shared pathways
###############################################################

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

###############################################################
# Display names
###############################################################

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

###############################################################
# Cell type order
###############################################################

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

###############################################################
# Overlap coefficient
###############################################################

overlap_coef <- function(x, y, n = 512){
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  if(length(x) < 5 | length(y) < 5){
    return(NA)
  }
  
  rng <- range(c(x, y))
  
  d1 <- density(
    x,
    from = rng[1],
    to = rng[2],
    n = n
  )
  
  d2 <- density(
    y,
    from = rng[1],
    to = rng[2],
    n = n
  )
  
  dx <- d1$x[2] - d1$x[1]
  
  sum(pmin(d1$y, d2$y)) * dx
}

###############################################################
# Output folder
###############################################################

outdir <- "../../results/07_snRNA/CrossSpecies_Overlap"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Generate one heatmap per pathway
###############################################################

for(path in score_order){
  
  overlap_matrix <- matrix(
    NA,
    nrow = length(mouse_levels),
    ncol = length(human_levels),
    dimnames = list(mouse_levels, human_levels)
  )
  
  for(mouse_cell in mouse_levels){
    
    mouse_scores <-
      mesenchymal@meta.data %>%
      filter(cell_type_refined == mouse_cell) %>%
      pull(!!sym(path))
    
    for(human_cell in human_levels){
      
      human_scores <-
        human@meta.data %>%
        filter(author_cell_type == human_cell) %>%
        pull(!!sym(path))
      
      overlap_matrix[mouse_cell, human_cell] <-
        overlap_coef(mouse_scores, human_scores)
      
    }
    
  }
  
  pheatmap(
    overlap_matrix,
    
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    
    color = colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(100),
    
    breaks = seq(0,1,length.out=101),
    
    display_numbers = round(overlap_matrix,2),
    
    number_color = "black",
    
    fontsize = 11,
    fontsize_number = 9,
    
    border_color = "grey80",
    
    main = pathway_names[path],
    
    filename = file.path(
      outdir,
      paste0(
        gsub("[/ ]","_",pathway_names[path]),
        "_overlap.png"
      )
    ),
    
    width = 6,
    height = 6
  )
  
}

# i7D mouse subset

mouse_i7d <- subset(
  mesenchymal,
  subset = condition == "I7D"
)

cell_counts <- table(mouse_i7d$cell_type_refined)

keep_cells <- names(cell_counts[cell_counts >= 50])

keep_cells

mouse_i7d <- subset(
  mouse_i7d,
  idents = keep_cells
)

mouse_avg <-
  mouse_i7d@meta.data %>%
  group_by(cell_type_refined) %>%
  summarise(across(all_of(score_order), mean, na.rm = TRUE))

human_rup <- subset(
  human,
  subset = condition == "Rupture"
)

human_avg <-
  human_rup@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(across(all_of(score_order), mean, na.rm = TRUE))

###############################################################
# Convert to matrices
###############################################################

mouse_matrix <- mouse_avg %>%
  column_to_rownames("cell_type_refined") %>%
  as.matrix()

human_matrix <- human_avg %>%
  column_to_rownames("author_cell_type") %>%
  as.matrix()

###############################################################
# Scale pathways within species
###############################################################

mouse_scaled <- scale(mouse_matrix)

human_scaled <- scale(human_matrix)

###############################################################
# Correlate mouse vs human
###############################################################

cor_matrix <- cor(
  t(mouse_scaled),
  t(human_scaled),
  method = "pearson"
)

###############################################################
# Output folder
###############################################################

dir.create(
  "../../results/07_snRNA/CrossSpecies_I7D_vs_Rupture",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Plot heatmap
###############################################################

pheatmap(
  cor_matrix,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  color = colorRampPalette(
    c("#313695",
      "#74add1",
      "white",
      "#f46d43",
      "#a50026")
  )(100),
  
  breaks = seq(-1, 1, length.out = 101),
  
  display_numbers = round(cor_matrix, 2),
  
  number_color = "black",
  
  border_color = "grey80",
  
  fontsize = 11,
  fontsize_number = 10,
  
  main = "Mouse I7D vs Human Rupture\nMechanotransduction Correlation",
  
  filename = "../../results/07_snRNA/CrossSpecies_I7D_vs_Rupture/I7D_vs_Rupture_Correlation.png",
  
  width = 7,
  height = 6
)
