library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

#-----------------------------
# Load scored fibroblast object
#-----------------------------
fib <- readRDS("../../data/human_fib_ucell.rds")

#-----------------------------
# Check conditions
#-----------------------------
table(fib$condition)

fib$condition <- factor(
  fib$condition,
  levels = c("Healthy", "Rupture")
)

#-----------------------------
# Pathway groups
#-----------------------------
mech_paths <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell"
)

activation_paths <- c(
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell"
)

ecm_paths <- c(
  "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
)

#-----------------------------
# Rename pathways
#-----------------------------
rename_paths <- function(x) {
  dplyr::recode(
    x,
    Integrin_FAK_UCell = "Integrin-FAK",
    RhoA_ROCK_UCell = "RhoA-ROCK",
    YAP_TEAD_UCell = "YAP/TEAD",
    HALLMARK_TGF_BETA_SIGNALING_UCell = "TGFβ",
    REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell = "SMAD",
    GOBP_FIBROBLAST_ACTIVATION_UCell = "Fibroblast Activation",
    Myofibroblast_UCell = "Myofibroblast",
    GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION_UCell = "ECM Organization",
    HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammation"
  )
}

#-----------------------------
# Plotting function
#-----------------------------
plot_ucell <- function(pathways, title, stats_filename = NULL){
  
  #-----------------------------
  # Extract data
  #-----------------------------
  
  df <- FetchData(
    fib,
    vars = c("condition", pathways)
  )
  
  df_long <- df %>%
    pivot_longer(
      cols = -condition,
      names_to = "Pathway",
      values_to = "Score"
    )
  
  #-----------------------------
  # Rename pathways
  #-----------------------------
  
  df_long$Pathway <- dplyr::recode(
    df_long$Pathway,
    
    Integrin_FAK_UCell = "Integrin-FAK",
    RhoA_ROCK_UCell = "RhoA-ROCK",
    YAP_TEAD_UCell = "YAP/TEAD",
    HALLMARK_TGF_BETA_SIGNALING_UCell = "TGFβ",
    REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell = "SMAD",
    
    GOBP_FIBROBLAST_ACTIVATION_UCell = "Fibroblast Activation",
    Myofibroblast_UCell = "Myofibroblast",
    
    GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION_UCell = "ECM Remodeling",
    HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammation"
  )
  
  desired_order <- c(
    "Integrin-FAK",
    "RhoA-ROCK",
    "YAP/TEAD",
    "TGFβ",
    "SMAD",
    "Fibroblast Activation",
    "Myofibroblast",
    "ECM Remodeling",
    "Inflammation"
  )
  
  df_long$Pathway <- factor(
    df_long$Pathway,
    levels = desired_order
  )
  
  #====================================================
  # Statistics for plotting
  #====================================================
  
  stats_plot <- df_long %>%
    group_by(Pathway) %>%
    wilcox_test(Score ~ condition) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "condition")
  
  #====================================================
  # Statistics for CSV
  #====================================================
  
  stats_export <- df_long %>%
    group_by(Pathway) %>%
    wilcox_test(Score ~ condition) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    ungroup()
  
  #-----------------------------
  # Save CSV
  #-----------------------------
  
  if(!is.null(stats_filename)){
    
    stats_export <- as.data.frame(stats_export)
    
    # Remove any list columns
    keep <- !sapply(stats_export, is.list)
    stats_export <- stats_export[, keep]
    
    write.csv(
      stats_export,
      stats_filename,
      row.names = FALSE
    )
    
  }
  
  #====================================================
  # Plot
  #====================================================
  
  p <- ggplot(
    df_long,
    aes(
      condition,
      Score,
      fill = condition
    )
  ) +
    
    geom_violin(
      trim = FALSE,
      alpha = 0.85,
      colour = "black",
      linewidth = 0.8
    ) +
    
    geom_boxplot(
      width = 0.12,
      fill = NA,
      colour = "black",
      outlier.shape = NA
    ) +
    
    stat_pvalue_manual(
      stats_plot,
      label = "p.adj.signif",
      tip.length = 0.01,
      size = 6
    ) +
    
    facet_wrap(
      ~Pathway,
      scales = "free_y",
      nrow = 1
    ) +
    
    scale_fill_manual(
      values = c(
        Healthy = "#6CC36C",
        Rupture = "#EF4444"
      )
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.03,0.15))
    ) +
    
    labs(
      title = title,
      x = NULL,
      y = "UCell Score"
    ) +
    
    theme_classic(base_size = 18) +
    
    theme(
      legend.position = "none",
      
      strip.background = element_rect(
        fill = "grey96",
        colour = "black"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 18
      ),
      
      axis.text = element_text(
        face = "bold",
        size = 16,
        colour = "black"
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 20
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 24,
        hjust = 0.5
      )
    )
  
  return(p)
  
}

#-----------------------------
# Create figures
#-----------------------------
p1 <- plot_ucell(
  mech_paths,
  "Mechanotransduction Signaling",
  "../../results/Figure5A_Mechanotransduction_Wilcoxon.csv"
)

p2 <- plot_ucell(
  activation_paths,
  "Fibroblast Activation",
  "../../results/Figure5B_FibroblastActivation_Wilcoxon.csv"
)

p3 <- plot_ucell(
  ecm_paths,
  "ECM Remodeling and Inflammation",
  "../../results/Figure5C_ECMInflammation_Wilcoxon.csv"
)

#-----------------------------
# Save
#-----------------------------
ggsave(
  "../../results/07_snRNA/Figure5A_Mechanotransduction_Violin.png",
  p1,
  width = 13,
  height = 6,
  dpi = 600
)

ggsave(
  "../../results/07_snRNA/Figure5B_FibroblastActivation_Violin.png",
  p2,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "../../results/07_snRNA/Figure5C_ECMInflammation_Violin.png",
  p3,
  width = 8,
  height = 6,
  dpi = 600
)
