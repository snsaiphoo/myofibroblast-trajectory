# ============================================================
# UCell Heatmaps:
# Mean UCell Score + Cell-Type Composition (%)
#
# Colour = Mean UCell score
# Number = % of cells within each CONDITION
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(stringr)


# ============================================================
# 1. Define condition and cell-type order
# ============================================================

condition_order <- c(
  "WT",
  "I1D",
  "I7D",
  "I30D"
)

celltype_order <- c(
  "Mature tenocytes",
  "Stromal progenitor-like cells",
  "Signaling stromal cells",
  "Fibrochondrocyte-like tenocytes",
  "Activated tenocytes",
  "Proinflammatory mesenchymal cells",
  "Proliferating mesenchymal cells",
  "Repair-activated stromal cells",
  "Repair fibroblasts",
  "ECM-remodelling tenocytes"
)


# ============================================================
# 2. UCell score labels
# ============================================================

score_label_map <- c(
  "Integrin_FAK_UCell" =
    "Integrin–FAK",
  
  "RhoA_ROCK_UCell" =
    "RhoA–ROCK",
  
  "YAP_TEAD_UCell" =
    "YAP/TEAD",
  
  "HALLMARK_TGF_BETA_SIGNALING_UCell" =
    "TGFβ Signaling",
  
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell" =
    "SMAD Signaling",
  
  "GOBP_FIBROBLAST_ACTIVATION_UCell" =
    "Fibroblast Activation",
  
  "Myofibroblast_UCell" =
    "Myofibroblast",
  
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell" =
    "ECM Organization",
  
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell" =
    "Inflammation",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_UCell" =
    "Mesenchymal Transition"
)


# ============================================================
# 3. Gene sets to plot
# ============================================================

score_order <- c(
  "Integrin_FAK_UCell",
  "RhoA_ROCK_UCell",
  "YAP_TEAD_UCell",
  "HALLMARK_TGF_BETA_SIGNALING_UCell",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS_UCell",
  "GOBP_FIBROBLAST_ACTIVATION_UCell",
  "Myofibroblast_UCell",
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell",
  "HALLMARK_INFLAMMATORY_RESPONSE_UCell",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_UCell"
)


# ============================================================
# 4. Check that scores exist
# ============================================================

missing_scores <- setdiff(
  score_order,
  colnames(mesenchymal@meta.data)
)

if (length(missing_scores) > 0) {
  
  cat("Missing UCell scores:\n")
  print(missing_scores)
  
}

score_order <- intersect(
  score_order,
  colnames(mesenchymal@meta.data)
)


# ============================================================
# 5. Create output directory
# ============================================================

output_dir <- "../../results/02_figures_refined/UCell_condition"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ============================================================
# 6. Count cells by condition and cell type
# ============================================================

celltype_condition_counts <- mesenchymal@meta.data %>%
  filter(
    !is.na(condition),
    !is.na(cell_type_refined),
    condition %in% condition_order,
    cell_type_refined %in% celltype_order
  ) %>%
  count(
    condition,
    cell_type_refined,
    name = "n_cells"
  )


# ============================================================
# 7. Calculate TOTAL number of cells in each condition
#
# This is the denominator.
#
# WT:
#     total WT cells = 300
#
# Repair fibroblasts:
#     3 / 300 = 1%
# ============================================================

condition_totals <- mesenchymal@meta.data %>%
  filter(
    !is.na(condition),
    condition %in% condition_order
  ) %>%
  count(
    condition,
    name = "total_cells"
  )


# ============================================================
# 8. Calculate percentage of each condition
# ============================================================

celltype_condition_percent <- celltype_condition_counts %>%
  left_join(
    condition_totals,
    by = "condition"
  ) %>%
  mutate(
    percent = n_cells / total_cells * 100
  )


# ============================================================
# 9. Make percentage matrix
# ============================================================

percent_matrix <- celltype_condition_percent %>%
  select(
    cell_type_refined,
    condition,
    percent
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = percent,
    values_fill = 0
  ) %>%
  column_to_rownames(
    "cell_type_refined"
  ) %>%
  as.matrix()


# Make sure all conditions exist

missing_conditions <- setdiff(
  condition_order,
  colnames(percent_matrix)
)

if (length(missing_conditions) > 0) {
  
  for (cond in missing_conditions) {
    percent_matrix[, cond] <- 0
  }
  
}


# Apply condition order

percent_matrix <- percent_matrix[
  ,
  condition_order,
  drop = FALSE
]


# Apply cell-type order

percent_matrix <- percent_matrix[
  intersect(
    celltype_order,
    rownames(percent_matrix)
  ),
  ,
  drop = FALSE
]


# ============================================================
# 10. Convert percentages into labels
# ============================================================

percent_labels <- matrix(
  "",
  nrow = nrow(percent_matrix),
  ncol = ncol(percent_matrix),
  dimnames = dimnames(percent_matrix)
)


for (i in seq_len(nrow(percent_matrix))) {
  
  for (j in seq_len(ncol(percent_matrix))) {
    
    value <- percent_matrix[i, j]
    
    if (value > 0 && value < 0.5) {
      
      percent_labels[i, j] <- "<0.5%"
      
    } else {
      
      percent_labels[i, j] <- paste0(
        round(value, 1),
        "%"
      )
      
    }
    
  }
  
}


# ============================================================
# 11. Create one heatmap per gene set
# ============================================================

for (score in score_order) {
  
  
  # ----------------------------------------------------------
  # Calculate mean UCell score
  # ----------------------------------------------------------
  
  score_data <- mesenchymal@meta.data %>%
    
    filter(
      !is.na(condition),
      !is.na(cell_type_refined),
      condition %in% condition_order,
      cell_type_refined %in% celltype_order
    ) %>%
    
    group_by(
      cell_type_refined,
      condition
    ) %>%
    
    summarise(
      mean_score = mean(
        .data[[score]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    
    pivot_wider(
      names_from = condition,
      values_from = mean_score
    ) %>%
    
    column_to_rownames(
      "cell_type_refined"
    ) %>%
    
    as.matrix()
  
  
  # ----------------------------------------------------------
  # Add missing conditions if necessary
  # ----------------------------------------------------------
  
  missing_conditions <- setdiff(
    condition_order,
    colnames(score_data)
  )
  
  if (length(missing_conditions) > 0) {
    
    for (cond in missing_conditions) {
      score_data[, cond] <- NA
    }
    
  }
  
  
  # ----------------------------------------------------------
  # Apply condition order
  # ----------------------------------------------------------
  
  score_data <- score_data[
    ,
    condition_order,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Apply cell-type order
  # ----------------------------------------------------------
  
  score_data <- score_data[
    intersect(
      celltype_order,
      rownames(score_data)
    ),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Match percentage labels to UCell matrix
  # ----------------------------------------------------------
  
  current_rows <- rownames(score_data)
  
  current_percent_labels <- percent_labels[
    current_rows,
    condition_order,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Get gene-set label
  # ----------------------------------------------------------
  
  if (score %in% names(score_label_map)) {
    
    gene_set_label <- score_label_map[[score]]
    
  } else {
    
    gene_set_label <- score
    
  }
  
  
  # ----------------------------------------------------------
  # Clean filename
  # ----------------------------------------------------------
  
  filename_clean <- score %>%
    str_remove("_UCell$") %>%
    str_replace_all(
      "[^A-Za-z0-9]+",
      "_"
    )
  
  
  # ----------------------------------------------------------
  # Create heatmap
  # ----------------------------------------------------------
  
  pheatmap(
    score_data,
    
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    
    # DO NOT SCALE
    scale = "none",
    
    border_color = NA,
    
    fontsize_row = 10,
    fontsize_col = 11,
    
    angle_col = 45,
    
    # Title
    main = paste0(
      gene_set_label,
      "\nColour = Mean UCell Score; Values = Cell-Type Composition (%) per Condition"
    ),
    
    # Percentage inside each cell
    display_numbers = current_percent_labels,
    
    number_format = "%s",
    
    fontsize_number = 9,
    
    # Save
    filename = paste0(
      output_dir,
      "/",
      filename_clean,
      "_Condition_CellType_Percent_Heatmap.png"
    ),
    
    width = 10,
    height = 7
  )
  
}


# ============================================================
# 12. Confirmation
# ============================================================

cat("\nHeatmaps successfully created!\n")
cat("Figures saved to:\n")
cat(output_dir, "\n")