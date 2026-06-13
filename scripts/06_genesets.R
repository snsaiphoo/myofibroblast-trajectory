library(Seurat)
library(msigdbr)
library(dplyr)
library(stringr)

# Load mesenchymal object
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")

# Load mouse msigdbr gene sets
msig_mouse <- msigdbr(db_species = "MM", species = "Mus musculus") 

# View available collections 
unique(msig_mouse$gs_collection) 

# Keep MH, M2, and M5
all_gs <- msig_mouse %>% filter(str_detect(gs_collection, "MH|M2|M5"))

# Group the msigdbr gene sets
selected_msigdb_sets <- c(
  # Hallmark MH
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  # Reactome MS:CP
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_INTEGRIN_SIGNALING",
  
  # GO Biological Process M5:GO:BP
  "GOBP_SMAD_PROTEIN_SIGNAL_TRANSDUCTION",
  "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION",
  "GOBP_FOCAL_ADHESION_ASSEMBLY",
  "GOBP_RESPONSE_TO_MECHANICAL_STIMULUS",
  "GOBP_WOUND_HEALING",
  "GOBP_FIBROBLAST_ACTIVATION",
  "GOBP_MYOFIBROBLAST_DIFFERENTIATION",
  "GOBP_FIBROBLAST_PROLIFERATION",
  "GOBP_TENDON_DEVELOPMENT",
  "GOBP_TENDON_FORMATION"
)

# find the listed gene sets in the database dataframe
msigdb_gene_sets <- all_gs %>%
  filter(gs_name %in% selected_msigdb_sets) %>%
  split(x = .$gene_symbol, f = .$gs_name)

missing_sets <- setdiff(selected_msigdb_sets, names(msigdb_gene_sets))
missing_sets

# Custom gene sets for YAP/TAZ, integrin/FAK, and myofibroblast
yap_taz <- c("Ctgf", "Cyr61", "Ankrd1", "Amotl2", "Ccnd1", "Birc5", "Serpine1", "Lats2")

integrin_fak <- c("Itga1", "Itga2", "Itga5", "Itgav", "Itgb1", "Ptk2", "Pxn", "Vcl", "Tln1", "Zyx")

myofibroblast <- c("Acta2", "Tagln", "Cnn1", "Myl9", "Cald1", "Fn1", "Col1a1", "Col1a2", "Postn", "Tnc", "Thbs1", "Ctgf", "Serpine1", "Timp1", "S100a4")

custom_gene_sets <- list(
  YAP_TAZ_CURATED = yap_taz,
  INTEGRIN_FAK_CURATED = integrin_fak,
  MYOFIBROBLAST_CURATED = myofibroblast
)

# combine both gene sets
gene_sets <- c(
  msigdb_gene_sets,
  custom_gene_sets
)

# Keep only genes present in the dataset
gene_sets_present <- lapply(
  gene_sets,
  function(x) intersect(unique(x), rownames(mesenchymal))
)

# Group gene sets into biology categories
gene_set_groups <- list(
  
  Mechanotransduction = c(
    "YAP_TAZ_CURATED",
    "REACTOME_INTEGRIN_SIGNALING",
    "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
    "GOBP_FOCAL_ADHESION_ASSEMBLY",
    "GOBP_RESPONSE_TO_MECHANICAL_STIMULUS",
    "INTEGRIN_FAK_CURATED"
  ),
  
  TGFB_SMAD = c(
    "HALLMARK_TGF_BETA_SIGNALING",
    "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
    "GOBP_SMAD_PROTEIN_SIGNAL_TRANSDUCTION"
  ),
  
  ECM_Fibrosis = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "REACTOME_ECM_PROTEOGLYCANS",
    "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
    "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
    "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION",
    "GOBP_MYOFIBROBLAST_DIFFERENTIATION",
    "MYOFIBROBLAST_CURATED"
  ),
  
  Inflammation = c(
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  ),
  
  Repair_Tendon_State = c(
    "GOBP_WOUND_HEALING",
    "GOBP_FIBROBLAST_ACTIVATION",
    "GOBP_FIBROBLAST_PROLIFERATION",
    "GOBP_TENDON_DEVELOPMENT",
    "GOBP_TENDON_FORMATION"
  )
)

# Keep only grouped sets that survived filtering
gene_set_groups_present <- lapply(
  gene_set_groups,
  function(x) intersect(x, names(gene_sets_present))
)


# Quality checks
gene_set_sizes <- sort(sapply(gene_sets_present, length))
gene_set_sizes

missing_grouped_sets <- setdiff(
  unlist(gene_set_groups),
  names(gene_sets_present)
)

ungrouped_sets <- setdiff(
  names(gene_sets_present),
  unlist(gene_set_groups)
)

missing_grouped_sets
ungrouped_sets


saveRDS(
  gene_sets_present,
  "../data/mechanotransduction_gene_sets.rds"
)

saveRDS(
  gene_set_groups_present,
  "../data/mechanotransduction_gene_set_groups.rds"
)

