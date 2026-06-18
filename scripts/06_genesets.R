library(Seurat)
library(msigdbr)
library(dplyr)
library(stringr)

# load data
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")

# load mouse msigdb gene sets
msig_mouse <- msigdbr(db_species = "MM", species = "Mus musculus")

# keep MH, M2, and M5 collections
all_gs <- msig_mouse %>% filter(str_detect(gs_collection, "MH|M2|M5"))

# select gene sets
selected_msigdb_sets <- c(
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_INTEGRIN_SIGNALING",
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

# extract gene sets
msigdb_gene_sets <- all_gs %>%
  filter(gs_name %in% selected_msigdb_sets) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# check for missing sets
missing_sets <- setdiff(selected_msigdb_sets, names(msigdb_gene_sets))
missing_sets

# custom gene sets
yap_taz <- c(
  "Ctgf", "Cyr61", "Ankrd1", "Amotl2", "Birc5", "Lats2", "Mmp3"
)

integrin_fak <- c(
  "Itgb1", "Ilk", "Itga1", "Itga2", "Itga5", "Itgav",
  "Ptk2", "Vcl", "Pxn", "Tln1", "Zyx", "Fblim1"
)

myofibroblast <- c(
  "Acta2", "Col1a1", "Col1a2", "Col3a1",
  "Tnc", "Postn", "Fmod", "Scx", "Fn1"
)

custom_gene_sets <- list(
  YAP_TAZ_CURATED       = yap_taz,
  INTEGRIN_FAK_CURATED  = integrin_fak,
  MYOFIBROBLAST_CURATED = myofibroblast
)

# universal fibroblast marker diagnostic
VlnPlot(
  mesenchymal,
  features = c("Pi16", "Dpp4", "Col15a1", "Cd55", "Procr"),
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
)

ggsave(
  "../figures/diagnostic_universal_fibroblast_markers.png",
  last_plot(),
  width = 14, height = 8, dpi = 300
)

# myofibroblast gene set expression diagnostic
VlnPlot(
  mesenchymal,
  features = c("Scx", "Fmod", "Thbs4", "Postn", "Tagln"),
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
)

ggsave(
  "../figures/diagnostic_myofibroblast_geneset_check.png",
  last_plot(),
  width = 14, height = 8, dpi = 300
)

# combine msigdb and custom gene sets
gene_sets <- c(msigdb_gene_sets, custom_gene_sets)

# filter to genes present in dataset
gene_sets_present <- lapply(
  gene_sets,
  function(x) intersect(unique(x), rownames(mesenchymal))
)

# group gene sets by biology
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
  
  ECM_Remodeling = c(
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

# filter groups to sets present in dataset
gene_set_groups_present <- lapply(
  gene_set_groups,
  function(x) intersect(x, names(gene_sets_present))
)

# quality checks
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

# save
saveRDS(gene_sets_present, "../data/mechanotransduction_gene_sets.rds")
saveRDS(gene_set_groups_present, "../data/mechanotransduction_gene_set_groups.rds")
