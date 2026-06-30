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

# reduced selected gene sets from msigdb
selected_msigdb_sets <- c(
  # TGFβ / SMAD
  "HALLMARK_TGF_BETA_SIGNALING",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  
  # ECM remodeling
  "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION",
  
  # Inflammation
  "HALLMARK_INFLAMMATORY_RESPONSE",
  
  # Repair / fibroblast activation
  "GOBP_FIBROBLAST_ACTIVATION"
)

# extract gene sets
msigdb_gene_sets <- all_gs %>%
  filter(gs_name %in% selected_msigdb_sets) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Customized gene sets
yap_taz_signature <- c("Ctgf", "Cyr61", "Ankrd1", "Amotl2", "Lats2", "Ccnd1", "Birc5")

Integrin_FAK_Core <- c("Itgb1", "Ilk", "Itga5", "Itgav", "Ptk2", "Vcl", "Pxn", "Tln1", "Zyx")

RhoA_ROCK_Actomyosin <- c("Rhoa", "Rock1", "Rock2", "Myl9", "Myl12a", "Myl12b", "Myh9", "Myh10", "Limk1", "Limk2", "Cfl1", "Diaph1", "Diaph2", "Arhgef1", "Arhgef2", "Arhgap1")

myofibroblast <- c("Acta2", "Tnc", "Postn", "Fn1")

custom_sets <- list(
  YAP_TEAD = yap_taz_signature,
  Integrin_FAK = Integrin_FAK_Core,
  RhoA_ROCK = RhoA_ROCK_Actomyosin,
  Myofibroblast = myofibroblast
)

# Check for overlap in gene sets 
all_gene_sets <- c(custom_sets, msigdb_gene_sets)
all_gene_sets <- lapply(all_gene_sets, unique)

repeat_genes <- data.frame(
  set1 = character(),
  set2 = character(),
  n_shared = integer(),
  shared_genes = character(),
  stringsAsFactors = FALSE
)

# for loop to check overlap in gene sets 
for(i in 1:(length(all_gene_sets) - 1)){
  for(j in (i+1):length(all_gene_sets)){
    shared <- intersect(all_gene_sets[[i]], all_gene_sets[[j]])
    new_row <- data.frame(set1 = names(all_gene_sets)[i], set2 = names(all_gene_sets)[j], n_shared = length(shared), shared_genes = paste(shared, collapse = ", ") )
    repeat_genes <- rbind(repeat_genes, new_row)
  }
}
  
repeat_genes <- repeat_genes[order(-repeat_genes$n_shared), ]


saveRDS(all_gene_sets, "../data/all_gene_sets_refined.rds")

