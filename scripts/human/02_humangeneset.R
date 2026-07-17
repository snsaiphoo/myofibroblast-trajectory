library(Seurat)
library(msigdbr)
library(dplyr)
library(stringr)

# load human msigdb gene sets
msig_human <- msigdbr(species = "Homo sapiens")

# keep Hallmark, Reactome, and GO Biological Process collections
all_gs <- msig_human %>%
  filter(str_detect(gs_collection, "H|C2|C5"))

# selected MSigDB gene sets
selected_msigdb_sets <- c(
  "HALLMARK_TGF_BETA_SIGNALING",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "GOBP_FIBROBLAST_ACTIVATION"
)

# extract MSigDB gene sets as Ensembl IDs
msigdb_gene_sets <- all_gs %>%
  filter(gs_name %in% selected_msigdb_sets) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

# ----------------------------
# Custom human gene sets
# ----------------------------

yap_taz_signature <- c(
  "CTGF", "CYR61", "ANKRD1",
  "AMOTL2", "LATS2",
  "CCND1", "BIRC5"
)

Integrin_FAK_Core <- c(
  "ITGB1", "ILK", "ITGA5", "ITGAV",
  "PTK2", "VCL", "PXN", "TLN1", "ZYX"
)

RhoA_ROCK_Actomyosin <- c(
  "RHOA", "ROCK1", "ROCK2",
  "MYL9", "MYL12A", "MYL12B",
  "MYH9", "MYH10",
  "LIMK1", "LIMK2",
  "CFL1",
  "DIAPH1", "DIAPH2",
  "ARHGEF1", "ARHGEF2",
  "ARHGAP1"
)

myofibroblast <- c(
  "ACTA2",
  "TNC",
  "POSTN",
  "FN1"
)

custom_sets <- list(
  YAP_TEAD = yap_taz_signature,
  Integrin_FAK = Integrin_FAK_Core,
  RhoA_ROCK = RhoA_ROCK_Actomyosin,
  Myofibroblast = myofibroblast
)

# ----------------------------
# Convert custom gene sets to Ensembl IDs
# ----------------------------

symbol_to_ensembl <- msig_human %>%
  distinct(gene_symbol, ensembl_gene)

custom_sets_ensembl <- lapply(custom_sets, function(genes) {
  symbol_to_ensembl %>%
    filter(gene_symbol %in% genes) %>%
    pull(ensembl_gene) %>%
    unique()
})

# ----------------------------
# Combine all gene sets
# ----------------------------

all_gene_sets <- c(custom_sets_ensembl, msigdb_gene_sets)

# remove duplicates and missing IDs
all_gene_sets <- lapply(all_gene_sets, function(x) {
  unique(na.omit(x))
})

# ----------------------------
# Check overlap between gene sets
# ----------------------------

repeat_genes <- data.frame(
  set1 = character(),
  set2 = character(),
  n_shared = integer(),
  shared_genes = character(),
  stringsAsFactors = FALSE
)

for(i in 1:(length(all_gene_sets)-1)){
  for(j in (i+1):length(all_gene_sets)){
    
    shared <- intersect(all_gene_sets[[i]], all_gene_sets[[j]])
    
    repeat_genes <- rbind(
      repeat_genes,
      data.frame(
        set1 = names(all_gene_sets)[i],
        set2 = names(all_gene_sets)[j],
        n_shared = length(shared),
        shared_genes = paste(shared, collapse = ", "),
        stringsAsFactors = FALSE
      )
    )
  }
}

repeat_genes <- repeat_genes[order(-repeat_genes$n_shared), ]

# ----------------------------
# Save
# ----------------------------

saveRDS(all_gene_sets, "../../data/all_gene_sets_human_ensembl.rds")
write.csv(
  repeat_genes,
  "../../results/repeat_genes_human.csv",
  row.names = FALSE
)
