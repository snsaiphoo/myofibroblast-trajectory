library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(scCustomize)
library(stringr)

# Load in data 
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")
gene_sets <- readRDS(all_gene_sets, "../data/all_gene_sets_refined.rds")