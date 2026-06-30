library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(stringr)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)

# load data
mesenchymal <- readRDS("../data/mesenchymal_refined_annotated.rds")

cds <- as.CellDataSet(mesenchymal)
cds