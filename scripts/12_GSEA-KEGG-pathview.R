## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(pathview)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
fig_dir <- "results/figures/02_GSEA/02_KEGG"


## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


## Define a function to plot pathways from list --------------------------------
plot_kegg_map <- function(hsa_id, fc_list, suffix){
  pathview(
    gene.data = fc_list,
    pathway.id = as.character(hsa_id),
    out.suffix = as.character(suffix),
    species = "hsa",
    limit = list(gene = c(-1,1), cpd = 1),
    low = list(gene = "blue", cpd = "green"))
}


## Tat.stable vs control -------------------------------------------------------
### Prepare gene list of log2FC values Tat.stable vs control
fc_list <- dgea_table_entrez %>% filter(Tat.stable_vs_control_padj < 0.05) %>% pull(Tat.stable_vs_control_log2FC)
names(fc_list) <- dgea_table_entrez %>% filter(Tat.stable_vs_control_padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)


### Plot KEGG pathway maps  
hsa_tatstable_vs_control <- c(
  "hsa05164", "hsa05162", "hsa04620", "hsa05135", "hsa04260", "hsa04972",
  "hsa04261")

walk(hsa_tatstable_vs_control, safely(plot_kegg_map), fc_list = fc_list, suffix = "Tatstable_vs_control")


## Tat.16h vs Tat.0h -----------------------------------------------------------
### Prepare gene list of log2FC values Tat.16h vs Tat.0h
fc_list <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05) %>% pull(Tat.16h_vs_Tat.0h_log2FC)
names(fc_list) <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)

### Plot KEGG pathway maps  
hsa_tat16h_vs_tat0h <- c(
  "hsa00562", "hsa00970", "hsa01522", "hsa03008", "hsa03010", "hsa03013",
  "hsa03020", "hsa03030", "hsa03040", "hsa04014", "hsa04015", "hsa04024",
  "hsa04062", "hsa04070", "hsa04071", "hsa04072", "hsa04140", "hsa04270",
  "hsa04350", "hsa04360", "hsa04371", "hsa04510", "hsa04530", "hsa04540",
  "hsa04611", "hsa04612", "hsa04613", "hsa04625", "hsa04666", "hsa04670",
  "hsa04672", "hsa04713", "hsa04724", "hsa04725", "hsa04728", "hsa04750",
  "hsa04810", "hsa04910", "hsa04912", "hsa04916", "hsa04921", "hsa04924",
  "hsa04928", "hsa04935", "hsa05131", "hsa05135", "hsa05146", "hsa05165",
  "hsa05231", "hsa05414")

walk(hsa_tat16h_vs_tat0h, safely(plot_kegg_map), fc_list = fc_list, suffix = "Tat16h")


## Tat.stable vs Tat.16h -------------------------------------------------------
### Prepare gene list of log2FC values Tat.stable vs Tat.16h
fc_list <- dgea_table_entrez %>% filter(Tat.stable_vs_Tat.16h_padj < 0.05) %>% pull(Tat.stable_vs_Tat.16h_log2FC)
names(fc_list) <- dgea_table_entrez %>% filter(Tat.stable_vs_Tat.16h_padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)


hsa_tatstable_vs_tat16h <- c(
  "hsa00562", "hsa00600", "hsa00970", "hsa01521", "hsa01522", "hsa03008", 
  "hsa03010", "hsa03013", "hsa03015", "hsa03020", "hsa03030", "hsa03040",
  "hsa04014", "hsa04015", "hsa04020", "hsa04022", "hsa04024", "hsa04070",
  "hsa04071", "hsa04072", "hsa04110", "hsa04142", "hsa04260", "hsa04261",
  "hsa04270", "hsa04350", "hsa04360", "hsa04370", "hsa04371", "hsa04510",
  "hsa04512", "hsa04520", "hsa04530", "hsa04540", "hsa04611", "hsa04612",
  "hsa04666", "hsa04670", "hsa04713", "hsa04723", "hsa04724", "hsa04725",
  "hsa04728", "hsa04730", "hsa04742", "hsa04750", "hsa04810", "hsa04910",
  "hsa04911", "hsa04912", "hsa04916", "hsa04921", "hsa04923", "hsa04924",
  "hsa04928", "hsa04929", "hsa04935", "hsa04961", "hsa04970", "hsa04971",
  "hsa04972", "hsa05131", "hsa05169", "hsa05231", "hsa05410", "hsa05412",
  "hsa05414")


walk(hsa_tatstable_vs_tat16h, safely(plot_kegg_map), fc_list = fc_list, suffix = "Tatstable_vs_Tat16h")
