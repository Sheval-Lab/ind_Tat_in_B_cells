## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
output_dir <- "results/tables/02_GSEA/02_KEGG"
fig_dir <- "results/figures/02_GSEA/02_KEGG"


## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


## Tat.stable vs Tat.16h -------------------------------------------------------
### Prepare pre-ranked gene list
gl_tatstable_tat16h <- 
  abs(log10(dgea_table_entrez$Tat.stable_vs_Tat.16h_padj)) * 
  sign(dgea_table_entrez$Tat.stable_vs_Tat.16h_log2FC)
names(gl_tatstable_tat16h) <- as.character(dgea_table_entrez$entrezID)
gl_tatstable_tat16h <- sort(gl_tatstable_tat16h, decreasing = TRUE)

### Run GSEA against KEGG database
### Use seed for reproducibility
set.seed(42) 
kegg_gsea_tatstable_tat16h <- gseKEGG(
  geneList = gl_tatstable_tat16h, 
  organism = "hsa", 
  # pvalueCutoff = 0.05, 
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Convert entrezID to symbol gene names
kegg_gsea_tatstable_tat16h <- setReadable(
  kegg_gsea_tatstable_tat16h,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID")

### Save GSEA results
saveRDS(kegg_gsea_tatstable_tat16h, file.path(output_dir, "kegg_gsea_tatstable_tat16h.rds"))

### Save enriched pathways
kegg_gsea_tatstable_tat16h@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.1) %>% 
  dplyr::select(ID, Description, NES, p.adjust, setSize, core_enrichment) %>%
  write_tsv(file.path(output_dir, "kegg_gsea_tatstable_tat16h_UP.padj01.tsv"))

kegg_gsea_tatstable_tat16h@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.1) %>% 
  dplyr::select(ID, Description, NES, p.adjust, setSize, core_enrichment) %>%
  write_tsv(file.path(output_dir, "kegg_gsea_tatstable_tat16h_DOWN.padj01.tsv"))


## Tat.16h vs Tat.0h -------------------------------------------------------
### Prepare pre-ranked gene list
gl_tat16h_tat0h <- 
  abs(log10(dgea_table_entrez$Tat.16h_vs_Tat.0h_padj)) * 
  sign(dgea_table_entrez$Tat.16h_vs_Tat.0h_log2FC)
names(gl_tat16h_tat0h) <- as.character(dgea_table_entrez$entrezID)
gl_tat16h_tat0h <- sort(gl_tat16h_tat0h, decreasing = TRUE)

### Run GSEA against KEGG database
### Use seed for reproducibility
set.seed(42) 
kegg_gsea_tat16h_tat0h <- gseKEGG(
  geneList = gl_tat16h_tat0h, 
  organism = "hsa", 
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Convert entrezID to symbol gene names
kegg_gsea_tat16h_tat0h <- setReadable(
  kegg_gsea_tat16h_tat0h,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID")

### Save GSEA results
saveRDS(kegg_gsea_tat16h_tat0h, file.path(output_dir, "kegg_gsea_tat16h_tat0h.rds"))

### Save enriched terms
kegg_gsea_tat16h_tat0h@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, Description, NES, p.adjust, setSize, core_enrichment) %>%
  write_tsv(file.path(output_dir, "kegg_gsea_tat16h_tat0h_UP.padj01.tsv"))

kegg_gsea_tat16h_tat0h@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, Description, NES, p.adjust, setSize, core_enrichment) %>%
  write_tsv(file.path(output_dir, "kegg_gsea_tat16h_tat0h_DOWN.padj01.tsv"))


## Tat.stable vs control -------------------------------------------------------
### Prepare pre-ranked gene list
gl_tatstable_control <- 
  abs(log10(dgea_table_entrez$Tat.stable_vs_control_padj)) * 
  sign(dgea_table_entrez$Tat.stable_vs_control_log2FC)
names(gl_tatstable_control) <- as.character(dgea_table_entrez$entrezID)
gl_tatstable_control <- sort(gl_tatstable_control, decreasing = TRUE)

### Run GSEA against KEGG database
### Use seed for reproducibility
set.seed(42) 
kegg_gsea_tatstable_control <- gseKEGG(
  geneList = gl_tatstable_control, 
  organism = "hsa", 
  pvalueCutoff = 1,
  verbose = FALSE,
  seed = TRUE)

### Convert entrezID to symbol gene names
kegg_gsea_tatstable_control <- setReadable(
  kegg_gsea_tatstable_control,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID")

### Save GSEA results - NO RESULTS
saveRDS(kegg_gsea_tatstable_control, file.path(output_dir, "kegg_gsea_tatstable_control.rds"))

### Save enriched terms
kegg_gsea_tatstable_control@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.1) %>% 
  dplyr::select(ID, Description, NES, p.adjust, setSize, core_enrichment) %>%
  write_tsv(file.path(output_dir, "kegg_gsea_tatstable_control_UP.padj01.tsv"))

kegg_gsea_tatstable_control@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.1) %>% 
  dplyr::select(ID, Description, NES, p.adjust, setSize, core_enrichment) %>%
  write_tsv(file.path(output_dir, "kegg_gsea_tatstable_control_DOWN.padj01.tsv"))

