## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
output_dir <- "results/tables/02_GSEA/03_Reactome"
fig_dir <- "results/figures/02_GSEA/03_Reactome"


## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


## Tat.stable vs Tat.16h -------------------------------------------------------
### Prepare pre-ranked gene list -----------------------------------------------
gl_tatstable_tat16h <- 
  abs(log10(dgea_table_entrez$Tat.stable_vs_Tat.16h_padj)) * 
  sign(dgea_table_entrez$Tat.stable_vs_Tat.16h_log2FC)
names(gl_tatstable_tat16h) <- as.character(dgea_table_entrez$entrezID)
gl_tatstable_tat16h <- sort(gl_tatstable_tat16h, decreasing = TRUE)

### Run GSEA against Reactome database -----------------------------------------
### Use seed for reproducibility
set.seed(42) 
reactome_gsea_tatstable_tat16h <- gsePathway(
  geneList = gl_tatstable_tat16h, 
  organism = "human", 
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Save GSEA results
saveRDS(reactome_gsea_tatstable_tat16h, file.path(output_dir, "reactome_gsea_tatstable_tat16h.rds"))

### Save enriched terms
reactome_gsea_tatstable_tat16h@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "reactome_gsea_tatstable_tat16h_UP.padj005.tsv"))

reactome_gsea_tatstable_tat16h@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "reactome_gsea_tatstable_tat16h_DOWN.padj005.tsv"))

### Emapplot -------------------------------------------------------------------
#### UP
up_reactome_gsea_tatstable_tat16h <- reactome_gsea_tatstable_tat16h %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES > 0, p.adjust < 0.01, nGenes >= 30) 

up_reactome_gsea_tatstable_tat16h_simm <- pairwise_termsim(up_reactome_gsea_tatstable_tat16h)
emapplot(up_reactome_gsea_tatstable_tat16h_simm, showCategory = 61, cex_label_category = 0.8)

ggsave(file.path(fig_dir, "REACTOME_GSEA_UP_Tat.stable_vs_Tat.16h.png"), units = "cm", width = 12, height = 10, dpi = 600, scaling = 0.35)

#### DOWN
down_reactome_gsea_tatstable_tat16h <- reactome_gsea_tatstable_tat16h %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES < 0, p.adjust < 0.05, nGenes >= 20) 

down_reactome_gsea_tatstable_tat16h_simm <- pairwise_termsim(down_reactome_gsea_tatstable_tat16h)
emapplot(down_reactome_gsea_tatstable_tat16h_simm, showCategory = 41, cex_label_category = 0.8)

ggsave(file.path(fig_dir, "REACTOME_GSEA_DOWN_Tat.stable_vs_Tat.16h.png"), units = "cm", width = 12, height = 10, dpi = 600, scaling = 0.35)


## Tat.16h vs Tat.0h -----------------------------------------------------------
### Prepare pre-ranked gene list -----------------------------------------------
gl_tat16h_tat0h <- 
  abs(log10(dgea_table_entrez$Tat.16h_vs_Tat.0h_padj)) * 
  sign(dgea_table_entrez$Tat.16h_vs_Tat.0h_log2FC)
names(gl_tat16h_tat0h) <- as.character(dgea_table_entrez$entrezID)
gl_tat16h_tat0h <- sort(gl_tat16h_tat0h, decreasing = TRUE)

### Run GSEA against Reactome database -----------------------------------------
### Use seed for reproducibility
set.seed(42) 
reactome_gsea_tat16h_tat0h <- gsePathway(
  geneList = gl_tat16h_tat0h, 
  organism = "human", 
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Save GSEA results
saveRDS(reactome_gsea_tat16h_tat0h, file.path(output_dir, "reactome_gsea_tat16h_tat0h.rds"))

### Save enriched terms
reactome_gsea_tat16h_tat0h@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "reactome_gsea_tat16h_tat0h_UP.padj005.tsv"))

reactome_gsea_tat16h_tat0h@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "reactome_gsea_tat16h_tat0h_DOWN.padj005.tsv"))

### Emapplot -------------------------------------------------------------------
#### UP
up_reactome_gsea_tat16h_tat0h <- reactome_gsea_tat16h_tat0h %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES > 0, p.adjust < 0.01) 

up_up_reactome_gsea_tat16h_tat0h_simm <- pairwise_termsim(up_reactome_gsea_tat16h_tat0h)
emapplot(up_up_reactome_gsea_tat16h_tat0h_simm, cex_label_category = 0.8)

ggsave(file.path(fig_dir, "REACTOME_GSEA_UP_Tat.16h_vs_Tat.0h.png"), units = "cm", width = 15, height = 10, dpi = 600, scaling = 0.7)

#### DOWN
down_reactome_gsea_tat16h_tat0h <- reactome_gsea_tat16h_tat0h %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES < 0, p.adjust < 0.05, nGenes >= 20) 

down_reactome_gsea_tat16h_tat0h_simm <- pairwise_termsim(down_reactome_gsea_tat16h_tat0h)
emapplot(down_reactome_gsea_tat16h_tat0h_simm, showCategory = 80, cex_label_category = 0.8)

ggsave(file.path(fig_dir, "REACTOME_GSEA_DOWN_Tat.16h_vs_Tat.0h.png"), units = "cm", width = 15, height = 10, dpi = 600, scaling = 0.35)



## Tat.stable vs control -------------------------------------------------------
### Prepare pre-ranked gene list -----------------------------------------------
gl_tatstable_control <- 
  abs(log10(dgea_table_entrez$Tat.stable_vs_control_padj)) * 
  sign(dgea_table_entrez$Tat.stable_vs_control_log2FC)
names(gl_tatstable_control) <- as.character(dgea_table_entrez$entrezID)
gl_tatstable_control <- sort(gl_tatstable_control, decreasing = TRUE)

### Run GSEA against Reactome database -----------------------------------------
### Use seed for reproducibility
set.seed(42) 
reactome_gsea_tatstable_control <- gsePathway(
  geneList = gl_tatstable_control, 
  organism = "human", 
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Save GSEA results
saveRDS(reactome_gsea_tatstable_control, file.path(output_dir, "reactome_gsea_tatstable_control.rds"))

### Save enriched terms
reactome_gsea_tatstable_control@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "reactome_gsea_tatstable_control_UP.padj005.tsv"))

reactome_gsea_tatstable_control@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "reactome_gsea_tatstable_control_DOWN.padj005.tsv"))


### Emapplot -------------------------------------------------------------------
# No enrichment...