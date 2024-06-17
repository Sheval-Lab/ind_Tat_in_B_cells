## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
output_dir <- "results/tables/02_GSEA/01_GOBP"
fig_dir <- "results/figures/02_GSEA/01_GOBP"


## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


## Tat.stable vs Tat.16h -------------------------------------------------------
### Prepare pre-ranked gene list
gl_tatstable_tat16h <- 
  abs(log10(dgea_table_entrez$Tat.stable_vs_Tat.16h_padj)) * 
  sign(dgea_table_entrez$Tat.stable_vs_Tat.16h_log2FC)
names(gl_tatstable_tat16h) <- as.character(dgea_table_entrez$entrezID)
gl_tatstable_tat16h <- sort(gl_tatstable_tat16h, decreasing = TRUE)

### Run GSEA against GO-BP database
### Use seed for reproducibility
set.seed(42) 
gobp_gsea_tatstable_tat16h <- gseGO(
  geneList = gl_tatstable_tat16h, 
  ont = "BP", 
  OrgDb = "org.Hs.eg.db",
  # pvalueCutoff = 0.05, 
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Save GSEA results
saveRDS(gobp_gsea_tatstable_tat16h, file.path(output_dir, "gobp_gsea_tatstable_tat16h.rds"))

### Simplify GO results
gobp_gsea_tatstable_tat16h_simple <- gobp_gsea_tatstable_tat16h %>% 
  # Filter GO terms by adjusted pvalue
  filter(p.adjust < 0.05) %>% 
  simplify()

saveRDS(gobp_gsea_tatstable_tat16h_simple, file.path(output_dir, "gobp_gsea_tatstable_tat16h_simple.rds"))

### Save enriched terms
gobp_gsea_tatstable_tat16h_simple@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "gobp_gsea_tatstable_tat16h_simple_UP.padj005.tsv"))

gobp_gsea_tatstable_tat16h_simple@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "gobp_gsea_tatstable_tat16h_simple_DOWN.padj005.tsv"))

### Treeplot
#### UP
up_gobp_gsea_tatstable_tat16h <- gobp_gsea_tatstable_tat16h_simple %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES > 0, p.adjust < 0.05, nGenes >= 20) 

up_gobp_gsea_tatstable_tat16h_simm <- pairwise_termsim(up_gobp_gsea_tatstable_tat16h)
treeplot(up_gobp_gsea_tatstable_tat16h_simm, nWords = 2, nCluster = 5, label_format_cladelab	= 15)

ggsave(file.path(fig_dir, "GOBP_GSEA_UP_Tat.stable_vs_Tat.16h.png"), units = "cm", width = 12, height = 8, dpi = 600, scaling = 0.35)

#### DOWN
down_gobp_gsea_tatstable_tat16h <- gobp_gsea_tatstable_tat16h_simple %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES < 0, p.adjust < 0.05, nGenes >= 20) 

down_gobp_gsea_tatstable_tat16h_simm <- pairwise_termsim(down_gobp_gsea_tatstable_tat16h)
treeplot(down_gobp_gsea_tatstable_tat16h_simm, nWords = 2, nCluster = 5, label_format_cladelab	= 15)

ggsave(file.path(fig_dir, "GOBP_GSEA_DOWN_Tat.stable_vs_Tat.16h.png"), units = "cm", width = 13, height = 6, dpi = 600, scaling = 0.35)


## Tat.16h vs Tat.0h -------------------------------------------------------
### Prepare pre-ranked gene list
gl_tat16h_tat0h <- 
  abs(log10(dgea_table_entrez$Tat.16h_vs_Tat.0h_padj)) * 
  sign(dgea_table_entrez$Tat.16h_vs_Tat.0h_log2FC)
names(gl_tat16h_tat0h) <- as.character(dgea_table_entrez$entrezID)
gl_tat16h_tat0h <- sort(gl_tat16h_tat0h, decreasing = TRUE)

### Run GSEA against GO-BP database
### Use seed for reproducibility
set.seed(42) 
gobp_gsea_tat16h_tat0h <- gseGO(
  geneList = gl_tat16h_tat0h, 
  ont = "BP", 
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Save GSEA results
saveRDS(gobp_gsea_tat16h_tat0h, file.path(output_dir, "gobp_gsea_tat16h_tat0h.rds"))

### Simplify GO results
gobp_gsea_tat16h_tat0h_simple <- gobp_gsea_tat16h_tat0h %>% 
  # Filter GO terms by adjusted pvalue
  filter(p.adjust < 0.05) %>% 
  simplify()

saveRDS(gobp_gsea_tat16h_tat0h_simple, file.path(output_dir, "gobp_gsea_tat16h_tat0h_simple.rds"))

### Save enriched terms
gobp_gsea_tat16h_tat0h_simple@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "gobp_gsea_tat16h_tat0h_simple_UP.padj005.tsv"))

gobp_gsea_tat16h_tat0h_simple@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.05) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "gobp_gsea_tat16h_tat0h_simple_DOWN.padj005.tsv"))

### Treeplot
#### UP
up_gobp_gsea_tat16h_tat0h <- gobp_gsea_tat16h_tat0h_simple %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES > 0, p.adjust < 0.05, nGenes >= 50) 

up_gobp_gsea_tat16h_tat0h_simm <- pairwise_termsim(up_gobp_gsea_tat16h_tat0h)
treeplot(up_gobp_gsea_tat16h_tat0h_simm, nWords = 2, nCluster = 5, label_format_cladelab	= 15)

ggsave(file.path(fig_dir, "GOBP_GSEA_UP_Tat.16h_vs_Tat.0h.png"), units = "cm", width = 12, height = 8, dpi = 600, scaling = 0.35)

#### DOWN
down_gobp_gsea_tat16h_tat0h <- gobp_gsea_tat16h_tat0h_simple %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES < 0, p.adjust < 0.05, nGenes >= 50) 

down_gobp_gsea_tat16h_tat0h_simm <- pairwise_termsim(down_gobp_gsea_tat16h_tat0h)
treeplot(down_gobp_gsea_tat16h_tat0h_simm, nWords = 2, nCluster = 5, label_format_cladelab	= 15)

ggsave(file.path(fig_dir, "GOBP_GSEA_DOWN_Tat.16h_vs_Tat.0h.png"), units = "cm", width = 12, height = 6, dpi = 600, scaling = 0.35)


## Tat.stable vs control -------------------------------------------------------
### Prepare pre-ranked gene list
gl_tatstable_control <- 
  abs(log10(dgea_table_entrez$Tat.stable_vs_control_padj)) * 
  sign(dgea_table_entrez$Tat.stable_vs_control_log2FC)
names(gl_tatstable_control) <- as.character(dgea_table_entrez$entrezID)
gl_tatstable_control <- sort(gl_tatstable_control, decreasing = TRUE)

### Run GSEA against GO-BP database
### Use seed for reproducibility
set.seed(42) 
gobp_gsea_tatstable_control <- gseGO(
  geneList = gl_tatstable_control, 
  ont = "BP", 
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 1, 
  verbose = FALSE,
  seed = TRUE)

### Save GSEA results
saveRDS(gobp_gsea_tatstable_control, file.path(output_dir, "gobp_gsea_tatstable_control.rds"))

### Simplify GO results
gobp_gsea_tatstable_control_simple <- gobp_gsea_tatstable_control %>% 
  # Filter GO terms by adjusted pvalue
  filter(p.adjust < 0.1) %>% 
  simplify()

saveRDS(gobp_gsea_tatstable_control_simple, file.path(output_dir, "gobp_gsea_tatstable_control_simple.rds"))

### Save enriched terms
gobp_gsea_tatstable_control_simple@result %>% 
  dplyr::filter(NES > 0, p.adjust < 0.1) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "gobp_gsea_tatstable_control_simple_UP.padj005.tsv"))

gobp_gsea_tatstable_control_simple@result %>% 
  dplyr::filter(NES < 0, p.adjust < 0.1) %>% 
  dplyr::select(ID, NES, p.adjust, Description) %>% 
  write_tsv(file.path(output_dir, "gobp_gsea_tatstable_control_simple_DOWN.padj005.tsv"))

### Treeplot
#### UP
up_gobp_gsea_tatstable_control <- gobp_gsea_tatstable_control_simple %>% 
  mutate(nGenes = str_count(core_enrichment, "/") + 1) %>% 
  filter(NES > 0, p.adjust < 0.1, nGenes >= 10) 

up_gobp_gsea_tatstable_control_simm <- pairwise_termsim(up_gobp_gsea_tatstable_control)
treeplot(up_gobp_gsea_tatstable_control_simm, nWords = 2, nCluster = 3, label_format_cladelab	= 15)

ggsave(file.path(fig_dir, "GOBP_GSEA_UP_Tat.stable_vs_control.png"), units = "cm", width = 10, height = 4, dpi = 300, scaling = 0.35)

#### DOWN - No enriched categories...
down_gobp_gsea_tatstable_control <- gobp_gsea_tatstable_control_simple %>% 
  filter(NES < 0, p.adjust < 0.1) 
