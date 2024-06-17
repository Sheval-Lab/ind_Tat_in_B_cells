## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
output_dir <- "results/tables/02_GSEA"


## Load DGEA results -----------------------------------------------------------
dgea_table <- read_tsv(file.path(input_dir, "LFC_DGEA_table.tsv"))


## Prepare input data for GSEA -------------------------------------------------
### Gene IDs conversion
geneIDs_table <- bitr(dgea_table$ensID, 
                      fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), 
                      OrgDb="org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID), !duplicated(ENSEMBL)) %>% 
  dplyr::rename(ensID = ENSEMBL, entrezID = ENTREZID, gene_symbol = SYMBOL)

dgea_table_entrez <- left_join(geneIDs_table, dgea_table, by = "ensID")


## Save table with Entrez IDs ---------------------------------------------------
write_tsv(dgea_table_entrez, file.path(output_dir, "LFC_DGEA_table_EntrezID.tsv"))
