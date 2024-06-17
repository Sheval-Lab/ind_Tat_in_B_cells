## Load libraries --------------------------------------------------------------
library(tidyverse)
library(org.Hs.eg.db)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load DGEA results -----------------------------------------------------------
degs_table <- read_tsv(file.path(input_dir, "LFC_DE_table.tsv"))


## Add GO categories -----------------------------------------------------------
keytypes(org.Hs.eg.db)
go <- select(org.Hs.eg.db, keys(org.Hs.eg.db), c("ENSEMBL", "GOALL"))

tf_go <- c(
  "GO:0140223", # general transcription factors
  "GO:0003700", # DNA-binding transcription factors
  "GO:0003712" # transcription co-regulators
  )


### Annotate all DE genes with GO categories
degs_go <- go %>% 
  dplyr::select(ensID = ENSEMBL, GO = GOALL, ONTOLOGY = ONTOLOGYALL) %>% 
  right_join(degs_table, by = "ensID") %>% 
  distinct()

degs_go_tf <- degs_go %>% 
  filter(GO %in% tf_go) %>% 
  dplyr::select(gene_name, comparison, log2FC) %>% 
  distinct()

degs_go_tf %>% 
  mutate(
    direction = if_else(log2FC > 0, "UP", "DOWN")) %>% 
  # count(comparison, direction)
  dplyr::select(gene_name, comparison, log2FC) %>% 
  distinct() %>% 
  ggplot(aes(x = log2FC, fill = comparison)) +
  geom_histogram() +
  facet_wrap(vars(comparison))

