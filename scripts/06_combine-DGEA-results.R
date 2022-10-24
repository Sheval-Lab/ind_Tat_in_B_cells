## Load libraries --------------------------------------------------------------
library(tidyverse)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA/deseq"
output_dir <- "results/tables/01_DGEA"


## Load DGEA results - DEGs ----------------------------------------------------
degs_table_files <- list.files(input_dir, pattern = "LFC.DE.tsv$", full.names = TRUE)

### Combine DEGs into one table
degs_table <- map_dfr(degs_table_files, read_tsv, id = "file")
degs_table <- degs_table %>% 
  # Keep comparison ID
  mutate(comparison = str_split(file, "/") %>% map_chr(last) %>% str_remove(".LFC.DE.tsv$")) %>% 
  select(comparison, everything()) %>% 
  select(-file)

### Save combines DEGs table
write_tsv(degs_table, paste(output_dir, "LFC_DE_table.tsv", sep = "/"))


## Load DGEA results - all genes -----------------------------------------------
dgea_table_files <- list.files(input_dir, pattern = "LFC.tsv$", full.names = TRUE)

dgea_table_read <- function(file){
  # Extract comparison name from file path
  comparison_name = str_split(file, "/") %>% 
    map_chr(last) %>% 
    str_remove(".LFC.tsv$")
  # Read in dataframe
  df = read_tsv(file)
  df %>% 
    select(geneID:log2FC) %>% 
    rename_with(~ str_c(comparison_name, .x, sep = "_"), padj:log2FC)
}

### Combine DGEA results into one table
dgea_table <- map(dgea_table_files, dgea_table_read) %>% purrr::reduce(full_join)

### Save combines DGEA results table
write_tsv(dgea_table, paste(output_dir, "LFC_DGEA_table.tsv", sep = "/"))

