## Load libraries --------------------------------------------------------------
library(tidyverse)

input_dir <- "data/counts"
output_dir <- "results/tables"


## Load counts data ------------------------------------------------------------
counts_files <- list.files(input_dir, full.names = TRUE, pattern = ".tab$")

### Extract sample names from file names
sample_names <- counts_files %>% 
  str_split("/") %>% 
  map_chr(last) %>% 
  str_remove(".tab$")

### Read multiple htseq count files into one dataframe
counts <- map_dfc(counts_files, read_tsv, col_names = FALSE) %>% 
  select(1, where(is.numeric)) %>% 
  magrittr::set_colnames(c("geneID", sample_names)) 


## Process gene count table ----------------------------------------------------
### Filter out summary rows
counts <- counts %>% 
  filter(!str_detect(geneID, "^__"))

### Count total counts and save as summary
unique_counts_summary <- counts %>% 
  select(where(is.numeric)) %>% 
  colSums() %>% as_tibble() %>% 
  mutate(sample = sample_names, .before = 1) %>% 
  dplyr::rename(unique = value)


## Save gene count tables ------------------------------------------------------
### Save combined raw counts data
write_tsv(counts, file.path(output_dir, "counts_raw.tsv"))

### Save total counts summary
write_tsv(unique_counts_summary, file.path(output_dir, "unique_counts_summary.tsv"))

