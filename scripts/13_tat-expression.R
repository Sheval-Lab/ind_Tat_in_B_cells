## Load libraries --------------------------------------------------------------
library(tidyverse)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")


## Load Tat gene plasmid read counts -------------------------------------------
tat <- read_tsv("data/plasmid_counts.tab", col_names = c("sample", "count"))

size_factors <- dget("results/tables/size_factors.txt")

tat_norm <- tat %>% 
  arrange(factor(sample, levels = names(size_factors))) %>% 
  mutate(
    size_fct = size_factors,
    ncount = count / size_fct,
    sample_group = str_remove(sample, ".\\d$"))
  

tat_norm %>% 
  group_by(sample_group) %>% 
  summarise(avg_ncount = mean(ncount))

