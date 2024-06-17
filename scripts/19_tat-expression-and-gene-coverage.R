## Load libraries --------------------------------------------------------------
library(tidyverse)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "data/tat_coverage"
fig_dir <- "results/figures"

## Tat gene plasmid read counts ------------------------------------------------
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


## Tat gene coverage -----------------------------------------------------------
cov_files <- list.files(input_dir, pattern = ".coverage.txt$", full.names = TRUE)

tat_cov <- map(
  cov_files, 
  read_tsv,
  col_names = c("chr", "position", "coverage"),
  id = "sample") %>% 
  list_rbind() %>% 
  mutate(
    sample = map_chr(str_split(sample, "/"), last),
    sample = str_remove(sample, ".coverage.txt"),
    sample_group = str_remove(sample, ".\\d+$"))


tat_cov %>% 
  ggplot(aes(x = position, y = coverage, fill = sample_group)) +
  geom_col() +
  facet_wrap(vars(sample), ncol = 3, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "tat_coverage.png"), units = "cm", width = 10, height = 8, scale = 1.5)
