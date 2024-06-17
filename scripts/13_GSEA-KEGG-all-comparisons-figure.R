## Load libraries --------------------------------------------------------------
library(tidyverse)
library(ggh4x)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA/02_KEGG"
fig_dir <- "results/figures/02_GSEA/02_KEGG"


## Load GSEA KEGG results ------------------------------------------------------
gsea_kegg_res_files <- list.files(input_dir, pattern = ".padj01.tsv$", full.names = TRUE)

read_gsea_kegg_res <- function(file) {
  file_name = str_split(file, "/") %>% 
    map_chr(last) %>% 
    str_remove(".padj01.tsv$")
  direction = file_name %>% str_extract("[A-Z]+$")
  comparison = file_name %>% str_remove("kegg_gsea_") %>% str_remove(paste0("_", direction))
  
  read_tsv(file) %>% 
    mutate(
      comparison = comparison,
      direction = direction)
}

gsea_kegg_res <- map_dfr(gsea_kegg_res_files, read_gsea_kegg_res)

### Add Count as number of genes in core enrichment
### And GeneRatio as (number of genes) / setSize
gsea_kegg_res <- gsea_kegg_res %>% 
  mutate(
    Count = str_count(core_enrichment, "/") + 1,
    GeneRatio = Count / setSize)


## Plots -----------------------------------------------------------------------
theme_dpl <- theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    strip.background = element_rect(fill = "transparent", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank())


### All 3 comparisons ----------------------------------------------------------
### Transform data into wide format to plot lollipop graph 
gsea_kegg_res_wd <- gsea_kegg_res %>% 
  # filter(p.adjust < 0.05) %>% # Filter KEGG pathways by p-value
  select(ID, Description, comparison, NES) %>% 
  pivot_wider(names_from = comparison, values_from = NES) %>% 
  mutate(group = case_when(
    (tat16h_tat0h > 0) & (tatstable_tat16h < 0) ~ "Tat up, cell down",
    (tat16h_tat0h < 0) & (tatstable_tat16h > 0) ~ "Tat down, cell up",
    is.na(tat16h_tat0h) & !is.na(tatstable_tat16h) ~ "cell",
    !is.na(tat16h_tat0h) & is.na(tatstable_tat16h) ~ "Tat"
  )) %>% 
  replace_na(list(
    tatstable_tat16h = 0, 
    tat16h_tat0h = 0,
    tatstable_control = 0))



ggplot(gsea_kegg_res_wd, aes(y = Description)) +
  geom_point(aes(x = tatstable_tat16h), color = cols4[1], alpha = 0.5) +
  geom_point(aes(x = tat16h_tat0h), color = cols4[2], alpha = 0.5) +
  geom_point(aes(x = tatstable_control), color = cols4[4], alpha = 0.5) +
  facet_wrap(vars(group), scales = "free") +
  theme_dpl 


### Add column for splitting into facets
gsea_kegg_tatstable_tat16h_gr <- gsea_kegg_tatstable_tat16h_wd %>% 
  mutate(group = case_when(
    (NES_tat16h_tat0h > 0) & (NES_tatstable_tat16h < 0) ~ "Tat up, cell down",
    (NES_tat16h_tat0h < 0) & (NES_tatstable_tat16h > 0) ~ "Tat down, cell up",
    is.na(NES_tat16h_tat0h) & !is.na(NES_tatstable_tat16h) ~ "cell",
    !is.na(NES_tat16h_tat0h) & is.na(NES_tatstable_tat16h) ~ "Tat"
  ))

### Plot
gsea_kegg_tatstable_tat16h_gr %>% 
  ggplot(aes(y = Description)) +
  geom_segment(aes(x = NES_tatstable_tat16h, xend = NES_tat16h_tat0h, yend = Description)) +
  geom_point(aes(x = NES_tatstable_tat16h, size = GeneRatio_tatstable_tat16h, color = "tatstable_tat16h")) +
  geom_point(aes(x = NES_tat16h_tat0h, size = GeneRatio_tat16h_tat0h, color = "tat16h_tat0h")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "DOWN"), alpha = 0.005) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "UP"), alpha = 0.005) +
  scale_fill_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  scale_color_manual(
    values = c("tat16h_tat0h" = cols4[2], "tatstable_tat16h" = cols4[1]),
    labels = comparison_names_format[2:3]) +
  guides(
    fill = guide_none()) +
  labs(x = "", y = "", size = "GeneRatio", color = "Comparison") +
  facet_wrap(
    ~ factor(group, levels = c("Tat down, cell up", "Tat", "Tat up, cell down", "cell")), 
    scales = "free_y") +
  force_panelsizes(rows = c(0.3, 1)) +
  theme_dpl +
  theme(
    legend.text = element_text(size = 12, hjust = 0),
    legend.title = element_text(size = 12))

ggsave(file.path(fig_dir, "ALL_Tat.stable_vs_Tat.16h.png"), units = "cm", width = 20, height = 15, scale = 2)
ggsave(file.path(fig_dir, "ALL_Tat.stable_vs_Tat.16h.pdf"), units = "cm", width = 20, height = 15, scale = 2)

