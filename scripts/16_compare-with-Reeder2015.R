## Load libraries --------------------------------------------------------------
library(tidyverse)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load DGEA results -----------------------------------------------------------
### Our results
degs_table <- read_tsv(file.path(input_dir, "LFC_DE_table.tsv"))
dgea_table <- read_tsv(file.path(input_dir, "LFC_DGEA_table.tsv"))

### Results from Reeder et al., 2015 (https://doi.org/10.7554/eLife.08955)
reeder <- read_csv("data/reeder_2015_from_GEO.gz")
# reeder <- read_csv("data/GSE65688_gene_exp.diff_GFP_TAT_Gene.csv") # unziped file
reeder <- reeder %>% 
  select(
    gene_name = gene, 
    q_value, 
    log2FC_reeder = log2.fold_change., 
    significant) %>% 
  # Remove too large log2FC values
  filter(abs(log2FC_reeder) < 1000)


## Compare our DGEA results with Reeder's --------------------------------------
### Combine tables (all genes)
dgea_combined <- inner_join(dgea_table, reeder, by = "gene_name") %>% 
  mutate(
    status = case_when(
      (Tat.16h_vs_Tat.0h_padj < 0.05) & (Tat.16h_vs_Tat.0h_log2FC >= log2(1.5)) & (q_value < 0.05) & (log2FC_reeder > 0) ~ "DE, signs match",
      (Tat.16h_vs_Tat.0h_padj < 0.05) & (Tat.16h_vs_Tat.0h_log2FC <= -log2(1.5)) & (q_value < 0.05) & (log2FC_reeder < 0) ~ "DE, signs match",
      (Tat.16h_vs_Tat.0h_padj < 0.05) & (Tat.16h_vs_Tat.0h_log2FC >= log2(1.5)) & (q_value < 0.05) & (log2FC_reeder < 0) ~ "DE, signs differ",
      (Tat.16h_vs_Tat.0h_padj < 0.05) & (Tat.16h_vs_Tat.0h_log2FC <= -log2(1.5)) & (q_value < 0.05) & (log2FC_reeder > 0) ~ "DE, signs differ",
      (Tat.16h_vs_Tat.0h_padj < 0.05) & (abs(Tat.16h_vs_Tat.0h_log2FC) >= log2(1.5)) & (q_value > 0.05) ~ "not DE in one experiment",
      (Tat.16h_vs_Tat.0h_padj > 0.05) & (q_value < 0.05) ~ "not DE in one experiment",
      TRUE ~ "not DE in one experiment"))

### Limits - X and Y axis
max_logFC <- max(abs(dgea_combined$Tat.16h_vs_Tat.0h_log2FC), na.rm = TRUE)
max_log2FC_reeder <- max(abs(dgea_combined$log2FC_reeder), na.rm = TRUE)

cols3 <- c("#006C67", "#F194B4", "#FFEBC6")

dgea_cor_pl <- dgea_combined %>% 
  ggplot(aes(x = Tat.16h_vs_Tat.0h_log2FC, y = log2FC_reeder, color = status)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = cols3) +
  coord_cartesian(
    xlim = c(-max_logFC, max_logFC),
    ylim = c(-max_log2FC_reeder, max_log2FC_reeder)) +
  labs(x = "log2FoldChange", y = "log2FoldChange (Reeder 2015)", color = "") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black", size = 10),
    legend.text = element_text(size = 12),
    legend.position = c(0.25, 0.87))


ragg::agg_png(file.path(fig_dir, "DGEA_comparison_with_Reeder2015.png"), units = "cm", width = 10, height = 8, res = 400, scaling = 0.65)
dgea_cor_pl
dev.off()

ggsave(file.path(fig_dir, "DGEA_comparison_with_Reeder2015.pdf"), dgea_cor_pl, units = "cm", width = 10, height = 8, scale = 3/2)


### Combine tables (only DEGs)
degs_combined <- full_join(
  filter(degs_table, comparison == "Tat.16h_vs_Tat.0h"), # Acute effect of Tat
  filter(reeder, q_value < 0.05), # DE genes from Reeder 2015
  by = "gene_name") 

degs_upset_pl <- degs_combined %>% 
  mutate(
    `Upregulated DEGs` = !is.na(log2FC) & log2FC > 0,
    `Downregulated DEGs` = !is.na(log2FC) & log2FC < 0,
    `Upregulated DEGs (Reeder 2015)` = !is.na(log2FC_reeder) & log2FC_reeder > 0,
    `Downregulated DEGs (Reeder 2015)` = !is.na(log2FC_reeder) & log2FC_reeder < 0) %>% 
  select(gene_name, `Upregulated DEGs`, `Downregulated DEGs`, `Upregulated DEGs (Reeder 2015)`, `Downregulated DEGs (Reeder 2015)`) %>% 
  pivot_longer(cols = -gene_name) %>% 
  distinct() %>%
  filter(value) %>% 
  group_by(gene_name) %>% 
  summarise(results = list(name)) %>% 
  ggplot(aes(x = results)) +
  geom_bar(fill = c(rep(cols3[3], 4), rep(cols3[2], 2), rep(cols3[1], 2))) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3) +
  ggupset::scale_x_upset(sets = c("Upregulated DEGs", "Downregulated DEGs", "Upregulated DEGs (Reeder 2015)", "Downregulated DEGs (Reeder 2015)")) +
  scale_y_continuous(breaks = NULL, name = "") +
  labs(x = "") +
  theme(
    panel.background = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = rel(1)))
  
  
ragg::agg_png(file.path(fig_dir, "DEGs_comparison_with_Reeder2015.png"), units = "cm", width = 10, height = 8, res = 400, scaling = 0.65)
degs_upset_pl
dev.off()
  

ggsave(file.path(fig_dir, "DEGs_comparison_with_Reeder2015.pdf"), degs_upset_pl, units = "cm", width = 10, height = 8, scale = 3/2)

