## Load libraries --------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(KEGGREST)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
output_dir <- "results/tables/02_GSEA/02_KEGG"
fig_dir <- "results/figures/02_GSEA/02_KEGG"

# hsa04120 - Ubiquitin mediated proteolysis

## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


## Get gene list associated with KEGG pathway ----------------------------------
query <- keggGet(c("hsa04120"))

hsa04120_entrez <- query[[1]]$GENE[c(TRUE, FALSE)]

### Save gene table
dgea_table_entrez %>% 
  filter(entrezID %in% hsa04120_entrez) %>% 
  write_tsv(file.path(output_dir, "hsa04120_genes.tsv"))


## Plot gene expression changes by subcategory of hsa04120 ---------------------
### Load gene table with sabcategory annotation --------------------------------
ubi <- read_tsv("results/tables/hsa04120_subcats.txt")
ubi <- ubi %>% 
  separate_rows(subcategory, sep = "; ")

ubi_lfc <- left_join(ubi, dgea_table_entrez, by = c("gene_name", "ensID"))


### Ridgeline plot -------------------------------------------------------------
subcategory_levels <- c(
  "Ubiquitins", 
  "Ubiquitin-activating enzymes (E1)",
  "Ubiquitin-conjugating enzymes (E2)",
  "HECT type E3",
  "U-box type E3",
  "Single Ring-finger type E3",
  "Multi subunit Ring-finger type E3")

comparison_levels <- comparison_order[2:4]

ubi_lfc %>% 
  select(
    gene_name, subcategory, 
    Tat.16h_vs_Tat.0h_log2FC, 
    Tat.stable_vs_control_log2FC, 
    Tat.stable_vs_Tat.16h_log2FC) %>% 
  mutate(subcategory = factor(subcategory, levels = rev(subcategory_levels))) %>% 
  pivot_longer(cols = ends_with("log2FC"), names_to = "comparison", values_to = "log2FC") %>% 
  mutate(
    comparison = str_remove(comparison, "_log2FC"),
    comparison =  factor(comparison, levels = comparison_levels)) %>% 
  ggplot(aes(x = log2FC, y = subcategory, fill = comparison)) +
  geom_density_ridges(
    alpha = 0.7) +
  scale_fill_manual(
    values = cols4[c(4,1,2)],
    labels = comparison_names_format[2:4]) +
  labs(x = expression(log[2]~Fold~Change), y = "", fill = "Comparisons") +
  theme_ridges() +
  theme(
    axis.title.x = element_text(hjust = 0.5),
    legend.text = element_text(size = 12, hjust = 0),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(fig_dir, "hsa04120_ridgeline_3.png"), units = "cm", width = 14, height = 6, scale = 2.2)
ggsave(file.path(fig_dir, "hsa04120_ridgeline_3.pdf"), units = "cm", width = 14, height = 6, scale = 2.2)


#### 2 comparisons -------------------------------------------------------------
ubi_lfc %>% 
  select(
    gene_name, subcategory, 
    Tat.16h_vs_Tat.0h_log2FC, 
    Tat.stable_vs_Tat.16h_log2FC) %>% 
  mutate(subcategory = factor(subcategory, levels = rev(subcategory_levels))) %>% 
  pivot_longer(cols = ends_with("log2FC"), names_to = "comparison", values_to = "log2FC") %>% 
  mutate(
    comparison = str_remove(comparison, "_log2FC"),
    comparison =  factor(comparison, levels = comparison_levels[2:3])) %>% 
  ggplot(aes(x = log2FC, y = subcategory, fill = comparison)) +
  geom_density_ridges(
    alpha = 0.7) +
  scale_fill_manual(
    values = cols4[c(1,2)],
    labels = comparison_names_format[3:4]) +
  labs(x = expression(log[2]~Fold~Change), y = "", fill = "Comparisons") +
  theme_ridges() +
  theme(
    axis.title.x = element_text(hjust = 0.5),
    legend.text = element_text(size = 12, hjust = 0),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(fig_dir, "hsa04120_ridgeline_2.png"), units = "cm", width = 14, height = 6, scale = 2.2)
ggsave(file.path(fig_dir, "hsa04120_ridgeline_2.pdf"), units = "cm", width = 14, height = 6, scale = 2.2)


### Volcano plot with HECT type E3 genes highlighted ---------------------------
ubi_lfc4vln <- ubi_lfc %>% 
  select(
    gene_name, subcategory, 
    Tat.16h_vs_Tat.0h_log2FC,
    Tat.16h_vs_Tat.0h_padj) 

ggplot() +
  geom_point(
    aes(x = Tat.16h_vs_Tat.0h_log2FC, y = -log10(Tat.16h_vs_Tat.0h_padj)),
    filter(ubi_lfc4vln, subcategory != "HECT type E3"),
    color = "grey", alpha = 0.7) +
  geom_point(
    aes(x = Tat.16h_vs_Tat.0h_log2FC, y = -log10(Tat.16h_vs_Tat.0h_padj)),
    filter(ubi_lfc4vln, subcategory == "HECT type E3"),
    color = "red", size = 2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(
      x = Tat.16h_vs_Tat.0h_log2FC, 
      y = -log10(Tat.16h_vs_Tat.0h_padj),
      label = gene_name),
    filter(ubi_lfc4vln, subcategory == "HECT type E3"),
    max.overlaps = 100, size = 3) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    limits = c(-1.23, 1.23)) +
  labs(
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~p.adjust)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    aspect.ratio = 1)

ggsave(file.path(fig_dir, "hsa04120_vln.png"), units = "cm", width = 6, height = 6, scale = 2.2)
ggsave(file.path(fig_dir, "hsa04120_vln.pdf"), units = "cm", width = 6, height = 6, scale = 2.2)

