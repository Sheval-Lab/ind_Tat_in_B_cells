## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
fig_dir <- "results/figures/02_GSEA/02_KEGG"


## Load GSEA KEGG results ------------------------------------------------------
kegg_gsea_tatstable_control <- readRDS(file.path(input_dir, "02_KEGG", "kegg_gsea_tatstable_control.rds"))


## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


### Prepare gene list of log2FC values Tat.stable vs control
fc_list <- dgea_table_entrez %>% filter(Tat.stable_vs_control_padj < 0.05) %>% pull(Tat.stable_vs_control_log2FC)
names(fc_list) <- dgea_table_entrez %>% filter(Tat.stable_vs_control_padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)


## Split GSEA KEGG results into UP- and DOWN-regulated pathways ----------------
### UP
up_kegg_gsea_tatstable_control <- kegg_gsea_tatstable_control %>% 
  filter(NES > 0,p.adjust < 0.1)

### DOWN
dn_kegg_gsea_tatstable_control <- kegg_gsea_tatstable_control %>% 
  filter(NES < 0, p.adjust < 0.1)


## Make cnetplot with UP-regulated pathways ------------------------------------
up_cnetpl <- cnetplot(up_kegg_gsea_tatstable_control, 
         showCategory = nrow(up_kegg_gsea_tatstable_control),
         foldChange = fc_list, #node_label = "none",
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") +
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(file.path(fig_dir, "up_KEGG_GSEA_Tat.stable_vs_control.png"), up_cnetpl, units = "cm", width = 8, height = 5, scale = 1.5)
ggsave(file.path(fig_dir, "up_KEGG_GSEA_Tat.stable_vs_control.pdf"), up_cnetpl, units = "cm", width = 8, height = 5, scale = 1.5)


## Make cnetplot with DOWN-regulated pathways ----------------------------------
dn_cnetpl <- cnetplot(dn_kegg_gsea_tatstable_control, 
         showCategory = nrow(dn_kegg_gsea_tatstable_control),
         foldChange = fc_list, #node_label = "none",
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") +
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(file.path(fig_dir, "down_KEGG_GSEA_Tat.stable_vs_control.png"), dn_cnetpl, units = "cm", width = 8, height = 5, scale = 1.5)
ggsave(file.path(fig_dir, "down_KEGG_GSEA_Tat.stable_vs_control.pdf"), dn_cnetpl, units = "cm", width = 8, height = 5, scale = 1.5)
