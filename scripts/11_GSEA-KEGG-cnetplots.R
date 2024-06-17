## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/02_GSEA"
fig_dir <- "results/figures/02_GSEA/02_KEGG"


## Load DGEA results -----------------------------------------------------------
dgea_table_entrez <- read_tsv(file.path(input_dir, "LFC_DGEA_table_EntrezID.tsv"))


## Tat.stable vs control -------------------------------------------------------
### Load GSEA KEGG results 
kegg_gsea_tatstable_control <- readRDS(file.path(input_dir, "02_KEGG", "kegg_gsea_tatstable_control.rds"))

### Prepare gene list of log2FC values Tat.stable vs control
fc_list <- dgea_table_entrez %>% filter(Tat.stable_vs_control_padj < 0.05) %>% pull(Tat.stable_vs_control_log2FC)
names(fc_list) <- dgea_table_entrez %>% filter(Tat.stable_vs_control_padj < 0.05) %>% pull(gene_symbol) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)

### Keep genes with low p-value in core enrichment
kegg_gsea_tatstable_control_flt <- kegg_gsea_tatstable_control 
kegg_gsea_tatstable_control_flt@result <- kegg_gsea_tatstable_control_flt@result %>% 
  separate_rows(core_enrichment) %>% 
  filter(core_enrichment %in% names(fc_list)) %>% 
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge) %>% 
  summarise(core_enrichment = paste(core_enrichment, collapse = "/")) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  magrittr::set_rownames(.$ID)

### Split GSEA KEGG results into UP- and DOWN-regulated pathways 
#### UP
up_kegg_gsea_tatstable_control <- kegg_gsea_tatstable_control_flt %>% 
  filter(NES > 0,p.adjust < 0.1)

#### DOWN
dn_kegg_gsea_tatstable_control <- kegg_gsea_tatstable_control_flt %>% 
  filter(NES < 0, p.adjust < 0.1)


### Make cnetplot with UP-regulated pathways
up_cnetpl <- cnetplot(
  up_kegg_gsea_tatstable_control, 
  showCategory = nrow(up_kegg_gsea_tatstable_control),
  foldChange = fc_list, #node_label = "none",
  colorEdge = TRUE, 
  cex_category = 1,
  cex_gene = 1,
  cex_label_gene = 0.6,
  cex_label_category = 0.7) +
  scale_color_gradient2(
    name = "log2FoldChange", 
    low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    limits = c(0, 15),
    breaks = seq(5, 15, 5)) +
  ggraph::scale_edge_color_manual(values = c("#6E495A", "#C0B298", "#A4778B", "#AA4586")) +
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4), edge_colour = "none") +
  theme(
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.6)))

ragg::agg_png(file.path(fig_dir, "UP_cnet_Tat.stable_vs_control.png"), units = "cm", width = 7, height = 4, res = 600, scaling = 0.5)
up_cnetpl
dev.off()

ggsave(file.path(fig_dir, "UP_cnet_Tat.stable_vs_control.pdf"), up_cnetpl, units = "cm", width = 7, height = 4, scale = 2)


### Make cnetplot with DOWN-regulated pathways
dn_cnetpl <- cnetplot(
  dn_kegg_gsea_tatstable_control, 
  showCategory = nrow(dn_kegg_gsea_tatstable_control),
  foldChange = fc_list, 
  colorEdge = TRUE, 
  cex_category = 1,
  cex_gene = 1,
  cex_label_gene = 0.6,
  cex_label_category = 0.7) +
  scale_color_gradient2(
    name = "log2FoldChange", 
    low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    limits = c(0, 15),
    breaks = seq(5, 15, 5)) +
  ggraph::scale_edge_color_manual(values = c("#A0CA92", "#75B09C", "#998650")) +
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4), edge_colour = "none") +
  theme(
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.6)))

ragg::agg_png(file.path(fig_dir, "DOWN_cnet_Tat.stable_vs_control.png"), units = "cm", width = 7, height = 4, res = 600, scaling = 0.5)
dn_cnetpl
dev.off()

ggsave(file.path(fig_dir, "DOWN_cnet_Tat.stable_vs_control.pdf"), dn_cnetpl, units = "cm", width = 7, height = 4, scale = 2)


## Tat.16h vs Tat.0h ----------------------------------------------------------
kegg_gsea_tat16h <- readRDS(file.path(input_dir, "02_KEGG", "kegg_gsea_tat16h_tat0h.rds"))

### Prepare gene list of log2FC values Tat.16h vs Tat.0h
fc_list <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05) %>% pull(Tat.16h_vs_Tat.0h_log2FC)
names(fc_list) <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05) %>% pull(gene_symbol) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)

gene_list_flt <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05) %>% pull(gene_symbol)

# fc_list_flt <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05, abs(Tat.16h_vs_Tat.0h_log2FC) >= log2(1.5)) %>% pull(Tat.16h_vs_Tat.0h_log2FC)
# names(fc_list_flt) <- dgea_table_entrez %>% filter(Tat.16h_vs_Tat.0h_padj < 0.05, abs(Tat.16h_vs_Tat.0h_log2FC) >= log2(1.5)) %>% pull(gene_symbol) %>% as.character()
# fc_list_flt <- sort(fc_list_flt, decreasing = TRUE)

### Filter genes with high p-value from core enrichment
kegg_gsea_tat16h_flt <- kegg_gsea_tat16h 
kegg_gsea_tat16h_flt@result <- kegg_gsea_tat16h_flt@result %>% 
  separate_rows(core_enrichment) %>% 
  filter(core_enrichment %in% gene_list_flt) %>% 
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge) %>% 
  summarise(core_enrichment = paste(core_enrichment, collapse = "/")) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  magrittr::set_rownames(.$ID)
  

### Split GSEA KEGG results into UP- and DOWN-regulated pathways 
#### UP
up_kegg_gsea_tat16h <- kegg_gsea_tat16h_flt %>% 
  filter(NES > 0, p.adjust < 0.05) %>% 
  # Select only top 20 pathways
  arrange(p.adjust) %>% 
  filter(row_number(p.adjust) <= 20)

#### DOWN
dn_kegg_gsea_tat16h <- kegg_gsea_tat16h_flt %>% 
  filter(NES < 0, p.adjust < 0.05)


### Make cnetplot with UP-regulated pathways
up_cnetpl <- cnetplot(up_kegg_gsea_tat16h, 
                      showCategory = nrow(up_kegg_gsea_tat16h),
                      foldChange = fc_list, #node_label = "none",
                      colorEdge = TRUE,
                      cex_category = 0.7,
                      cex_label_gene = 0.5,
                      cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") +
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(file.path(fig_dir, "top20_UP_cnet_Tat.16h_vs_Tat.0h.png"), up_cnetpl, units = "cm", width = 8, height = 5, scale = 1.5)
ggsave(file.path(fig_dir, "top20_UP_cnet_Tat.16h_vs_Tat.0h.pdf"), up_cnetpl, units = "cm", width = 8, height = 5, scale = 1.5)


### Make cnetplot with DOWN-regulated pathways 

dn_kegg_gsea_tat16h@result %>% dplyr::select(ID, Description)
up_kegg_gsea_tat16h@result %>% dplyr::select(ID, Description)


### Make cnetplot by pathway groups
make_cnetplot <- function(kegg_gsea, kegg_ids, fc_list, colors, max_size) {
  cnetplot(
    filter(kegg_gsea, ID %in% kegg_ids), 
    showCategory = nrow(kegg_gsea),
    foldChange = fc_list, 
    colorEdge = TRUE, 
    # node_label = "gene",
    cex_category = 1,
    cex_gene = 5,
    cex_label_gene = 0.5,
    cex_label_category = 0.8) +
    scale_color_gradient2(
      name = "log2FoldChange", 
      low = "blue", mid = "white", high = "red",
      limits = c(-2, 3.6),
      breaks = seq(-2, 3, 1)) +
    scale_size_continuous(
      limits = c(0, max_size),
      breaks = seq(10, 30, 10)) +
    ggraph::scale_edge_color_manual(values = colors) +
    guides(color = guide_colourbar(barwidth = 0.7, barheight = 4), edge_colour = "none") +
    theme(
      legend.title = element_text(size = rel(0.7)),
      legend.text = element_text(size = rel(0.6)))
}

#### Replication
replication_ids <- c("hsa03030")

cnet_replication <- make_cnetplot(kegg_gsea_tat16h_flt, replication_ids, fc_list, c("#903953"), 30)

ragg::agg_png(file.path(fig_dir, "Tat.16h_vs_Tat.0h_replication.png"), units = "cm", width = 4, height = 3, res = 600, scaling = 0.5)
cnet_replication
dev.off()

ggsave(file.path(fig_dir, "Tat.16h_vs_Tat.0h_replication.pdf"), cnet_replication, units = "cm", width = 4, height = 3, scale = 2)


#### Transcription & mRNA processing
trancription_ids <- c("hsa03020", "hsa03040")

cnet_transcription <- make_cnetplot(kegg_gsea_tat16h_flt, trancription_ids, fc_list, c("#0a9396", "#ca6702"), 64)

ragg::agg_png(file.path(fig_dir, "Tat.16h_vs_Tat.0h_transcription.png"), units = "cm", width = 5, height = 4, res = 750, scaling = 0.4)
cnet_transcription
dev.off()

ggsave(file.path(fig_dir, "Tat.16h_vs_Tat.0h_transcription.pdf"), cnet_transcription, units = "cm", width = 5, height = 4, scale = 2.5)


#### Translation
translation_ids <- c("hsa00970", "hsa03008", "hsa03010", "hsa03013")

cnet_translation <- make_cnetplot(kegg_gsea_tat16h_flt, translation_ids, fc_list, c("#94d2bd", "#005f73", "#bb3e03", "#b56576"), 101)

ragg::agg_png(file.path(fig_dir, "Tat.16h_vs_Tat.0h_translation.png"), units = "cm", width = 5, height = 4, res = 1000, scaling = 0.25)
cnet_translation
dev.off()

ggsave(file.path(fig_dir, "Tat.16h_vs_Tat.0h_translation.pdf"), cnet_translation, units = "cm", width = 5, height = 4, scale = 4)


#### Antipathogenic effects
antipathogen_ids <- c("hsa04612", "hsa04672", "hsa04062", "hsa04666", "hsa05131", "hsa04625", "hsa05135", "hsa05146", "hsa04613", "hsa05165")

cnet_antipathogen <- make_cnetplot(kegg_gsea_tat16h_flt, antipathogen_ids, fc_list, c("#C22D2D", "#D5AE4D", "#CDEBA9", "#A7D6A3", "#3C978D", "#b46c10", "#6897D9", "#6C4A8F", "#5A4082", "#E89EB1", "#EAB2A6"), 60)

ragg::agg_png(file.path(fig_dir, "Tat.16h_vs_Tat.0h_antipathogen.png"), units = "cm", width = 6, height = 4, res = 1000, scaling = 0.25)
cnet_antipathogen
dev.off()

ggsave(file.path(fig_dir, "Tat.16h_vs_Tat.0h_antipathogen.pdf"), cnet_antipathogen, units = "cm", width = 6, height = 4, scale = 4)


#### Cytoskeleton, contacts, movement
cytoskeleton_ids <- c("hsa04015", "hsa04014", "hsa04810", "hsa04072", "hsa04670", "hsa04530", "hsa04510", "hsa04540") # hsa04666 (Fc - cytoskeleton rearrangements)

cnet_cytoskeleton <- make_cnetplot(kegg_gsea_tat16h_flt, cytoskeleton_ids, fc_list, c("#2a1a1f", "#023c40", "#009fb7", "#c9e4ca", "#f4f1de", "#ceae69", "#b46c10", "#6e8898"), 33)

ragg::agg_png(file.path(fig_dir, "Tat.16h_vs_Tat.0h_cytoskeleton.png"), units = "cm", width = 5, height = 4, res = 1000, scaling = 0.25)
cnet_cytoskeleton
dev.off()


ggsave(file.path(fig_dir, "Tat.16h_vs_Tat.0h_cytoskeleton.pdf"), cnet_cytoskeleton, units = "cm", width = 5, height = 4, scale = 4)


#### Signal
# hsa04350 - TGF-beta signaling pathway
# hsa04611 - Platelet activation
other_ids <- c("hsa04350", "hsa04611")
