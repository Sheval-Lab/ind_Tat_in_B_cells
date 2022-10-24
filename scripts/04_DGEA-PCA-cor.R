## Load libraries --------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load DESeq object -----------------------------------------------------------
dds <- readRDS(paste(input_dir, "dds.rds", sep = "/"))


## Visualize expression data with PCA ------------------------------------------
pca <- plotPCA(rlog(dds), intgroup = "sample", ntop = 500, returnData = TRUE)
pca_percent_var <- round(100 * attr(pca, "percentVar"))


pca_pl <- ggplot(pca, aes(PC1, PC2, color = factor(sample, levels = sample_order))) + 
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", pca_percent_var[1], "% variance")) + 
  ylab(paste0("PC2: ", pca_percent_var[2], "% variance")) +
  scale_color_manual(
    name = "",
    values = cols4_sample,
    labels = sample_names_format) +
  theme(text = element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text.align = 0,
        aspect.ratio = 1)

ggsave(paste(fig_dir, "PCA.png", sep = "/"), pca_pl, units = "cm", width = 5, height = 5, scale = 2, dpi = 600)
ggsave(paste(fig_dir, "PCA.pdf", sep = "/"), pca_pl, units = "cm", width = 5, height = 5, scale = 2)


## Visualize cross-sample correlation with heatmap -----------------------------
### Calculate pairwise correlation
sample_cor <- cor(counts(dds, normalized = TRUE),method = "spearman" )

### Column annotation - samples by color
ha <- HeatmapAnnotation(
  sample = gsub(".\\d+$", "", colnames(sample_cor)),
  col = list(sample = cols4_sample),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    sample = list(
      title = "Sample", 
      at = sample_order,
      labels = gt_render(c(
        "RPMI^Tat",
        "RPMI^Tat-ind Doxy+", 
        "RPMI^Tat-ind Doxy-",
        "RPMI")))),
  gp = gpar(col = "white"))

### Row annotation - same as column annotation
ra <- rowAnnotation(
  sample = gsub(".\\d+$", "", colnames(sample_cor)),
  col = list(sample = cols4_sample),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_legend_param = list(
    sample = list(
      title = "", 
      at = sample_order,
      labels = gt_render(c(
        "RPMI^Tat",
        "RPMI^Tat-ind Doxy+", 
        "RPMI^Tat-ind Doxy-",
        "RPMI")))),
  gp = gpar(col = "white"))

### Heatmap body
hm <- Heatmap(
  sample_cor, 
  col = colorRampPalette(c("#F0FADA", "#E2737C", "#DB3844"), bias = 0.5)(50),
  heatmap_legend_param = list(
    title = "Spearman correlation",
    legend_direction = "horizontal",
    legend_height = unit(5, "cm")),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3g", sample_cor[i, j]), x, y, gp = gpar(fontsize = 8))},
  column_labels = paste("rep.", str_extract(colnames(sample_cor), "\\d+$")),
  row_labels = paste("rep.", str_extract(colnames(sample_cor), "\\d+$")),
  top_annotation = ha,
  left_annotation = ra)

### Save heatmap
ragg::agg_png(paste(fig_dir, "sample_cor_hm.png", sep = "/"), units = "cm", width = 15, height = 10, res = 300, scaling = 3/4)
# pdf(paste(fig_dir, "sample_cor_hm.pdf", sep = "/"), height = 5)
draw(
  hm,
  heatmap_legend_side = "top",
  annotation_legend_side = "right")
dev.off()


