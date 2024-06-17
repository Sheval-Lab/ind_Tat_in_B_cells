## Load libraries --------------------------------------------------------------
library(tidyverse)
library(DESeq2)
#library(ggh4x)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load DESeq object -----------------------------------------------------------
dds <- readRDS(file.path(input_dir, "dds.rds"))


## Visualize expression data with PCA ------------------------------------------
pca <- plotPCA(rlog(dds), intgroup = "sample", ntop = 500, returnData = TRUE)
pca_percent_var <- round(100 * attr(pca, "percentVar"))

pca %>% 
  mutate(
    sample = factor(sample, levels = sample_order),
    control = if_else(sample %in% c("control", "Tat.0h"), "control", "experiment")) %>% 
  ggplot(aes(PC1, PC2, color = sample)) + 
  geom_point(size = 4, alpha = 0.8) +
  # ggforce::geom_ellipse(aes(x0 = -2.5, y0 = -2.7, a = 3, b = 0.9, angle =  3.8 * pi / 4), show.legend = FALSE) +
  ggforce::geom_mark_ellipse(
    aes(filter = control == "control", color = control, label = "Control samples"), 
    label.fontface = "plain", label.fontsize = 10, show.legend = FALSE) +
  # geom_segment(
  #   aes(x = -4, y = -1.5, xend = -5, yend = 3.7), arrow = arrow(length = unit(0.5, "cm")), 
  #   color = cols4[1], show.legend = FALSE) +
  geom_segment(
    aes(x = 0, y = -2.5, xend = 9.3, yend = 0.7), arrow = arrow(length = unit(0.5, "cm")), 
    color = cols4[2], show.legend = FALSE) +
  geom_segment(
    aes(x = 9.1, y = 1.3, xend = -4.3, yend = 4.2), arrow = arrow(length = unit(0.5, "cm")), 
    color = cols4[1], show.legend = FALSE) +
  annotate("text", x = 3, y = 2.9, label = "transition from acute to chronic effects", angle = -24, size = 3) +
  annotate("text", x = 4.8, y = -0.5, label = "acute Tat exposure", angle = 34.7, size = 3) +
  labs(x = "acute effects", y = "chronic effects") +
  scale_color_manual(
    name = "",
    values = cols4_sample,
    labels = sample_names_format) +
  expand_limits(y = -3.5) +
  theme_classic() +
  theme(text = element_text(size = 12, color = "black"), 
        axis.title.x = element_text(hjust = 1),
        axis.title.y = element_text(hjust = 1),
        # axis.text = element_blank(),
        # axis.ticks = element_blank(),
        legend.text.align = 0,
        # legend.key = element_rect(color = NA, fill = NA),
        aspect.ratio = 1)

ggsave(file.path(fig_dir, "PCA_effects.png"))
# ggsave(file.path(fig_dir, "PCA_effects.png"), units = "cm", width = 5, height = 5, scale = 2, dpi = 300)
# ggsave(file.path(fig_dir, "PCA_effects.pdf"), units = "cm", width = 5, height = 5, scale = 2)

