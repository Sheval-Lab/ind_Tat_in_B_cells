## Load libraries --------------------------------------------------------------
library(tidyverse)
library(colorspace)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load DGEA results -----------------------------------------------------------
degs_table <- read_tsv(file.path(input_dir, "LFC_DE_table.tsv"))


## Calculate DEGs numbers ------------------------------------------------------
### Calculate number of all types of DEGs
degs_number_all <- degs_table %>% 
  mutate(direction = if_else(log2FC > 0, "UP", "DOWN")) %>% 
  dplyr::count(comparison, direction) %>% 
  mutate(type = "ALL") %>% 
  ungroup()


### Calculate number of protein-coding DEGs
degs_number_pc <- degs_table %>% 
  filter(gene_type == "protein_coding") %>% 
  mutate(direction = if_else(log2FC > 0, "UP", "DOWN")) %>% 
  dplyr::count(comparison, direction) %>% 
  mutate(type = "protein_coding") %>% 
  ungroup()


### Combine dataframes with DEGs numbers
degs_number <- bind_rows(degs_number_all, degs_number_pc) %>% 
  mutate(
    n = if_else(direction == "DOWN", -n, n),
    text_position = ifelse(direction == "DOWN", -270, 320),
    comparison = factor(comparison, levels = comparison_order))


## Plot DEGs number as barplot -------------------------------------------------
facet_names <- c(ALL = "All genes", protein_coding = "Protein-coding genes")
y_axis <- seq(-500, 1000, 500)

degs_pl <- degs_number %>% 
  filter(comparison != "Tat.stable_vs_Tat.16h") %>% 
  ggplot(aes(x = comparison, y = n)) + 
  geom_col(
    aes(color = direction, fill = after_scale(lighten(color, 0.3))),
    position = "identity") + 
  labs(x = "", y = "# DEGs", color = "") +
  geom_text(aes(x = comparison, y = text_position, label = abs(n)), size = 3) +
  facet_wrap(~ type, labeller = labeller(type = facet_names), scales = "free_x") +
  # Edit color scheme
  scale_color_manual(
    breaks = c("UP", "DOWN"),
    values = c("UP" = "red", "DOWN" = "blue"),
    labels = c("Upregulated", "Downregulated")) +
  # Edit sample names
  scale_x_discrete(labels = rep(comparison_names_format, times = 2)) +
  # Edit Y-axis
  scale_y_continuous(breaks = y_axis, labels = abs(y_axis)) +
  # Edit theme
  theme(text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white", colour = "black", size = 0.5, linetype = "solid"),
        strip.text = element_text(face="bold", size = 10), 
        legend.position = "top",
        aspect.ratio = 1)



ragg::agg_png(file.path(fig_dir, "DEGs_number_barplot.png"), units = "cm", width = 10, height = 8, res = 400, scaling = 0.75)
degs_pl
dev.off()

ggsave(file.path(fig_dir, "DEGs_number_barplot.pdf"), degs_pl, units = "cm", width = 10, height = 8, scale = 4/3)
