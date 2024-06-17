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

### Tat.stable vs control ------------------------------------------------------
gsea_kegg_tatstable <- gsea_kegg_res %>% 
  filter(comparison == "tatstable_control")

size_limits <- c(min(gsea_kegg_tatstable$Count), max(gsea_kegg_tatstable$Count))
size_breaks <- seq(size_limits[1], size_limits[2], 3)

color_limits <- c(round(min(gsea_kegg_tatstable$p.adjust), 3) - 0.001, round(max(gsea_kegg_tatstable$p.adjust), 3) + 0.001)
color_breaks <- seq(0.07, 0.09, 0.01)

### UP
dpl_gsea_kegg_tatstable_up <- gsea_kegg_tatstable %>% 
  filter(direction == "UP") %>% 
  ggplot(aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(color = p.adjust, size = Count)) +
  facet_wrap(~ direction, labeller = as_labeller(c("UP" = "Activated\nKEGG pathways"))) +
  labs(y = "") +
  scale_x_continuous(breaks = c(0.09, 0.11, 0.13), expand = expansion(add = 0.005)) +
  scale_size_continuous(breaks = size_breaks, limits = size_limits) +
  scale_color_gradient(
    low = "red", high = "blue",
    breaks = color_breaks, limits = color_limits) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(reverse = TRUE, order = 2)) +
  theme_dpl
  
ragg::agg_png(file.path(fig_dir, "UP_Tat.stable_vs_control.png"), units = "cm", width = 7, height = 4, res = 600, scaling = 0.5)
dpl_gsea_kegg_tatstable_up
dev.off()

ggsave(file.path(fig_dir, "UP_Tat.stable_vs_control.pdf"), dpl_gsea_kegg_tatstable_up, units = "cm", width = 7, height = 4, scale = 2)

### DOWN
dpl_gsea_kegg_tatstable_dn <- gsea_kegg_tatstable %>% 
  filter(direction == "DOWN") %>% 
  ggplot(aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(color = p.adjust, size = Count)) +
  facet_wrap(~ direction, labeller = as_labeller(c("DOWN" = "Suppressed\nKEGG pathways"))) +
  labs(y = "") +
  scale_x_continuous(breaks = c(0.11, 0.13, 0.15), expand = expansion(add = 0.005)) +
  scale_size_continuous(breaks = size_breaks, limits = size_limits) +
  scale_color_gradient(
    low = "red", high = "blue",
    breaks = color_breaks, limits = color_limits) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(reverse = TRUE, order = 2)) +
  theme_dpl

ragg::agg_png(file.path(fig_dir, "DOWN_Tat.stable_vs_control.png"), units = "cm", width = 7.5, height = 4, res = 600, scaling = 0.5)
dpl_gsea_kegg_tatstable_dn
dev.off()

ggsave(file.path(fig_dir, "DOWN_Tat.stable_vs_control.pdf"), dpl_gsea_kegg_tatstable_dn, units = "cm", width = 7.5, height = 4, scale = 2)


### Tat.16h vs Tat.0h ----------------------------------------------------------
gsea_kegg_tat16h <- gsea_kegg_res %>% 
  filter(comparison == "tat16h_tat0h") %>% 
  # Filter KEGG pathways by p-value
  filter(p.adjust < 0.05)

size_limits <- c(min(gsea_kegg_tat16h$Count), max(gsea_kegg_tat16h$Count))
size_breaks <- seq(
  plyr::round_any(size_limits[1], 10, ceiling), 
  plyr::round_any(size_limits[2], 10, floor), 30)

color_limits <- c(round(min(gsea_kegg_tat16h$p.adjust), 3) , round(max(gsea_kegg_tat16h$p.adjust), 3) + 0.006)
color_breaks <- seq(0.01, 0.05, 0.01)

### UP - all
dpl_gsea_kegg_tat16h_up <- gsea_kegg_tat16h %>% 
  filter(direction == "UP") %>% 
  ggplot(aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(color = p.adjust, size = Count)) +
  facet_wrap(~ direction, labeller = as_labeller(c("UP" = "Activated\nKEGG pathways"))) +
  labs(y = "") +
  # scale_x_continuous(breaks = c(0.09, 0.11, 0.13), expand = expansion(add = 0.005)) +
  scale_size_continuous(breaks = size_breaks, limits = size_limits) +
  scale_color_gradient(
    low = "red", high = "blue",
    breaks = color_breaks, limits = color_limits) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(reverse = TRUE, order = 2)) +
  theme_dpl

ggsave(file.path(fig_dir, "UP_Tat.16h_vs_Tat.0h.png"), dpl_gsea_kegg_tat16h_up, units = "cm", width = 7, height = 7, scale = 2.5)
ggsave(file.path(fig_dir, "UP_Tat.16h_vs_Tat.0h.pdf"), dpl_gsea_kegg_tat16h_up, units = "cm", width = 7, height = 7, scale = 2.5)


### UP - top20
dpl_gsea_kegg_tat16h_up_top20 <- gsea_kegg_tat16h %>% 
  filter(direction == "UP") %>% 
  slice_max(order_by = Count, n = 20) %>% 
  ggplot(aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(color = p.adjust, size = Count)) +
  facet_wrap(~ direction, labeller = as_labeller(c("UP" = "Activated\nKEGG pathways"))) +
  labs(y = "") +
  # scale_x_continuous(breaks = c(0.09, 0.11, 0.13), expand = expansion(add = 0.005)) +
  scale_size_continuous(breaks = size_breaks, limits = size_limits) +
  scale_color_gradient(
    low = "red", high = "blue",
    breaks = color_breaks, limits = color_limits) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(reverse = TRUE, order = 2)) +
  theme_dpl

ragg::agg_png(file.path(fig_dir, "top20count_UP_Tat.16h_vs_Tat.0h.png"), units = "cm", width = 9, height = 6, res = 600, scaling = 0.5)
dpl_gsea_kegg_tat16h_up_top20
dev.off()

ggsave(file.path(fig_dir, "top20count_UP_Tat.16h_vs_Tat.0h.pdf"), dpl_gsea_kegg_tat16h_up_top20, units = "cm", width = 9, height = 6, scale = 2)

### DOWN
dpl_gsea_kegg_tat16h_dn <- gsea_kegg_tat16h %>% 
  filter(direction == "DOWN") %>% 
  ggplot(aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(color = p.adjust, size = Count)) +
  facet_wrap(~ direction, labeller = as_labeller(c("DOWN" = "Suppressed\nKEGG pathways"))) +
  labs(y = "") +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75), expand = expansion(add = 0.07)) +
  scale_size_continuous(breaks = size_breaks, limits = size_limits) +
  scale_color_gradient(
    low = "red", high = "blue",
    breaks = color_breaks, limits = color_limits) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(reverse = TRUE, order = 2)) +
  theme_dpl


ragg::agg_png(file.path(fig_dir, "DOWN_Tat.16h_vs_Tat.0h.png"), units = "cm", width = 8, height = 6, res = 600, scaling = 0.5)
dpl_gsea_kegg_tat16h_dn
dev.off()

ggsave(file.path(fig_dir, "DOWN_Tat.16h_vs_Tat.0h.pdf"), dpl_gsea_kegg_tat16h_dn, units = "cm", width = 8, height = 6, scale = 2)



### Tat.stable vs Tat.16h & Tat.16h vs Tat.0h ----------------------------------
gsea_kegg_tatstable_tat16h <- gsea_kegg_res %>% 
  filter(comparison %in% c("tatstable_tat16h", "tat16h_tat0h"))

### Transform data into wide format to plot lollipop graph 
gsea_kegg_tatstable_tat16h_wd <- gsea_kegg_tatstable_tat16h %>% 
  filter(p.adjust < 0.05) %>% # Filter KEGG pathways by p-value
  select(-c(core_enrichment, direction)) %>% 
  pivot_wider(names_from = comparison, values_from = c(NES, p.adjust, Count, GeneRatio))

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

