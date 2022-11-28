## Load libraries --------------------------------------------------------------
library(tidyverse)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load DESeq object -----------------------------------------------------------
dds <- readRDS(file.path(input_dir, "dds.rds"))


## Visualization of normalization effect on count distribution -----------------
### Import raw and normalized counts
counts_flt <- as.data.frame(counts(dds, normalized = FALSE))
counts_norm <- as.data.frame(counts(dds, normalized = TRUE))

counts_flt_long <- counts_flt %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "count") %>% 
  mutate(
    count = log2(count + 1),
    state = "Before normalization")

counts_norm_long <- counts_norm %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "count") %>% 
  mutate(
    count = log2(count + 1),
    state = "After normalization")

counts_long <- rbind(counts_flt_long, counts_norm_long) %>% 
  mutate(
    sample = factor(sample, levels = sample_names),
    state = factor(state, levels = c("Before normalization", "After normalization")))


norm_pl <- ggplot(counts_long, aes(x = count, y = fct_rev(sample))) + 
  geom_boxplot() + 
  facet_wrap(~state) +
  labs(x = expression("log"[10]*"(count+1)"), y = "") +
  scale_y_discrete(labels = rev(sample_names_format_all)) +
  ## Edit theme
  theme(text = element_text(size = 12, color = "black"), 
        axis.text = element_text(color = "black"),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white", colour = "black", size = 0.5, linetype = "solid"),
        strip.text = element_text(face="bold", size = 10), 
        legend.position = "top",
        aspect.ratio = 1)

ggsave(file.path(fig_dir, "DESeq_normalization.png"), norm_pl, units = "cm", width = 10, scale = 3/2)
ggsave(file.path(fig_dir, "DESeq_normalization.pdf"), norm_pl, units = "cm", width = 10, scale = 3/2)

