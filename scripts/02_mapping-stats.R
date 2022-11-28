## Load libraries --------------------------------------------------------------
library(tidyverse)
### Import sample names notations and color palette
source("scripts/00_sample-names.R")

input_dir <- "data/summaries"
fig_dir <- "results/figures"


## Load data -------------------------------------------------------------------
summary_df <- read_tsv(file.path(input_dir, "mapping_stats.txt"))

### Transform dataframe to long format for plotting
summary_lng <- summary_df %>% 
  pivot_longer(cols = -sample) %>% 
  mutate(
    sample = factor(sample, levels = sample_names),
    name = factor(name, levels = c("raw", "trimmed", "unique_aligned", "mapped_gencode")))


## Plot ------------------------------------------------------------------------
summary_pl <- ggplot(summary_lng, aes(x = value, y = fct_rev(sample), fill = fct_rev(name))) +
  geom_col(width = 0.7, alpha = 0.7, position = position_dodge()) +
  scale_x_continuous(labels = scales::format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)) +
  scale_y_discrete(labels = rev(sample_names_format_all)) +
  scale_fill_manual(
    values = c(
      "raw" = cols4[1], 
      "trimmed" = cols4[2], 
      "unique_aligned" = cols4[3], 
      "mapped_gencode" = cols4[4]),
    labels = c(
      "Total reads", 
      "After trimming", 
      "Uniquely aligned to genome",
      "Uniquely mapped to GENCODE annotation")) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "# reads", y = "", fill = "") +
  ## Edit theme
  theme(text = element_text(size = 12, color = "black"), 
        axis.text = element_text(color = "black"),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top",
        aspect.ratio = 2/3)


ggsave(file.path(fig_dir, "mapping_stats.png"), summary_pl, units = "cm", width = 15)
ggsave(file.path(fig_dir, "mapping_stats.pdf"), summary_pl, units = "cm", width = 15)


