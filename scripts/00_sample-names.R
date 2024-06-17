## Sample names ----------------------------------------------------------------
sample_order <- c("control", "Tat.stable", "Tat.0h", "Tat.16h")

sample_names <- paste(rep(sample_order, each = 3), 1:3, sep = ".")

sample_names_format <- c(
  expression("RPMI 8866"),
  expression("RPMI"^{"Tat"}), 
  expression("RPMI"^{"iTat dox-"}),
  expression("RPMI"^{"iTat dox+"}))

sample_names_format_all <- c(
  expression("RPMI 8866"~1), 
  expression("RPMI 8866"~2), 
  expression("RPMI 8866"~3),
  expression("RPMI"^{"Tat"}~1), 
  expression("RPMI"^{"Tat"}~2), 
  expression("RPMI"^{"Tat"}~3),  
  expression("RPMI"^{"iTat dox-"}~1), 
  expression("RPMI"^{"iTat dox-"}~2), 
  expression("RPMI"^{"iTat dox-"}~3),
  expression("RPMI"^{"iTat dox+"}~1), 
  expression("RPMI"^{"iTat dox+"}~2), 
  expression("RPMI"^{"iTat dox+"}~3))


## Comparison names ------------------------------------------------------------
comparison_order <- c(
  "Tat.0h_vs_control",
  "Tat.stable_vs_control",
  "Tat.16h_vs_Tat.0h",
  "Tat.stable_vs_Tat.16h")

comparison_names_format <- c(
  expression("RPMI"^{"iTat dox-"}~"vs"~"RPMI 8866"),
  expression("RPMI"^{"Tat"}~"vs"~"RPMI 8866"), 
  expression("RPMI"^{"iTat dox+"}~"vs"~"RPMI"^{"iTat dox-"}),  
  expression("RPMI"^{"Tat"}~"vs"~"RPMI"^{"iTat dox+"}))


## Color palette ---------------------------------------------------------------
cols4 <- c("#FE5E41", "#00A878", "#D8F1A0", "#F3C178")
cols4_sample <- c(
  "Tat.stable" = "#FE5E41", 
  "Tat.16h" = "#00A878", 
  "Tat.0h" = "#D8F1A0", 
  "control" = "#F3C178")

