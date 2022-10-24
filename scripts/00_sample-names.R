## Sample names ----------------------------------------------------------------
sample_order <- c("Tat.stable", "Tat.16h", "Tat.0h", "control")

sample_names <- paste(rep(sample_order, each = 3), 1:3, sep = ".")

sample_names_format <- c(
  expression("RPMI"^{"Tat"}), 
  expression("RPMI"^{"Tat-ind"}~"Doxy+"),  
  expression("RPMI"^{"Tat-ind"}~"Doxy-"),
  expression("RPMI"))

sample_names_format_all <- c(
  expression("RPMI"^{"Tat"}~1), 
  expression("RPMI"^{"Tat"}~2), 
  expression("RPMI"^{"Tat"}~3),  
  expression("RPMI"^{"Tat-ind"}~"Doxy+"~1), 
  expression("RPMI"^{"Tat-ind"}~"Doxy+"~2), 
  expression("RPMI"^{"Tat-ind"}~"Doxy+"~3),
  expression("RPMI"^{"Tat-ind"}~"Doxy-"~1), 
  expression("RPMI"^{"Tat-ind"}~"Doxy-"~2), 
  expression("RPMI"^{"Tat-ind"}~"Doxy-"~3),
  expression("RPMI"~1), 
  expression("RPMI"~2), 
  expression("RPMI"~3))


## Comparison names ------------------------------------------------------------
comparison_order <- c(
  "Tat.stable_vs_control",
  "Tat.16h_vs_Tat.0h",
  "Tat.stable_vs_Tat.16h",
  "Tat.0h_vs_control")

comparison_names_format <- c(
  expression("RPMI"^{"Tat"}~"vs"~"RPMI"), 
  expression("RPMI"^{"Tat-ind"}~"Doxy+"~"vs"~"RPMI"^{"Tat-ind"}~"Doxy-"),  
  expression("RPMI"^{"Tat"}~"vs"~"RPMI"^{"Tat-ind"}~"Doxy+"),
  expression("RPMI"^{"Tat-ind"}~"Doxy-"~"vs"~"RPMI"))


## Color palette ---------------------------------------------------------------
cols4 <- c("#FE5E41", "#00A878", "#D8F1A0", "#F3C178")
cols4_sample <- c(
  "Tat.stable" = "#FE5E41", 
  "Tat.16h" = "#00A878", 
  "Tat.0h" = "#D8F1A0", 
  "control" = "#F3C178")

