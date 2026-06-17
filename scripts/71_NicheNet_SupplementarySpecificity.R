#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

paper_font <- "Arial"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}
theme_set(theme_bw(base_size = 10.5, base_family = paper_font))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  project_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  project_root <- getwd()
}

target_dir <- file.path(project_root, "results", "human_nichenet_targetset_validation")
out_dir <- file.path(project_root, "results", "human_nichenet_targetset_validation", "supplementary_specificity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tests_path <- file.path(target_dir, "gse30122_nichenet_targetset_group_tests.csv")
if (!file.exists(tests_path)) {
  stop("Missing required file: ", tests_path)
}

tests <- readr::read_csv(tests_path, show_col_types = FALSE) %>%
  mutate(
    score_label = recode(
      score_type,
      mean_z_module_score = "Mean z-score",
      rank_ssgsea_like_score = "Rank-based ssGSEA-like score"
    ),
    analysis_label = recode(
      analysis,
      primary_tubules_only = "Primary tubules",
      all_tubular_including_control_reanalysis = "All tubular",
      glomerulus_specificity = "Glomerulus"
    ),
    analysis_label = factor(analysis_label, levels = c("Primary tubules", "All tubular", "Glomerulus")),
    score_label = factor(score_label, levels = c("Mean z-score", "Rank-based ssGSEA-like score")),
    label = sprintf("delta=%.3f\nP=%.3g\nFDR=%.3g", dkd_minus_control, wilcox_p, wilcox_fdr)
  )

readr::write_csv(tests, file.path(out_dir, "Supplementary_Table_GSE30122_targetset_specificity_stats.csv"))

p <- ggplot(tests, aes(x = analysis_label, y = score_label, fill = dkd_minus_control)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.0, lineheight = 0.95) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "GSE30122 compartment specificity of the NicheNet SPP1 target set",
    subtitle = "Positive delta indicates higher target-set score in DKD than control",
    x = NULL,
    y = NULL,
    fill = "DKD - Control"
  ) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(size = 13.5),
    plot.subtitle = element_text(size = 9.5)
  )

ggsave(file.path(out_dir, "Supplementary_Figure_S3_GSE30122_NicheNet_targetset_specificity.png"), p, width = 8.4, height = 4.2, dpi = 300)

summary_lines <- c(
  "GSE30122 compartment specificity summary for the NicheNet SPP1 target set",
  paste0("Date: ", Sys.Date()),
  "",
  "Primary tubules:",
  paste0("- Mean z-score delta=", sprintf("%.3f", tests$dkd_minus_control[tests$analysis == "primary_tubules_only" & tests$score_type == "mean_z_module_score"]),
         "; P=", signif(tests$wilcox_p[tests$analysis == "primary_tubules_only" & tests$score_type == "mean_z_module_score"], 3),
         "; FDR=", signif(tests$wilcox_fdr[tests$analysis == "primary_tubules_only" & tests$score_type == "mean_z_module_score"], 3)),
  paste0("- Rank-based ssGSEA-like delta=", sprintf("%.3f", tests$dkd_minus_control[tests$analysis == "primary_tubules_only" & tests$score_type == "rank_ssgsea_like_score"]),
         "; P=", signif(tests$wilcox_p[tests$analysis == "primary_tubules_only" & tests$score_type == "rank_ssgsea_like_score"], 3),
         "; FDR=", signif(tests$wilcox_fdr[tests$analysis == "primary_tubules_only" & tests$score_type == "rank_ssgsea_like_score"], 3)),
  "",
  "Glomerulus specificity check:",
  paste0("- Mean z-score delta=", sprintf("%.3f", tests$dkd_minus_control[tests$analysis == "glomerulus_specificity" & tests$score_type == "mean_z_module_score"]),
         "; P=", signif(tests$wilcox_p[tests$analysis == "glomerulus_specificity" & tests$score_type == "mean_z_module_score"], 3),
         "; FDR=", signif(tests$wilcox_fdr[tests$analysis == "glomerulus_specificity" & tests$score_type == "mean_z_module_score"], 3)),
  paste0("- Rank-based ssGSEA-like delta=", sprintf("%.3f", tests$dkd_minus_control[tests$analysis == "glomerulus_specificity" & tests$score_type == "rank_ssgsea_like_score"]),
         "; P=", signif(tests$wilcox_p[tests$analysis == "glomerulus_specificity" & tests$score_type == "rank_ssgsea_like_score"], 3),
         "; FDR=", signif(tests$wilcox_fdr[tests$analysis == "glomerulus_specificity" & tests$score_type == "rank_ssgsea_like_score"], 3)),
  "",
  "Interpretation:",
  "- The target-set projection was strongest in GSE30122 primary tubules.",
  "- Glomerular samples did not show a comparable significant target-set elevation.",
  "- This is a compartment-specificity check, not mechanistic proof of SPP1/CD44 signaling."
)
writeLines(enc2utf8(summary_lines), con = file.path(out_dir, "Supplementary_Figure_S3_specificity_summary.txt"), useBytes = TRUE)

message("Done. Outputs written to: ", out_dir)

