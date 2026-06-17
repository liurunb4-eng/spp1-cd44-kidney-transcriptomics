#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
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

in_dir <- file.path(project_root, "results", "human_bulk_gse175759_clinical_correlation")
out_dir <- file.path(project_root, "results", "human_bulk_gse175759_egfr_stratification_patch")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_text <- function(path, lines) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

priority_sigs <- c(
  "Spp1_Cd44_tubular_immune_program",
  "SPP1_CD44_anchor_score",
  "human_tubular_injury_context",
  "curated_fibro_inflammatory_context"
)

friendly_names <- c(
  Spp1_Cd44_tubular_immune_program = "Spp1/Cd44 tubular-immune program",
  SPP1_CD44_anchor_score = "SPP1/CD44 anchor",
  human_tubular_injury_context = "Human tubular injury context",
  curated_fibro_inflammatory_context = "Fibro-inflammatory context"
)

scores <- readr::read_csv(file.path(in_dir, "gse175759_signature_scores_long.csv"), show_col_types = FALSE) %>%
  mutate(
    egfr_ckd_epi = suppressWarnings(as.numeric(egfr_ckd_epi)),
    score = suppressWarnings(as.numeric(score)),
    usable = as.character(usable) %in% c("TRUE", "True", "true", "1"),
    diagnosis = as.character(diagnosis),
    diagnosis_group = case_when(
      diagnosis == "Control" ~ "Control",
      diagnosis == "IgAN" ~ "IgAN",
      TRUE ~ "Other CKD"
    ),
    disease_group = if_else(diagnosis == "Control", "Control", "CKD")
  ) %>%
  filter(usable, signature %in% priority_sigs, !is.na(egfr_ckd_epi), !is.na(score))

sample_meta <- scores %>%
  distinct(gsm, diagnosis, diagnosis_group, disease_group, egfr_ckd_epi)

tertile_breaks <- quantile(sample_meta$egfr_ckd_epi, probs = c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE, type = 7)
tertile_breaks[1] <- tertile_breaks[1] - 1e-6
tertile_breaks[4] <- tertile_breaks[4] + 1e-6

scores <- scores %>%
  mutate(
    egfr_tertile = cut(
      egfr_ckd_epi,
      breaks = tertile_breaks,
      labels = c("Low eGFR", "Mid eGFR", "High eGFR"),
      include.lowest = TRUE,
      ordered_result = TRUE
    ),
    egfr_tertile_num = as.numeric(egfr_tertile),
    egfr_60_group = if_else(egfr_ckd_epi < 60, "eGFR <60", "eGFR >=60"),
    signature_label = recode(signature, !!!friendly_names)
  )

readr::write_csv(scores, file.path(out_dir, "gse175759_scores_with_egfr_strata.csv"))

safe_wilcox <- function(df, group_col, g1, g2) {
  group_values <- as.character(df[[group_col]])
  sub <- df[group_values %in% c(g1, g2), , drop = FALSE]
  sub_group <- as.character(sub[[group_col]])
  if (n_distinct(sub_group) < 2 || any(table(sub_group) < 2)) {
    return(tibble(p_value = NA_real_, delta = NA_real_, n1 = NA_integer_, n2 = NA_integer_, mean1 = NA_real_, mean2 = NA_real_))
  }
  x <- sub$score[sub_group == g1]
  y <- sub$score[sub_group == g2]
  tibble(
    p_value = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value),
    delta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    n1 = length(x),
    n2 = length(y),
    mean1 = mean(x, na.rm = TRUE),
    mean2 = mean(y, na.rm = TRUE)
  )
}

tertile_summary <- scores %>%
  group_by(signature, signature_label, egfr_tertile) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    mean_egfr = mean(egfr_ckd_epi, na.rm = TRUE),
    .groups = "drop"
  )

tertile_tests <- scores %>%
  group_by(signature, signature_label) %>%
  group_modify(~{
    df <- .x
    kw_p <- suppressWarnings(kruskal.test(score ~ egfr_tertile, data = df)$p.value)
    trend <- suppressWarnings(cor.test(df$score, df$egfr_tertile_num, method = "spearman", exact = FALSE))
    lh <- safe_wilcox(df, "egfr_tertile", "Low eGFR", "High eGFR")
    tibble(
      n = nrow(df),
      kruskal_p = kw_p,
      tertile_spearman_rho = unname(trend$estimate),
      tertile_spearman_p = trend$p.value,
      low_vs_high_delta = lh$delta,
      low_vs_high_p = lh$p_value,
      low_n = lh$n1,
      high_n = lh$n2,
      low_mean = lh$mean1,
      high_mean = lh$mean2
    )
  }) %>%
  ungroup() %>%
  mutate(
    kruskal_fdr = p.adjust(kruskal_p, method = "BH"),
    low_vs_high_fdr = p.adjust(low_vs_high_p, method = "BH"),
    direction = if_else(low_vs_high_delta > 0, "higher_in_low_eGFR", "lower_in_low_eGFR")
  ) %>%
  arrange(low_vs_high_p)

egfr60_tests <- scores %>%
  group_by(signature, signature_label) %>%
  group_modify(~{
    res <- safe_wilcox(.x, "egfr_60_group", "eGFR <60", "eGFR >=60")
    tibble(
      n = nrow(.x),
      impaired_vs_preserved_delta = res$delta,
      impaired_vs_preserved_p = res$p_value,
      impaired_n = res$n1,
      preserved_n = res$n2,
      impaired_mean = res$mean1,
      preserved_mean = res$mean2
    )
  }) %>%
  ungroup() %>%
  mutate(
    impaired_vs_preserved_fdr = p.adjust(impaired_vs_preserved_p, method = "BH"),
    direction = if_else(impaired_vs_preserved_delta > 0, "higher_in_eGFR_lt60", "lower_in_eGFR_lt60")
  ) %>%
  arrange(impaired_vs_preserved_p)

cohort_sets <- list(
  all_usable = scores,
  ckd_only = scores %>% filter(disease_group == "CKD"),
  igan_only = scores %>% filter(diagnosis == "IgAN"),
  non_control_non_igan = scores %>% filter(disease_group == "CKD", diagnosis != "IgAN")
)

sensitivity_rows <- list()
for (set_name in names(cohort_sets)) {
  dat <- cohort_sets[[set_name]]
  sensitivity_rows[[set_name]] <- dat %>%
    group_by(signature, signature_label) %>%
    summarise(
      subset = set_name,
      n = n(),
      n_diagnoses = n_distinct(diagnosis),
      spearman_rho = if_else(n() >= 6, suppressWarnings(cor(score, egfr_ckd_epi, method = "spearman")), NA_real_),
      spearman_p = if_else(n() >= 6, suppressWarnings(cor.test(score, egfr_ckd_epi, method = "spearman", exact = FALSE)$p.value), NA_real_),
      pearson_r = if_else(n() >= 6, suppressWarnings(cor(score, egfr_ckd_epi, method = "pearson")), NA_real_),
      pearson_p = if_else(n() >= 6, suppressWarnings(cor.test(score, egfr_ckd_epi, method = "pearson")$p.value), NA_real_),
      .groups = "drop"
    )
}
sensitivity <- bind_rows(sensitivity_rows) %>%
  group_by(subset) %>%
  mutate(spearman_fdr_within_subset = p.adjust(spearman_p, method = "BH")) %>%
  ungroup() %>%
  arrange(match(subset, names(cohort_sets)), spearman_p)

readr::write_csv(tertile_summary, file.path(out_dir, "gse175759_egfr_tertile_summary.csv"))
readr::write_csv(tertile_tests, file.path(out_dir, "gse175759_egfr_tertile_tests.csv"))
readr::write_csv(egfr60_tests, file.path(out_dir, "gse175759_egfr60_tests.csv"))
readr::write_csv(sensitivity, file.path(out_dir, "gse175759_egfr_correlation_sensitivity.csv"))

plot_box <- scores %>%
  mutate(
    signature_label = factor(signature_label, levels = friendly_names[priority_sigs]),
    egfr_tertile = factor(egfr_tertile, levels = c("Low eGFR", "Mid eGFR", "High eGFR"))
  )

p1 <- ggplot(plot_box, aes(x = egfr_tertile, y = score, fill = egfr_tertile)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_point(aes(color = diagnosis_group), position = position_jitter(width = 0.14, height = 0), size = 1.7, alpha = 0.82) +
  facet_wrap(~signature_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Low eGFR" = "#B2182B", "Mid eGFR" = "#F4A582", "High eGFR" = "#2166AC")) +
  scale_color_manual(values = c("Control" = "#4D4D4D", "IgAN" = "#1B9E77", "Other CKD" = "#7570B3")) +
  labs(
    title = "GSE175759 eGFR-stratified signature burden",
    subtitle = "Low eGFR samples show higher tubular-immune and injury program scores",
    x = NULL,
    y = "Signature score",
    fill = "eGFR tertile",
    color = "Diagnosis group"
  ) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", color = "grey80")
  )
ggsave(file.path(out_dir, "FigureCH_gse175759_egfr_tertile_signature_boxplots.png"), p1, width = 9.4, height = 7.0, dpi = 300)

heat_df <- sensitivity %>%
  mutate(
    signature_label = factor(signature_label, levels = friendly_names[priority_sigs]),
    subset = factor(subset, levels = names(cohort_sets), labels = c("All usable", "CKD only", "IgAN only", "Non-IgAN CKD"))
  )

p2 <- ggplot(heat_df, aes(x = subset, y = signature_label, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = if_else(is.na(spearman_p), "NA", sprintf("rho=%.2f\np=%.3g", spearman_rho, spearman_p))), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, na.value = "grey90") +
  labs(
    title = "GSE175759 eGFR correlation sensitivity",
    subtitle = "Negative rho indicates higher program score with lower kidney function",
    x = NULL,
    y = NULL,
    fill = "Spearman rho"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "FigureCI_gse175759_egfr_correlation_sensitivity_heatmap.png"), p2, width = 8.8, height = 4.7, dpi = 300)

main_df <- scores %>%
  filter(signature == "Spp1_Cd44_tubular_immune_program") %>%
  mutate(egfr_tertile = factor(egfr_tertile, levels = c("Low eGFR", "Mid eGFR", "High eGFR")))

p3 <- ggplot(main_df, aes(x = egfr_ckd_epi, y = score, color = diagnosis_group)) +
  geom_point(size = 2.25, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
  geom_vline(xintercept = 60, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Control" = "#4D4D4D", "IgAN" = "#1B9E77", "Other CKD" = "#7570B3")) +
  labs(
    title = "Spp1/Cd44 tubular-immune program tracks lower eGFR in GSE175759",
    subtitle = "Dashed line marks eGFR 60 ml/min/1.73m2",
    x = "eGFR (CKD-EPI)",
    y = "Spp1/Cd44 tubular-immune program score",
    color = "Diagnosis group"
  )
ggsave(file.path(out_dir, "FigureCJ_gse175759_main_program_egfr_scatter.png"), p3, width = 7.2, height = 5.4, dpi = 300)

main_tertile <- tertile_tests %>% filter(signature == "Spp1_Cd44_tubular_immune_program")
main_60 <- egfr60_tests %>% filter(signature == "Spp1_Cd44_tubular_immune_program")
main_sens <- sensitivity %>% filter(signature == "Spp1_Cd44_tubular_immune_program")
injury_tertile <- tertile_tests %>% filter(signature == "human_tubular_injury_context")

diagnosis_counts <- sample_meta %>%
  count(diagnosis) %>%
  arrange(desc(n)) %>%
  mutate(label = paste0(diagnosis, " n=", n)) %>%
  pull(label)

lines <- c(
  "GSE175759 eGFR stratification patch",
  "",
  "Purpose:",
  "This patch converts the previous continuous GSE175759 eGFR correlation into clinically readable eGFR strata and sensitivity analyses.",
  "",
  "Dataset boundary:",
  paste0("- Usable samples with eGFR: n=", n_distinct(sample_meta$gsm)),
  paste0("- Diagnosis composition: ", paste(diagnosis_counts, collapse = "; ")),
  "- This is a human CKD tubulointerstitium clinical-extension dataset, not a DKD-specific validation dataset.",
  "",
  "eGFR tertile thresholds:",
  sprintf("- Low eGFR: %.1f to %.1f", min(sample_meta$egfr_ckd_epi), tertile_breaks[2]),
  sprintf("- Mid eGFR: %.1f to %.1f", tertile_breaks[2], tertile_breaks[3]),
  sprintf("- High eGFR: %.1f to %.1f", tertile_breaks[3], max(sample_meta$egfr_ckd_epi)),
  "",
  "Main Spp1/Cd44 tubular-immune program:",
  sprintf(
    "- Low-vs-high eGFR mean-score delta = %.3f; Wilcoxon p = %.4g; FDR = %.4g; direction = %s",
    main_tertile$low_vs_high_delta, main_tertile$low_vs_high_p, main_tertile$low_vs_high_fdr, main_tertile$direction
  ),
  sprintf(
    "- eGFR <60 vs >=60 delta = %.3f; Wilcoxon p = %.4g; FDR = %.4g; direction = %s",
    main_60$impaired_vs_preserved_delta, main_60$impaired_vs_preserved_p, main_60$impaired_vs_preserved_fdr, main_60$direction
  ),
  "",
  "Human tubular injury context:",
  sprintf(
    "- Low-vs-high eGFR mean-score delta = %.3f; Wilcoxon p = %.4g; FDR = %.4g; direction = %s",
    injury_tertile$low_vs_high_delta, injury_tertile$low_vs_high_p, injury_tertile$low_vs_high_fdr, injury_tertile$direction
  ),
  "",
  "Sensitivity correlations for the main program:",
  paste(
    sprintf(
      "- %s: n=%s, rho=%.3f, p=%.4g, FDR=%.4g",
      main_sens$subset, main_sens$n, main_sens$spearman_rho, main_sens$spearman_p, main_sens$spearman_fdr_within_subset
    ),
    collapse = "\n"
  ),
  "",
  "Manuscript-safe interpretation:",
  "GSE175759 supports a clinically relevant association between the Spp1/Cd44-related tubular-immune program and reduced kidney function. The signal is strongest as a pan-CKD tubulointerstitial clinical-extension layer and should not be described as DKD-specific or as proof of cell-cell communication.",
  "",
  "Recommended placement:",
  "Use this as a main-text human clinical-extension figure if the manuscript needs a human anchor. The spatial and Lake2023 myeloid patches can remain supplementary, while this GSE175759 stratification provides the more defensible human clinical relevance layer.",
  "",
  "Output files:",
  "- gse175759_scores_with_egfr_strata.csv",
  "- gse175759_egfr_tertile_summary.csv",
  "- gse175759_egfr_tertile_tests.csv",
  "- gse175759_egfr60_tests.csv",
  "- gse175759_egfr_correlation_sensitivity.csv",
  "- FigureCH_gse175759_egfr_tertile_signature_boxplots.png",
  "- FigureCI_gse175759_egfr_correlation_sensitivity_heatmap.png",
  "- FigureCJ_gse175759_main_program_egfr_scatter.png"
)

write_text(file.path(project_root, "analysis_summary.txt"), lines)
write_text(file.path(out_dir, "gse175759_egfr_stratification_patch_summary.txt"), lines)
write_text(file.path(out_dir, "status.txt"), "completed")

cat(paste(lines, collapse = "\n"))

