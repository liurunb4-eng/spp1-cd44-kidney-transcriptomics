#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
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

out_dir <- file.path(project_root, "results", "human_nichenet_targetset_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nichenet_targets_path <- file.path(
  project_root,
  "results", "gse233078_nichenet_ligand_target",
  "nichenet_SPP1_predicted_MonoMacro_targets_top100.csv"
)
gse175759_dir <- file.path(project_root, "results", "human_bulk_gse175759_clinical_correlation")
gse30122_data_dir <- file.path(project_root, "data", "public_datasets", "DKD_bulk_candidate", "GSE30122")
gse30122_out_dir <- file.path(project_root, "results", "human_dkd_bulk_gse30122_signature_validation")

write_text <- function(path, lines) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

score_mean_z <- function(mat, genes) {
  genes <- unique(toupper(trimws(genes)))
  rn_upper <- toupper(rownames(mat))
  idx <- match(intersect(genes, rn_upper), rn_upper)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    return(list(score = rep(NA_real_, ncol(mat)), present = character(0)))
  }
  sub <- mat[idx, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  list(score = colMeans(z), present = rownames(mat)[idx])
}

score_rank_ssgsea_like <- function(mat, genes, alpha = 0.25) {
  genes <- unique(toupper(trimws(genes)))
  rn_upper <- toupper(rownames(mat))
  hit_index <- rn_upper %in% genes
  present <- rownames(mat)[hit_index]
  if (!any(hit_index)) {
    return(list(score = rep(NA_real_, ncol(mat)), present = character(0)))
  }
  scores <- apply(mat, 2, function(x) {
    ok <- is.finite(x)
    x <- x[ok]
    hits <- hit_index[ok]
    n <- length(x)
    nh <- sum(hits)
    nm <- n - nh
    if (nh == 0 || nm == 0) {
      return(NA_real_)
    }
    ranks <- rank(x, ties.method = "average")
    ord <- order(x, decreasing = TRUE)
    weights <- ranks[ord]^alpha
    hit_weights <- ifelse(hits[ord], weights, 0)
    phit <- cumsum(hit_weights) / sum(hit_weights)
    pmiss <- cumsum(!hits[ord]) / nm
    sum(phit - pmiss) / n
  })
  list(score = as.numeric(scores), present = present)
}

make_score_rows <- function(mat, genes, dataset) {
  mean_z <- score_mean_z(mat, genes)
  ss <- score_rank_ssgsea_like(mat, genes)
  sample_ids <- colnames(mat)
  if (is.null(sample_ids) || length(sample_ids) != ncol(mat)) {
    sample_ids <- paste0("sample_", seq_len(ncol(mat)))
  }
  scores <- bind_rows(
    tibble(dataset = dataset, gsm = sample_ids, score_type = "mean_z_module_score", score = mean_z$score),
    tibble(dataset = dataset, gsm = sample_ids, score_type = "rank_ssgsea_like_score", score = ss$score)
  )
  presence <- tibble(
    dataset = dataset,
    score_type = c("mean_z_module_score", "rank_ssgsea_like_score"),
    requested_n = length(unique(toupper(trimws(genes)))),
    present_n = c(length(mean_z$present), length(ss$present)),
    present_genes = c(paste(mean_z$present, collapse = ";"), paste(ss$present, collapse = ";")),
    missing_genes = c(
      paste(setdiff(unique(toupper(trimws(genes))), toupper(mean_z$present)), collapse = ";"),
      paste(setdiff(unique(toupper(trimws(genes))), toupper(ss$present)), collapse = ";")
    )
  )
  list(scores = scores, presence = presence)
}

targets <- readr::read_csv(nichenet_targets_path, show_col_types = FALSE) %>%
  mutate(gene_upper = toupper(gene_upper)) %>%
  filter(!is.na(gene_upper), gene_upper != "") %>%
  distinct(gene_upper, .keep_all = TRUE)

target_genes <- targets$gene_upper
readr::write_csv(
  tibble(signature = "SPP1_NicheNet_top100_targets", gene = target_genes),
  file.path(out_dir, "nichenet_spp1_top100_targetset_genes.csv")
)

message("Loading GSE175759 count-derived expression matrix...")
counts_path <- file.path(gse175759_dir, "gse175759_merged_counts_by_ensembl.csv.gz")
map_path <- file.path(gse175759_dir, "gse175759_ensembl_to_symbol_map.csv")
meta175_path <- file.path(gse175759_dir, "gse175759_metadata_clean.csv")
if (!file.exists(counts_path) || !file.exists(map_path) || !file.exists(meta175_path)) {
  stop("Missing GSE175759 processed files. Run 63_GSE175759_ClinicalBulkCorrelation.R first.")
}

counts <- readr::read_csv(counts_path, show_col_types = FALSE)
gene_map <- readr::read_csv(map_path, show_col_types = FALSE)
meta175 <- readr::read_csv(meta175_path, show_col_types = FALSE)
gse175_sample_cols <- setdiff(names(counts), "ensembl_version")
meta175 <- meta175 %>%
  mutate(usable = !is.na(egfr_ckd_epi) & technical_outlier == "No") %>%
  arrange(egfr_ckd_epi) %>%
  mutate(
    egfr_tertile = if_else(
      usable,
      c("Low eGFR", "Mid eGFR", "High eGFR")[dplyr::ntile(egfr_ckd_epi, 3)],
      NA_character_
    )
  ) %>%
  arrange(gsm)

expr175 <- counts %>%
  mutate(ensembl = sub("\\..*$", "", ensembl_version)) %>%
  left_join(gene_map, by = "ensembl") %>%
  filter(!is.na(symbol), symbol != "") %>%
  select(symbol, all_of(gse175_sample_cols)) %>%
  group_by(symbol) %>%
  summarise(across(everything(), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

mat175 <- as.matrix(expr175[, -1, drop = FALSE])
rownames(mat175) <- expr175$symbol
storage.mode(mat175) <- "numeric"
lib_size <- colSums(mat175)
cpm <- t(t(mat175) / lib_size * 1e6)
log_cpm175 <- log2(cpm + 1)

scores175 <- make_score_rows(log_cpm175, target_genes, "GSE175759")
gse175_scores <- scores175$scores %>%
  left_join(meta175, by = "gsm") %>%
  mutate(usable = !is.na(egfr_ckd_epi) & technical_outlier == "No")

cor_summary_one <- function(df) {
  df <- df %>% filter(is.finite(score), is.finite(egfr_ckd_epi))
  if (nrow(df) < 3) {
    return(tibble(
      n = nrow(df),
      spearman_rho = NA_real_,
      spearman_p = NA_real_,
      pearson_r = NA_real_,
      pearson_p = NA_real_
    ))
  }
  tibble(
    n = nrow(df),
    spearman_rho = suppressWarnings(cor(df$score, df$egfr_ckd_epi, method = "spearman")),
    spearman_p = suppressWarnings(cor.test(df$score, df$egfr_ckd_epi, method = "spearman", exact = FALSE)$p.value),
    pearson_r = suppressWarnings(cor(df$score, df$egfr_ckd_epi, method = "pearson")),
    pearson_p = suppressWarnings(cor.test(df$score, df$egfr_ckd_epi, method = "pearson")$p.value)
  )
}

gse175_cor <- gse175_scores %>%
  filter(usable) %>%
  group_by(score_type) %>%
  group_modify(~cor_summary_one(.x)) %>%
  ungroup() %>%
  mutate(
    spearman_fdr = p.adjust(spearman_p, method = "BH"),
    pearson_fdr = p.adjust(pearson_p, method = "BH")
  )

gse175_tertile_summary <- gse175_scores %>%
  filter(usable, !is.na(score), !is.na(egfr_tertile)) %>%
  group_by(score_type, egfr_tertile) %>%
  summarise(
    n = n(),
    mean_score = mean(score),
    median_score = median(score),
    sd_score = sd(score),
    mean_egfr = mean(egfr_ckd_epi),
    .groups = "drop"
  )

tertile_test_one <- function(df) {
  df <- df %>% filter(is.finite(score), egfr_tertile %in% c("Low eGFR", "High eGFR"))
  low <- df$score[df$egfr_tertile == "Low eGFR"]
  high <- df$score[df$egfr_tertile == "High eGFR"]
  p <- if (length(low) > 0 && length(high) > 0) {
    suppressWarnings(wilcox.test(low, high, exact = FALSE)$p.value)
  } else {
    NA_real_
  }
  tibble(
    low_n = length(low),
    high_n = length(high),
    low_mean = ifelse(length(low) > 0, mean(low), NA_real_),
    high_mean = ifelse(length(high) > 0, mean(high), NA_real_),
    low_minus_high = low_mean - high_mean,
    wilcox_p = p
  )
}

gse175_tertile_tests <- gse175_scores %>%
  filter(usable) %>%
  group_by(score_type) %>%
  group_modify(~tertile_test_one(.x)) %>%
  ungroup() %>%
  mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH"))

message("Loading GSE30122 expression matrix...")
series_path <- file.path(gse30122_data_dir, "GSE30122_series_matrix.txt.gz")
meta301_path <- file.path(gse30122_out_dir, "gse30122_metadata_parsed.csv")
if (!file.exists(series_path) || !file.exists(meta301_path)) {
  stop("Missing GSE30122 processed files. Run 68_GSE30122_HumanDKDBulkValidation.R first.")
}

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) {
    stop("Cannot locate series matrix table in ", path)
  }
  expr_text <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
  expr_df <- read.delim(text = expr_text, check.names = FALSE, stringsAsFactors = FALSE)
  names(expr_df)[1] <- "probe_id"
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(expr_mat) <- expr_df$probe_id
  storage.mode(expr_mat) <- "numeric"
  colnames(expr_mat) <- gsub('^"|"$', "", colnames(expr_mat))
  expr_mat
}

probe_mat <- parse_series_matrix(series_path)
meta301 <- readr::read_csv(meta301_path, show_col_types = FALSE)

message("Loading GPL571 annotation...")
gpl <- GEOquery::getGEO("GPL571", destdir = gse30122_data_dir)
gpl_table <- GEOquery::Table(gpl) %>%
  as_tibble() %>%
  transmute(
    probe_id = ID,
    gene_symbol_raw = `Gene Symbol`
  ) %>%
  filter(!is.na(gene_symbol_raw), gene_symbol_raw != "")

probe_map <- gpl_table %>%
  mutate(gene_symbol_raw = gsub(" /// ", ";", gene_symbol_raw, fixed = TRUE)) %>%
  separate_rows(gene_symbol_raw, sep = ";") %>%
  mutate(symbol = toupper(trimws(gene_symbol_raw))) %>%
  filter(symbol != "", !grepl("///", symbol, fixed = TRUE)) %>%
  distinct(probe_id, symbol)

gene_expr301 <- as.data.frame(probe_mat) %>%
  rownames_to_column("probe_id") %>%
  inner_join(probe_map, by = "probe_id") %>%
  select(symbol, all_of(colnames(probe_mat))) %>%
  group_by(symbol) %>%
  summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

mat301 <- as.matrix(gene_expr301[, -1, drop = FALSE])
rownames(mat301) <- gene_expr301$symbol
storage.mode(mat301) <- "numeric"

scores301 <- make_score_rows(mat301, target_genes, "GSE30122")
gse301_scores <- scores301$scores %>%
  left_join(meta301, by = "gsm") %>%
  mutate(
    disease_binary = case_when(
      disease_state == "diabetic kidney disease (DKD)" ~ "DKD",
      disease_state == "control" ~ "Control",
      TRUE ~ NA_character_
    ),
    is_primary_tubules = tissue_subregion == "tubules",
    is_tubular_any = tissue_subregion %in% c("tubules", "tubulus")
  )

run_group_tests <- function(df, analysis_name) {
  group_test_one <- function(x) {
    x <- x %>% filter(disease_binary %in% c("Control", "DKD"), is.finite(score))
    dkd <- x$score[x$disease_binary == "DKD"]
    ctrl <- x$score[x$disease_binary == "Control"]
    pooled <- if (length(dkd) > 1 && length(ctrl) > 1) {
      sqrt(((length(dkd) - 1) * var(dkd) + (length(ctrl) - 1) * var(ctrl)) / (length(dkd) + length(ctrl) - 2))
    } else {
      NA_real_
    }
    tibble(
      analysis = analysis_name,
      control_n = length(ctrl),
      dkd_n = length(dkd),
      control_mean = ifelse(length(ctrl) > 0, mean(ctrl), NA_real_),
      dkd_mean = ifelse(length(dkd) > 0, mean(dkd), NA_real_),
      dkd_minus_control = dkd_mean - control_mean,
      wilcox_p = ifelse(length(ctrl) > 0 && length(dkd) > 0, suppressWarnings(wilcox.test(dkd, ctrl, exact = FALSE)$p.value), NA_real_),
      cohen_d = ifelse(is.finite(pooled) && pooled > 0, (mean(dkd) - mean(ctrl)) / pooled, NA_real_)
    )
  }
  df %>%
    group_by(score_type) %>%
    group_modify(~group_test_one(.x)) %>%
    ungroup() %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH"))
}

gse301_tests <- bind_rows(
  run_group_tests(gse301_scores %>% filter(is_primary_tubules), "primary_tubules_only"),
  run_group_tests(gse301_scores %>% filter(is_tubular_any), "all_tubular_including_control_reanalysis"),
  run_group_tests(gse301_scores %>% filter(tissue_subregion == "glomerulus"), "glomerulus_specificity")
)

presence <- bind_rows(scores175$presence, scores301$presence)
readr::write_csv(presence, file.path(out_dir, "nichenet_spp1_targetset_gene_presence.csv"))
readr::write_csv(gse175_scores, file.path(out_dir, "gse175759_nichenet_targetset_scores.csv"))
readr::write_csv(gse175_cor, file.path(out_dir, "gse175759_nichenet_targetset_egfr_correlations.csv"))
readr::write_csv(gse175_tertile_summary, file.path(out_dir, "gse175759_nichenet_targetset_egfr_tertile_summary.csv"))
readr::write_csv(gse175_tertile_tests, file.path(out_dir, "gse175759_nichenet_targetset_egfr_tertile_tests.csv"))
readr::write_csv(gse301_scores, file.path(out_dir, "gse30122_nichenet_targetset_scores.csv"))
readr::write_csv(gse301_tests, file.path(out_dir, "gse30122_nichenet_targetset_group_tests.csv"))

score_labels <- c(
  mean_z_module_score = "Mean z-score",
  rank_ssgsea_like_score = "Rank-based ssGSEA-like score"
)

p1 <- gse175_scores %>%
  filter(usable, !is.na(score)) %>%
  mutate(score_type = factor(score_type, levels = names(score_labels), labels = score_labels)) %>%
  ggplot(aes(x = score, y = egfr_ckd_epi, color = diagnosis)) +
  geom_point(size = 2.1, alpha = 0.82) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
  facet_wrap(~score_type, scales = "free_x", ncol = 2) +
  labs(
    title = "GSE175759: NicheNet SPP1 target-set score versus eGFR",
    subtitle = "Top 100 SPP1-prioritized Mono/Macro receiver targets projected into human CKD tubulointerstitium",
    x = "Target-set score",
    y = "eGFR (CKD-EPI)",
    color = "Diagnosis"
  )
ggsave(file.path(out_dir, "Figure_NicheNet_targetset_GSE175759_eGFR_scatter.png"), p1, width = 9.0, height = 4.8, dpi = 300)

p2 <- gse175_scores %>%
  filter(usable, !is.na(score), !is.na(egfr_tertile)) %>%
  mutate(
    egfr_tertile = factor(egfr_tertile, levels = c("Low eGFR", "Mid eGFR", "High eGFR")),
    score_type = factor(score_type, levels = names(score_labels), labels = score_labels)
  ) %>%
  ggplot(aes(x = egfr_tertile, y = score, fill = egfr_tertile)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.8, alpha = 0.82) +
  facet_wrap(~score_type, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Low eGFR" = "#B2182B", "Mid eGFR" = "#777777", "High eGFR" = "#2166AC")) +
  labs(
    title = "GSE175759: NicheNet SPP1 target-set score across eGFR tertiles",
    x = NULL,
    y = "Target-set score",
    fill = NULL
  ) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "Figure_NicheNet_targetset_GSE175759_eGFR_tertiles.png"), p2, width = 8.2, height = 4.8, dpi = 300)

p3 <- gse301_scores %>%
  filter(is_primary_tubules, disease_binary %in% c("Control", "DKD"), !is.na(score)) %>%
  mutate(
    disease_binary = factor(disease_binary, levels = c("Control", "DKD")),
    score_type = factor(score_type, levels = names(score_labels), labels = score_labels)
  ) %>%
  ggplot(aes(x = disease_binary, y = score, fill = disease_binary)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.0, alpha = 0.86) +
  facet_wrap(~score_type, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Control" = "#4D4D4D", "DKD" = "#B2182B")) +
  labs(
    title = "GSE30122: NicheNet SPP1 target-set score in human DKD tubules",
    subtitle = "Primary tubules comparison: DKD n=10 vs control n=12",
    x = NULL,
    y = "Target-set score",
    fill = NULL
  ) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "Figure_NicheNet_targetset_GSE30122_DKD_tubules.png"), p3, width = 8.2, height = 4.8, dpi = 300)

heat_df <- bind_rows(
  gse175_cor %>%
    transmute(dataset_analysis = "GSE175759 eGFR", score_type, effect = spearman_rho, p = spearman_p, label = sprintf("rho=%.2f\np=%.3g", effect, p)),
  gse301_tests %>%
    filter(analysis == "primary_tubules_only") %>%
    transmute(dataset_analysis = "GSE30122 DKD tubules", score_type, effect = dkd_minus_control, p = wilcox_p, label = sprintf("delta=%.2f\np=%.3g", effect, p))
) %>%
  mutate(
    score_type = factor(score_type, levels = names(score_labels), labels = score_labels),
    dataset_analysis = factor(dataset_analysis, levels = c("GSE175759 eGFR", "GSE30122 DKD tubules"))
  )

p4 <- ggplot(heat_df, aes(x = dataset_analysis, y = score_type, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Human projection of the NicheNet SPP1 target set",
    subtitle = "Negative rho: higher score with lower eGFR; positive delta: higher score in DKD",
    x = NULL,
    y = NULL,
    fill = "Effect"
  ) +
  theme(
    axis.text.x = element_text(angle = 12, hjust = 1),
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(size = 9.5)
  )
ggsave(file.path(out_dir, "Figure_NicheNet_targetset_human_validation_summary.png"), p4, width = 8.8, height = 3.9, dpi = 300)

top_cor <- gse175_cor %>% arrange(spearman_p)
top_dkd <- gse301_tests %>% filter(analysis == "primary_tubules_only") %>% arrange(wilcox_p)
summary_lines <- c(
  "Human validation of the NicheNet SPP1-prioritized Mono/Macro target set",
  paste0("Date: ", Sys.Date()),
  "",
  "Input target set:",
  paste0("- NicheNet SPP1-prioritized top targets requested: ", length(target_genes)),
  paste0("- GSE175759 present genes: ", scores175$presence$present_n[1], "/", scores175$presence$requested_n[1]),
  paste0("- GSE30122 present genes: ", scores301$presence$present_n[1], "/", scores301$presence$requested_n[1]),
  "",
  "GSE175759 human CKD tubulointerstitium/eGFR association:",
  paste0(
    "- Mean z-score: n=", top_cor$n[top_cor$score_type == "mean_z_module_score"],
    "; Spearman rho=", sprintf("%.3f", top_cor$spearman_rho[top_cor$score_type == "mean_z_module_score"]),
    "; P=", signif(top_cor$spearman_p[top_cor$score_type == "mean_z_module_score"], 3),
    "; FDR=", signif(top_cor$spearman_fdr[top_cor$score_type == "mean_z_module_score"], 3)
  ),
  paste0(
    "- Rank-based ssGSEA-like score: n=", top_cor$n[top_cor$score_type == "rank_ssgsea_like_score"],
    "; Spearman rho=", sprintf("%.3f", top_cor$spearman_rho[top_cor$score_type == "rank_ssgsea_like_score"]),
    "; P=", signif(top_cor$spearman_p[top_cor$score_type == "rank_ssgsea_like_score"], 3),
    "; FDR=", signif(top_cor$spearman_fdr[top_cor$score_type == "rank_ssgsea_like_score"], 3)
  ),
  "",
  "GSE30122 human DKD primary tubules:",
  paste0(
    "- Mean z-score: DKD-Control delta=", sprintf("%.3f", top_dkd$dkd_minus_control[top_dkd$score_type == "mean_z_module_score"]),
    "; Wilcoxon P=", signif(top_dkd$wilcox_p[top_dkd$score_type == "mean_z_module_score"], 3),
    "; FDR=", signif(top_dkd$wilcox_fdr[top_dkd$score_type == "mean_z_module_score"], 3)
  ),
  paste0(
    "- Rank-based ssGSEA-like score: DKD-Control delta=", sprintf("%.3f", top_dkd$dkd_minus_control[top_dkd$score_type == "rank_ssgsea_like_score"]),
    "; Wilcoxon P=", signif(top_dkd$wilcox_p[top_dkd$score_type == "rank_ssgsea_like_score"], 3),
    "; FDR=", signif(top_dkd$wilcox_fdr[top_dkd$score_type == "rank_ssgsea_like_score"], 3)
  ),
  "",
  "Interpretation guardrail:",
  "- This analysis projects the NicheNet-prioritized SPP1 receiver-target set into human bulk cohorts.",
  "- It should be interpreted as cross-dataset target-program contextualization, not as validation of causal SPP1/CD44 signaling.",
  "- Rank-based scores are reported as ssGSEA-like enrichment because the GSVA package was not available in the local R environment."
)
write_text(file.path(out_dir, "human_nichenet_targetset_validation_summary.txt"), summary_lines)

write_text(
  file.path(out_dir, "status.txt"),
  c("human_nichenet_targetset_validation completed", paste0("timestamp: ", Sys.time()))
)

message("Done. Outputs written to: ", out_dir)

