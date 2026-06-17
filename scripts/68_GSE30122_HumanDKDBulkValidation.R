#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
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

data_dir <- file.path(project_root, "data", "public_datasets", "DKD_bulk_candidate", "GSE30122")
out_dir <- file.path(project_root, "results", "human_dkd_bulk_gse30122_signature_validation")
sig_path <- file.path(project_root, "results", "bulk_signature_projection_spp1_cd44_program", "signature_definitions.csv")
series_path <- file.path(data_dir, "GSE30122_series_matrix.txt.gz")
gpl_dir <- data_dir
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_text <- function(path, lines) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(invisible(TRUE))
  }
  cmd <- sprintf('curl.exe -L "%s" -o "%s"', url, dest)
  status <- system(cmd)
  if (!file.exists(dest) || file.info(dest)$size == 0 || status != 0) {
    stop("Failed to download: ", url)
  }
  invisible(TRUE)
}

if (!file.exists(series_path) || file.info(series_path)$size == 0) {
  download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30122/matrix/GSE30122_series_matrix.txt.gz",
    series_path
  )
}

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) {
    stop("Cannot locate series matrix table in ", path)
  }

  sample_rows <- lines[grep("^!Sample_", lines)]
  row_names <- sub("\t.*$", "", sample_rows)
  row_names <- sub("^!Sample_", "", row_names)
  row_counts <- ave(seq_along(row_names), row_names, FUN = seq_along)
  unique_names <- ifelse(duplicated(row_names) | duplicated(row_names, fromLast = TRUE), paste0(row_names, "_", row_counts), row_names)

  sample_meta <- NULL
  for (i in seq_along(sample_rows)) {
    parts <- strsplit(sample_rows[[i]], "\t", fixed = TRUE)[[1]]
    values <- gsub('^"|"$', "", parts[-1])
    if (is.null(sample_meta)) {
      sample_meta <- tibble(sample_index = seq_along(values))
    }
    sample_meta[[unique_names[[i]]]] <- values
  }
  sample_meta <- sample_meta %>%
    rename(gsm = geo_accession) %>%
    mutate(
      disease_state = sub("^disease state: ", "", characteristics_ch1_3),
      tissue_subregion = sub("^tissue subregion: ", "", characteristics_ch1_2),
      tissue = sub("^tissue: ", "", characteristics_ch1_1),
      source_name = source_name_ch1,
      relation_clean = relation
    )

  expr_text <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
  expr_df <- read.delim(text = expr_text, check.names = FALSE, stringsAsFactors = FALSE)
  names(expr_df)[1] <- "probe_id"
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(expr_mat) <- expr_df$probe_id
  storage.mode(expr_mat) <- "numeric"
  colnames(expr_mat) <- gsub('^"|"$', "", colnames(expr_mat))

  list(meta = sample_meta, expr = expr_mat)
}

message("Parsing GSE30122 series matrix...")
parsed <- parse_series_matrix(series_path)
pheno <- parsed$meta
probe_mat <- parsed$expr

message("Loading GPL571 annotation...")
gpl <- GEOquery::getGEO("GPL571", destdir = gpl_dir)
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

message("Collapsing probes to gene symbols...")
probe_df <- as.data.frame(probe_mat) %>%
  rownames_to_column("probe_id") %>%
  inner_join(probe_map, by = "probe_id") %>%
  select(symbol, all_of(colnames(probe_mat)))

gene_expr <- probe_df %>%
  group_by(symbol) %>%
  summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

gene_mat <- as.matrix(gene_expr[, -1, drop = FALSE])
rownames(gene_mat) <- gene_expr$symbol
storage.mode(gene_mat) <- "numeric"

signature_defs <- readr::read_csv(sig_path, show_col_types = FALSE) %>%
  mutate(
    genes = as.character(genes),
    human_genes = vapply(strsplit(genes, ";"), function(gs) paste(unique(toupper(trimws(gs))), collapse = ";"), character(1))
  )

extra_signatures <- tibble::tribble(
  ~signature, ~focus_group, ~n_genes, ~genes, ~human_genes,
  "SPP1_CD44_anchor_score", "anchor", 2, "SPP1;CD44", "SPP1;CD44",
  "human_tubular_injury_context", "curated_human_context", 5, "SPP1;CD44;HAVCR1;LCN2;VCAM1", "SPP1;CD44;HAVCR1;LCN2;VCAM1"
)

signature_defs <- bind_rows(
  signature_defs %>% select(signature, focus_group, n_genes, genes, human_genes),
  extra_signatures
)

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

score_signature <- function(mat, genes) {
  genes <- unique(toupper(trimws(genes)))
  present <- intersect(genes, toupper(rownames(mat)))
  idx <- match(present, toupper(rownames(mat)))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    return(list(score = rep(NA_real_, ncol(mat)), present = character(0)))
  }
  sub <- mat[idx, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[is.na(z)] <- 0
  list(score = colMeans(z), present = rownames(mat)[idx])
}

score_rows <- list()
presence_rows <- list()
for (i in seq_len(nrow(signature_defs))) {
  sig <- signature_defs[i, ]
  genes <- strsplit(sig$human_genes, ";")[[1]]
  res <- score_signature(gene_mat, genes)
  score_rows[[sig$signature]] <- tibble(
    gsm = colnames(gene_mat),
    signature = sig$signature,
    score = as.numeric(res$score)
  )
  presence_rows[[sig$signature]] <- tibble(
    signature = sig$signature,
    requested_n = length(unique(toupper(trimws(genes)))),
    present_n = length(res$present),
    present_genes = paste(res$present, collapse = ";"),
    missing_genes = paste(setdiff(unique(toupper(trimws(genes))), toupper(res$present)), collapse = ";")
  )
}

scores <- bind_rows(score_rows) %>%
  left_join(pheno, by = "gsm") %>%
  mutate(
    disease_binary = case_when(
      disease_state == "diabetic kidney disease (DKD)" ~ "DKD",
      disease_state == "control" ~ "Control",
      TRUE ~ NA_character_
    ),
    subregion_clean = case_when(
      tissue_subregion %in% c("tubules", "tubulus") ~ "tubular",
      tissue_subregion == "glomerulus" ~ "glomerulus",
      TRUE ~ tissue_subregion
    ),
    is_primary_tubules = tissue_subregion == "tubules",
    is_tubular_any = tissue_subregion %in% c("tubules", "tubulus"),
    is_reanalysis = grepl("Reanalysis of", relation_clean, fixed = TRUE)
  )

presence <- bind_rows(presence_rows)
readr::write_csv(pheno, file.path(out_dir, "gse30122_metadata_parsed.csv"))
readr::write_csv(presence, file.path(out_dir, "gse30122_signature_gene_presence.csv"))
readr::write_csv(scores, file.path(out_dir, "gse30122_signature_scores_long.csv"))

run_group_tests <- function(df, analysis_name) {
  df %>%
    filter(signature %in% priority_sigs, disease_binary %in% c("Control", "DKD"), !is.na(score)) %>%
    group_by(signature) %>%
    summarise(
      analysis = analysis_name,
      control_n = sum(disease_binary == "Control"),
      dkd_n = sum(disease_binary == "DKD"),
      control_mean = mean(score[disease_binary == "Control"], na.rm = TRUE),
      dkd_mean = mean(score[disease_binary == "DKD"], na.rm = TRUE),
      dkd_minus_control = dkd_mean - control_mean,
      wilcox_p = suppressWarnings(wilcox.test(score ~ disease_binary, data = cur_data(), exact = FALSE)$p.value),
      cohen_d = {
        x <- score[disease_binary == "DKD"]
        y <- score[disease_binary == "Control"]
        pooled <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
        ifelse(is.finite(pooled) && pooled > 0, (mean(x) - mean(y)) / pooled, NA_real_)
      },
      .groups = "drop"
    ) %>%
    mutate(
      wilcox_fdr = p.adjust(wilcox_p, method = "BH"),
      direction = if_else(dkd_minus_control > 0, "higher_in_DKD", "lower_in_DKD"),
      signature_label = recode(signature, !!!friendly_names)
    ) %>%
    arrange(wilcox_p)
}

primary_df <- scores %>%
  filter(is_primary_tubules, disease_binary %in% c("Control", "DKD"))

all_tubular_df <- scores %>%
  filter(is_tubular_any, disease_binary %in% c("Control", "DKD"))

glomerulus_df <- scores %>%
  filter(tissue_subregion == "glomerulus", disease_binary %in% c("Control", "DKD"))

tests <- bind_rows(
  run_group_tests(primary_df, "primary_tubules_only"),
  run_group_tests(all_tubular_df, "all_tubular_including_control_reanalysis"),
  run_group_tests(glomerulus_df, "glomerulus_specificity")
)
readr::write_csv(tests, file.path(out_dir, "gse30122_signature_group_tests.csv"))

plot_df <- primary_df %>%
  filter(signature %in% priority_sigs) %>%
  mutate(
    signature_label = factor(recode(signature, !!!friendly_names), levels = friendly_names[priority_sigs]),
    disease_binary = factor(disease_binary, levels = c("Control", "DKD"))
  )

p1 <- ggplot(plot_df, aes(x = disease_binary, y = score, fill = disease_binary)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.0, alpha = 0.86) +
  facet_wrap(~signature_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Control" = "#4D4D4D", "DKD" = "#B2182B")) +
  labs(
    title = "GSE30122 human DKD tubules: signature validation",
    subtitle = "Primary tubules comparison: DKD n=10 vs control n=12",
    x = NULL,
    y = "Signature score",
    fill = NULL
  ) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = "grey80")
  )
ggsave(file.path(out_dir, "FigureCK_gse30122_human_dkd_tubules_signature_boxplots.png"), p1, width = 8.8, height = 6.8, dpi = 300)

heat_df <- tests %>%
  filter(signature %in% priority_sigs) %>%
  mutate(
    signature_label = factor(signature_label, levels = friendly_names[priority_sigs]),
    analysis = factor(
      analysis,
      levels = c("primary_tubules_only", "all_tubular_including_control_reanalysis", "glomerulus_specificity"),
      labels = c("Primary tubules", "All tubular", "Glomerulus")
    )
  )

p2 <- ggplot(heat_df, aes(x = analysis, y = signature_label, fill = dkd_minus_control)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("delta=%.2f\np=%.3g", dkd_minus_control, wilcox_p)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "GSE30122 DKD-control signature effect summary",
    subtitle = "Positive delta indicates higher score in DKD",
    x = NULL,
    y = NULL,
    fill = "DKD - Control"
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out_dir, "FigureCL_gse30122_dkd_signature_effect_heatmap.png"), p2, width = 8.4, height = 4.8, dpi = 300)

gene_focus <- c("SPP1", "CD44", "HAVCR1", "LCN2", "VCAM1")
gene_scores <- gene_mat[intersect(gene_focus, rownames(gene_mat)), , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "gsm", values_to = "expression") %>%
  left_join(pheno, by = "gsm") %>%
  mutate(
    disease_binary = case_when(
      disease_state == "diabetic kidney disease (DKD)" ~ "DKD",
      disease_state == "control" ~ "Control",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(tissue_subregion == "tubules", disease_binary %in% c("Control", "DKD"))
readr::write_csv(gene_scores, file.path(out_dir, "gse30122_focus_gene_expression_long.csv"))

gene_tests <- gene_scores %>%
  group_by(gene) %>%
  summarise(
    control_n = sum(disease_binary == "Control"),
    dkd_n = sum(disease_binary == "DKD"),
    control_mean = mean(expression[disease_binary == "Control"], na.rm = TRUE),
    dkd_mean = mean(expression[disease_binary == "DKD"], na.rm = TRUE),
    dkd_minus_control = dkd_mean - control_mean,
    wilcox_p = suppressWarnings(wilcox.test(expression ~ disease_binary, data = cur_data(), exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH")) %>%
  arrange(wilcox_p)
readr::write_csv(gene_tests, file.path(out_dir, "gse30122_focus_gene_tests.csv"))

main_primary <- tests %>% filter(analysis == "primary_tubules_only", signature == "Spp1_Cd44_tubular_immune_program")
injury_primary <- tests %>% filter(analysis == "primary_tubules_only", signature == "human_tubular_injury_context")
anchor_primary <- tests %>% filter(analysis == "primary_tubules_only", signature == "SPP1_CD44_anchor_score")
fib_primary <- tests %>% filter(analysis == "primary_tubules_only", signature == "curated_fibro_inflammatory_context")
main_all_tub <- tests %>% filter(analysis == "all_tubular_including_control_reanalysis", signature == "Spp1_Cd44_tubular_immune_program")
main_glom <- tests %>% filter(analysis == "glomerulus_specificity", signature == "Spp1_Cd44_tubular_immune_program")

sample_counts <- pheno %>%
  count(`disease_state`, `tissue_subregion`) %>%
  mutate(label = paste0(disease_state, " / ", tissue_subregion, " n=", n)) %>%
  pull(label)

lines <- c(
  "GSE30122 human DKD tubulointerstitial bulk validation",
  "",
  "Purpose:",
  "This patch addresses the main reviewer concern that GSE175759 is IgAN-predominant by adding a DKD-specific human kidney dataset.",
  "",
  "Dataset:",
  "GSE30122, human diabetic kidney disease microarray dataset, GPL571.",
  paste0("Sample composition parsed from GEO: ", paste(sample_counts, collapse = "; ")),
  "",
  "Primary analysis design:",
  "Primary comparison uses tissue_subregion == tubules only, avoiding the additional control tubulus reanalysis samples.",
  sprintf("- Primary tubules: DKD n=%s vs Control n=%s", main_primary$dkd_n, main_primary$control_n),
  "",
  "Main signature results in primary tubules:",
  sprintf(
    "- Spp1_Cd44_tubular_immune_program: DKD mean=%.3f, Control mean=%.3f, delta=%.3f, Wilcoxon p=%.4g, FDR=%.4g, direction=%s",
    main_primary$dkd_mean, main_primary$control_mean, main_primary$dkd_minus_control, main_primary$wilcox_p, main_primary$wilcox_fdr, main_primary$direction
  ),
  sprintf(
    "- human_tubular_injury_context: DKD mean=%.3f, Control mean=%.3f, delta=%.3f, Wilcoxon p=%.4g, FDR=%.4g, direction=%s",
    injury_primary$dkd_mean, injury_primary$control_mean, injury_primary$dkd_minus_control, injury_primary$wilcox_p, injury_primary$wilcox_fdr, injury_primary$direction
  ),
  sprintf(
    "- SPP1_CD44_anchor_score: DKD mean=%.3f, Control mean=%.3f, delta=%.3f, Wilcoxon p=%.4g, FDR=%.4g, direction=%s",
    anchor_primary$dkd_mean, anchor_primary$control_mean, anchor_primary$dkd_minus_control, anchor_primary$wilcox_p, anchor_primary$wilcox_fdr, anchor_primary$direction
  ),
  sprintf(
    "- curated_fibro_inflammatory_context: DKD mean=%.3f, Control mean=%.3f, delta=%.3f, Wilcoxon p=%.4g, FDR=%.4g, direction=%s",
    fib_primary$dkd_mean, fib_primary$control_mean, fib_primary$dkd_minus_control, fib_primary$wilcox_p, fib_primary$wilcox_fdr, fib_primary$direction
  ),
  "",
  "Sensitivity:",
  sprintf(
    "- All tubular including control reanalysis for main program: DKD n=%s vs Control n=%s, delta=%.3f, p=%.4g, FDR=%.4g",
    main_all_tub$dkd_n, main_all_tub$control_n, main_all_tub$dkd_minus_control, main_all_tub$wilcox_p, main_all_tub$wilcox_fdr
  ),
  sprintf(
    "- Glomerulus specificity for main program: DKD n=%s vs Control n=%s, delta=%.3f, p=%.4g, FDR=%.4g",
    main_glom$dkd_n, main_glom$control_n, main_glom$dkd_minus_control, main_glom$wilcox_p, main_glom$wilcox_fdr
  ),
  "",
  "Interpretation boundary:",
  "GSE30122 is an older microarray dataset, not RNA-seq. It is useful as a disease-matched human DKD support layer, but should be written as orthogonal external support rather than definitive clinical validation.",
  "The primary value is to repair the etiology gap left by GSE175759, which is IgAN-predominant.",
  "",
  "Recommended manuscript placement:",
  "If the main Spp1/Cd44 tubular-immune program is higher in DKD tubules, use this as a small but important human DKD-matched validation panel alongside GSE175759.",
  "If the signal is weak or inconsistent, keep it in supplementary material and avoid claiming DKD-specific validation.",
  "",
  "Output files:",
  "- gse30122_metadata_parsed.csv",
  "- gse30122_signature_gene_presence.csv",
  "- gse30122_signature_scores_long.csv",
  "- gse30122_signature_group_tests.csv",
  "- gse30122_focus_gene_expression_long.csv",
  "- gse30122_focus_gene_tests.csv",
  "- FigureCK_gse30122_human_dkd_tubules_signature_boxplots.png",
  "- FigureCL_gse30122_dkd_signature_effect_heatmap.png"
)

write_text(file.path(project_root, "analysis_summary.txt"), lines)
write_text(file.path(out_dir, "gse30122_human_dkd_bulk_validation_summary.txt"), lines)
write_text(file.path(out_dir, "status.txt"), "completed")

cat(paste(lines, collapse = "\n"))

