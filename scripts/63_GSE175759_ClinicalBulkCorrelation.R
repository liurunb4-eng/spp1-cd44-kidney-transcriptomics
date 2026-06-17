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

data_dir <- file.path(project_root, "data", "public_datasets", "Human_CKD_bulk_clinical_screen", "GSE175759")
out_dir <- file.path(project_root, "results", "human_bulk_gse175759_clinical_correlation")
sig_path <- file.path(project_root, "results", "bulk_signature_projection_spp1_cd44_program", "signature_definitions.csv")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_text <- function(path, lines) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

message("Loading GEO metadata...")
gse <- GEOquery::getGEO("GSE175759", GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
gse <- gse[[1]]
pheno <- pData(gse) %>%
  rownames_to_column("gsm") %>%
  as_tibble()

pheno_clean <- pheno %>%
  transmute(
    gsm = geo_accession,
    title = title,
    diagnosis = `diagnosis:ch1`,
    egfr_ckd_epi = suppressWarnings(as.numeric(`estimated gfr (ckd-epi):ch1`)),
    technical_outlier = `technical outlier:ch1`,
    tissue = `tissue:ch1`,
    supplementary_file = supplementary_file_1
  )

readr::write_csv(pheno_clean, file.path(out_dir, "gse175759_metadata_clean.csv"))

download_one <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(invisible(TRUE))
  }
  https_url <- sub("^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov", url)
  download.file(https_url, destfile = dest, mode = "wb", quiet = TRUE, timeout = 600)
  invisible(TRUE)
}

message("Downloading sample count files if missing...")
for (i in seq_len(nrow(pheno_clean))) {
  url <- pheno_clean$supplementary_file[i]
  gsm <- pheno_clean$gsm[i]
  sample_label <- sub(".*_(sample[0-9]+)\\.txt\\.gz$", "\\1", basename(url))
  dest <- file.path(data_dir, paste0(gsm, "_", sample_label, ".txt.gz"))
  tryCatch(download_one(url, dest), error = function(e) {
    message("Download failed for ", gsm, ": ", conditionMessage(e))
  })
}

message("Merging count matrix...")
count_files <- file.path(
  data_dir,
  paste0(
    pheno_clean$gsm,
    "_",
    sub(".*_(sample[0-9]+)\\.txt\\.gz$", "\\1", basename(pheno_clean$supplementary_file)),
    ".txt.gz"
  )
)
names(count_files) <- pheno_clean$gsm
ok_files <- file.exists(count_files) & file.info(count_files)$size > 0
if (sum(ok_files) < 50) {
  stop("Too few sample count files were downloaded: ", sum(ok_files))
}

read_count <- function(path, sample_id) {
  readr::read_tsv(path, show_col_types = FALSE, col_names = TRUE) %>%
    setNames(c("ensembl_version", sample_id))
}

count_list <- Map(read_count, count_files[ok_files], names(count_files)[ok_files])
counts <- Reduce(function(x, y) full_join(x, y, by = "ensembl_version"), count_list)
counts[is.na(counts)] <- 0
readr::write_csv(counts, file.path(out_dir, "gse175759_merged_counts_by_ensembl.csv.gz"))

ensembl_clean <- sub("\\..*$", "", counts$ensembl_version)
symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensembl_clean),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
gene_map <- tibble(
  ensembl = names(symbols),
  symbol = as.character(symbols)
) %>%
  filter(!is.na(symbol), symbol != "")
readr::write_csv(gene_map, file.path(out_dir, "gse175759_ensembl_to_symbol_map.csv"))

expr <- counts %>%
  mutate(ensembl = ensembl_clean) %>%
  left_join(gene_map, by = "ensembl") %>%
  filter(!is.na(symbol), symbol != "") %>%
  select(symbol, all_of(names(count_files)[ok_files])) %>%
  group_by(symbol) %>%
  summarise(across(everything(), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

expr_mat <- as.matrix(expr[, -1, drop = FALSE])
rownames(expr_mat) <- expr$symbol
mode(expr_mat) <- "numeric"

lib_size <- colSums(expr_mat)
cpm <- t(t(expr_mat) / lib_size * 1e6)
log_cpm <- log2(cpm + 1)

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

score_signature <- function(mat, genes) {
  present <- intersect(toupper(genes), toupper(rownames(mat)))
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
  genes <- unique(toupper(trimws(strsplit(sig$human_genes, ";")[[1]])))
  res <- score_signature(log_cpm, genes)
  score_rows[[sig$signature]] <- tibble(
    gsm = colnames(log_cpm),
    signature = sig$signature,
    score = as.numeric(res$score)
  )
  presence_rows[[sig$signature]] <- tibble(
    signature = sig$signature,
    requested_n = length(genes),
    present_n = length(res$present),
    present_genes = paste(res$present, collapse = ";"),
    missing_genes = paste(setdiff(genes, toupper(res$present)), collapse = ";")
  )
}

scores <- bind_rows(score_rows) %>%
  left_join(pheno_clean, by = "gsm") %>%
  mutate(
    usable = !is.na(egfr_ckd_epi) & technical_outlier == "No"
  )
presence <- bind_rows(presence_rows)
readr::write_csv(scores, file.path(out_dir, "gse175759_signature_scores_long.csv"))
readr::write_csv(presence, file.path(out_dir, "gse175759_signature_gene_presence.csv"))

cor_rows <- scores %>%
  filter(usable, !is.na(score)) %>%
  group_by(signature) %>%
  summarise(
    n = n(),
    spearman_rho = suppressWarnings(cor(score, egfr_ckd_epi, method = "spearman")),
    spearman_p = suppressWarnings(cor.test(score, egfr_ckd_epi, method = "spearman", exact = FALSE)$p.value),
    pearson_r = suppressWarnings(cor(score, egfr_ckd_epi, method = "pearson")),
    pearson_p = suppressWarnings(cor.test(score, egfr_ckd_epi, method = "pearson")$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    spearman_fdr = p.adjust(spearman_p, method = "BH"),
    pearson_fdr = p.adjust(pearson_p, method = "BH")
  ) %>%
  arrange(spearman_p)

readr::write_csv(cor_rows, file.path(out_dir, "gse175759_signature_egfr_correlations.csv"))

priority_sigs <- c(
  "Spp1_Cd44_tubular_immune_program",
  "SPP1_CD44_anchor_score",
  "human_tubular_injury_context",
  "curated_fibro_inflammatory_context"
)

plot_df <- scores %>%
  filter(usable, signature %in% priority_sigs, !is.na(score)) %>%
  mutate(
    signature = factor(signature, levels = priority_sigs),
    diagnosis = factor(diagnosis)
  )

p1 <- ggplot(plot_df, aes(x = score, y = egfr_ckd_epi, color = diagnosis)) +
  geom_point(size = 2.1, alpha = 0.82) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
  facet_wrap(~signature, scales = "free_x", ncol = 2) +
  labs(
    title = "GSE175759 human kidney tubulointerstitium: signature score versus eGFR",
    x = "Signature score (mean gene-level z-score)",
    y = "eGFR (CKD-EPI)",
    color = "Diagnosis"
  )

ggsave(file.path(out_dir, "FigureBY_gse175759_signature_vs_egfr_scatter.png"), p1, width = 9.2, height = 7.2, dpi = 300)

heat_df <- cor_rows %>%
  filter(signature %in% priority_sigs) %>%
  mutate(signature = factor(signature, levels = priority_sigs))

p2 <- ggplot(heat_df, aes(x = "eGFR", y = signature, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("rho=%.2f\np=%.3g", spearman_rho, spearman_p)), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0) +
  labs(
    title = "Clinical correlation summary in GSE175759",
    x = NULL,
    y = NULL,
    fill = "Spearman rho"
  )

ggsave(file.path(out_dir, "FigureBZ_gse175759_egfr_correlation_heatmap.png"), p2, width = 6.6, height = 3.9, dpi = 300)

score_wide <- scores %>%
  filter(usable, signature %in% priority_sigs) %>%
  select(gsm, diagnosis, egfr_ckd_epi, signature, score) %>%
  pivot_wider(names_from = signature, values_from = score)
readr::write_csv(score_wide, file.path(out_dir, "gse175759_priority_signature_scores_wide.csv"))

main_sig <- cor_rows %>% filter(signature == "Spp1_Cd44_tubular_immune_program")
anchor_sig <- cor_rows %>% filter(signature == "SPP1_CD44_anchor_score")
injury_sig <- cor_rows %>% filter(signature == "human_tubular_injury_context")
usable_n <- scores %>% filter(usable) %>% distinct(gsm) %>% nrow()
diagnosis_counts <- pheno_clean %>%
  filter(technical_outlier == "No") %>%
  count(diagnosis) %>%
  mutate(label = paste0(diagnosis, " n=", n)) %>%
  pull(label)

lines <- c(
  "GSE175759 human clinical bulk correlation analysis",
  "",
  "Dataset:",
  "GSE175759 contains human kidney tubulointerstitium RNA-seq samples with estimated GFR (CKD-EPI) metadata.",
  paste0("Usable non-technical-outlier samples with eGFR: n=", usable_n),
  paste0("Diagnosis composition: ", paste(diagnosis_counts, collapse = "; ")),
  "",
  "Key signature correlation with eGFR:",
  sprintf(
    "- Spp1_Cd44_tubular_immune_program: Spearman rho=%.3f, p=%.4g, FDR=%.4g; Pearson r=%.3f, p=%.4g",
    main_sig$spearman_rho, main_sig$spearman_p, main_sig$spearman_fdr, main_sig$pearson_r, main_sig$pearson_p
  ),
  sprintf(
    "- SPP1_CD44_anchor_score: Spearman rho=%.3f, p=%.4g, FDR=%.4g; Pearson r=%.3f, p=%.4g",
    anchor_sig$spearman_rho, anchor_sig$spearman_p, anchor_sig$spearman_fdr, anchor_sig$pearson_r, anchor_sig$pearson_p
  ),
  sprintf(
    "- human_tubular_injury_context: Spearman rho=%.3f, p=%.4g, FDR=%.4g; Pearson r=%.3f, p=%.4g",
    injury_sig$spearman_rho, injury_sig$spearman_p, injury_sig$spearman_fdr, injury_sig$pearson_r, injury_sig$pearson_p
  ),
  "",
  "Interpretation boundary:",
  "This is a human clinical bulk correlation layer. It can support clinical relevance of the transcriptomic program, but it cannot prove cell-cell communication or causality.",
  "Because GSE175759 is predominantly IgAN tubulointerstitium rather than DKD-only, it should be described as human CKD clinical-extension evidence, not as a DKD-specific validation.",
  "",
  "Output files:",
  "- gse175759_signature_scores_long.csv",
  "- gse175759_signature_egfr_correlations.csv",
  "- gse175759_signature_gene_presence.csv",
  "- FigureBY_gse175759_signature_vs_egfr_scatter.png",
  "- FigureBZ_gse175759_egfr_correlation_heatmap.png"
)

write_text(file.path(project_root, "analysis_summary.txt"), lines)
write_text(file.path(out_dir, "status.txt"), "completed")

