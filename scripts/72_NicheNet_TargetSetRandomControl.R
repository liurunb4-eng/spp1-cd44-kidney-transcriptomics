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
if (basename(project_root) == "scripts") {
  project_root <- dirname(project_root)
}

trailing <- commandArgs(trailingOnly = TRUE)
n_perm <- 5000L
seed <- 20260622L
for (arg in trailing) {
  if (grepl("^--n_perm=", arg)) {
    n_perm <- as.integer(sub("^--n_perm=", "", arg))
  }
  if (grepl("^--seed=", arg)) {
    seed <- as.integer(sub("^--seed=", "", arg))
  }
}
set.seed(seed)

out_dir <- file.path(project_root, "results", "human_nichenet_targetset_random_control")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gse30122_data_dir <- file.path(project_root, "data", "public_datasets", "DKD_bulk_candidate", "GSE30122")
meta_path <- file.path(project_root, "results", "human_dkd_bulk_gse30122_signature_validation", "gse30122_metadata_parsed.csv")
target_path <- file.path(project_root, "results", "human_nichenet_targetset_validation", "nichenet_spp1_top100_targetset_genes.csv")
observed_tests_path <- file.path(project_root, "results", "human_nichenet_targetset_validation", "gse30122_nichenet_targetset_group_tests.csv")

required <- c(
  meta_path,
  target_path,
  observed_tests_path,
  file.path(gse30122_data_dir, "GSE30122_series_matrix.txt.gz"),
  file.path(gse30122_data_dir, "GPL571.soft.gz")
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"))
}

parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
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

parse_gpl571_symbols <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) {
    stop("Cannot locate GPL table in ", path)
  }
  gpl_text <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
  gpl <- read.delim(text = gpl_text, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("ID", "Gene Symbol") %in% colnames(gpl))) {
    stop("GPL571 table does not contain required columns ID and Gene Symbol")
  }
  gpl %>%
    transmute(probe_id = ID, gene_symbol_raw = `Gene Symbol`) %>%
    filter(!is.na(gene_symbol_raw), gene_symbol_raw != "") %>%
    mutate(gene_symbol_raw = gsub(" /// ", ";", gene_symbol_raw, fixed = TRUE)) %>%
    tidyr::separate_rows(gene_symbol_raw, sep = ";") %>%
    mutate(symbol = toupper(trimws(gene_symbol_raw))) %>%
    filter(symbol != "", !grepl("///", symbol, fixed = TRUE)) %>%
    distinct(probe_id, symbol)
}

precompute_rank_context <- function(mat, alpha = 0.25) {
  n_genes <- nrow(mat)
  n_samples <- ncol(mat)
  weight <- matrix(NA_real_, nrow = n_genes, ncol = n_samples)
  weight_tail <- matrix(NA_real_, nrow = n_genes, ncol = n_samples)
  tail <- matrix(NA_real_, nrow = n_genes, ncol = n_samples)
  total_tail <- numeric(n_samples)
  for (j in seq_len(n_samples)) {
    x <- mat[, j]
    ranks <- rank(x, ties.method = "average")
    ord <- order(x, decreasing = TRUE)
    pos_tail <- rev(seq_len(n_genes))
    w <- ranks[ord]^alpha
    weight[ord, j] <- w
    weight_tail[ord, j] <- w * pos_tail
    tail[ord, j] <- pos_tail
    total_tail[j] <- sum(pos_tail)
  }
  list(weight = weight, weight_tail = weight_tail, tail = tail, total_tail = total_tail, n_genes = n_genes)
}

score_rank_fast <- function(rank_ctx, idx) {
  nh <- length(idx)
  n <- rank_ctx$n_genes
  nm <- n - nh
  if (nh == 0 || nm <= 0) {
    return(rep(NA_real_, ncol(rank_ctx$weight)))
  }
  weight_sum <- colSums(rank_ctx$weight[idx, , drop = FALSE])
  hit_weight_tail_sum <- colSums(rank_ctx$weight_tail[idx, , drop = FALSE])
  hit_tail_sum <- colSums(rank_ctx$tail[idx, , drop = FALSE])
  phit_sum <- hit_weight_tail_sum / weight_sum
  pmiss_sum <- (rank_ctx$total_tail - hit_tail_sum) / nm
  as.numeric((phit_sum - pmiss_sum) / n)
}

effect_delta <- function(scores, dkd_idx, ctrl_idx) {
  mean(scores[dkd_idx], na.rm = TRUE) - mean(scores[ctrl_idx], na.rm = TRUE)
}

effect_p <- function(scores, dkd_idx, ctrl_idx) {
  suppressWarnings(wilcox.test(scores[dkd_idx], scores[ctrl_idx], exact = FALSE)$p.value)
}

message("Loading GSE30122 expression matrix and GPL571 annotation...")
probe_mat <- parse_series_matrix(file.path(gse30122_data_dir, "GSE30122_series_matrix.txt.gz"))
probe_map <- parse_gpl571_symbols(file.path(gse30122_data_dir, "GPL571.soft.gz"))
meta301 <- readr::read_csv(meta_path, show_col_types = FALSE)

gene_expr <- as.data.frame(probe_mat) %>%
  rownames_to_column("probe_id") %>%
  inner_join(probe_map, by = "probe_id") %>%
  select(symbol, all_of(colnames(probe_mat))) %>%
  group_by(symbol) %>%
  summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

mat301 <- as.matrix(gene_expr[, -1, drop = FALSE])
rownames(mat301) <- toupper(gene_expr$symbol)
storage.mode(mat301) <- "numeric"
rn <- rownames(mat301)

zmat301 <- t(scale(t(mat301)))
zmat301[is.na(zmat301)] <- 0
rank_ctx <- precompute_rank_context(mat301)

targets <- readr::read_csv(target_path, show_col_types = FALSE) %>%
  mutate(gene = toupper(gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)
target_present <- intersect(targets$gene, rn)
target_idx <- match(target_present, rn)
target_n <- length(target_idx)
if (target_n < 10) {
  stop("Too few target genes present in GSE30122: ", target_n)
}

meta_short <- tibble(gsm = colnames(mat301)) %>%
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

sample_sets <- list(
  primary_tubules_only = list(
    dkd = which(meta_short$disease_binary == "DKD" & meta_short$is_primary_tubules),
    ctrl = which(meta_short$disease_binary == "Control" & meta_short$is_primary_tubules)
  ),
  all_tubular_including_control_reanalysis = list(
    dkd = which(meta_short$disease_binary == "DKD" & meta_short$is_tubular_any),
    ctrl = which(meta_short$disease_binary == "Control" & meta_short$is_tubular_any)
  ),
  glomerulus_specificity = list(
    dkd = which(meta_short$disease_binary == "DKD" & meta_short$tissue_subregion == "glomerulus"),
    ctrl = which(meta_short$disease_binary == "Control" & meta_short$tissue_subregion == "glomerulus")
  )
)

gene_stats <- tibble(
  gene = rn,
  row_index = seq_along(rn),
  mean_expr = rowMeans(mat301, na.rm = TRUE),
  sd_expr = apply(mat301, 1, sd, na.rm = TRUE)
) %>%
  filter(is.finite(mean_expr), is.finite(sd_expr), sd_expr > 0) %>%
  mutate(
    expr_bin = ntile(mean_expr, 20L),
    is_target = gene %in% target_present
  )

target_bins <- gene_stats %>%
  filter(is_target) %>%
  count(expr_bin, name = "n")
background <- gene_stats %>% filter(!is_target)
pool_by_bin <- split(background$row_index, background$expr_bin)

sample_matched_idx <- function() {
  sampled <- unlist(lapply(seq_len(nrow(target_bins)), function(i) {
    bin <- as.character(target_bins$expr_bin[i])
    n <- target_bins$n[i]
    pool <- pool_by_bin[[bin]]
    if (length(pool) == 0) {
      stop("No background genes available for expression bin ", bin)
    }
    sample(pool, size = n, replace = length(pool) < n)
  }), use.names = FALSE)
  sampled <- unique(sampled)
  if (length(sampled) < target_n) {
    extra_pool <- setdiff(background$row_index, sampled)
    sampled <- c(sampled, sample(extra_pool, target_n - length(sampled), replace = FALSE))
  }
  sampled
}

score_types <- c("mean_z_module_score", "rank_ssgsea_like_score")
analyses <- names(sample_sets)

observed_mean <- colMeans(zmat301[target_idx, , drop = FALSE])
observed_rank <- score_rank_fast(rank_ctx, target_idx)
observed_calc <- bind_rows(lapply(score_types, function(score_type) {
  scores <- if (score_type == "mean_z_module_score") observed_mean else observed_rank
  bind_rows(lapply(analyses, function(a) {
    idx <- sample_sets[[a]]
    tibble(
      score_type = score_type,
      analysis = a,
      observed_delta_recalculated = effect_delta(scores, idx$dkd, idx$ctrl),
      observed_wilcox_p_recalculated = effect_p(scores, idx$dkd, idx$ctrl)
    )
  }))
}))

observed_reported <- readr::read_csv(observed_tests_path, show_col_types = FALSE) %>%
  filter(analysis %in% analyses) %>%
  select(score_type, analysis, observed_delta = dkd_minus_control, observed_wilcox_p = wilcox_p)

message("Running fast matched random gene-set control: n_perm=", n_perm, ", target_present_n=", target_n)
perm_rows <- vector("list", n_perm * length(score_types) * length(analyses))
k <- 0L
for (perm_id in seq_len(n_perm)) {
  idx_genes <- sample_matched_idx()
  mean_scores <- colMeans(zmat301[idx_genes, , drop = FALSE])
  rank_scores <- score_rank_fast(rank_ctx, idx_genes)
  for (score_type in score_types) {
    scores <- if (score_type == "mean_z_module_score") mean_scores else rank_scores
    for (a in analyses) {
      idx <- sample_sets[[a]]
      k <- k + 1L
      perm_rows[[k]] <- data.frame(
        perm_id = perm_id,
        score_type = score_type,
        analysis = a,
        random_delta = effect_delta(scores, idx$dkd, idx$ctrl),
        gene_set_size = length(idx_genes),
        stringsAsFactors = FALSE
      )
    }
  }
  if (perm_id %% 1000L == 0L) {
    message("  completed ", perm_id, "/", n_perm)
  }
}
perm <- bind_rows(perm_rows)
readr::write_csv(perm, file.path(out_dir, "gse30122_nichenet_targetset_matched_random_permutations.csv.gz"))

summary <- perm %>%
  left_join(observed_reported, by = c("score_type", "analysis")) %>%
  group_by(score_type, analysis) %>%
  summarise(
    n_perm = n(),
    target_present_n = target_n,
    observed_delta = first(observed_delta),
    observed_wilcox_p = first(observed_wilcox_p),
    random_delta_mean = mean(random_delta, na.rm = TRUE),
    random_delta_sd = sd(random_delta, na.rm = TRUE),
    random_delta_q025 = quantile(random_delta, 0.025, na.rm = TRUE),
    random_delta_q500 = quantile(random_delta, 0.500, na.rm = TRUE),
    random_delta_q975 = quantile(random_delta, 0.975, na.rm = TRUE),
    empirical_p_right = (sum(random_delta >= first(observed_delta), na.rm = TRUE) + 1) / (sum(is.finite(random_delta)) + 1),
    empirical_p_two_sided_abs = (sum(abs(random_delta) >= abs(first(observed_delta)), na.rm = TRUE) + 1) / (sum(is.finite(random_delta)) + 1),
    .groups = "drop"
  ) %>%
  left_join(observed_calc, by = c("score_type", "analysis"))
readr::write_csv(summary, file.path(out_dir, "gse30122_nichenet_targetset_matched_random_summary.csv"))
readr::write_csv(target_bins, file.path(out_dir, "target_gene_expression_bin_counts.csv"))

plot_df <- perm %>%
  filter(analysis %in% c("primary_tubules_only", "glomerulus_specificity")) %>%
  left_join(observed_reported, by = c("score_type", "analysis")) %>%
  mutate(
    score_label = recode(
      score_type,
      mean_z_module_score = "Mean z-score",
      rank_ssgsea_like_score = "Rank-based ssGSEA-like"
    ),
    analysis_label = recode(
      analysis,
      primary_tubules_only = "Primary tubules",
      glomerulus_specificity = "Glomerulus"
    ),
    score_label = factor(score_label, levels = c("Mean z-score", "Rank-based ssGSEA-like")),
    analysis_label = factor(analysis_label, levels = c("Primary tubules", "Glomerulus"))
  )

p <- ggplot(plot_df, aes(x = random_delta)) +
  geom_histogram(bins = 60, fill = "#D9D9D9", color = "white", linewidth = 0.2) +
  geom_vline(aes(xintercept = observed_delta), color = "#B2182B", linewidth = 0.75) +
  facet_grid(score_label ~ analysis_label, scales = "free_x") +
  labs(
    title = "Matched random gene-set control for the NicheNet SPP1 target set",
    subtitle = paste0("Observed DKD-Control delta shown in red; n=", n_perm, " expression-bin-matched random gene sets"),
    x = "Random gene-set DKD - Control delta",
    y = "Number of random sets"
  )
ggsave(file.path(out_dir, "Supplementary_Figure_S4_NicheNet_targetset_matched_random_control.png"), p, width = 9.2, height = 5.8, dpi = 300)
ggsave(file.path(out_dir, "Supplementary_Figure_S4_NicheNet_targetset_matched_random_control.svg"), p, width = 9.2, height = 5.8)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg"))
}
primary <- summary %>% filter(analysis == "primary_tubules_only")
glom <- summary %>% filter(analysis == "glomerulus_specificity")
summary_lines <- c(
  "Matched random gene-set control for the NicheNet SPP1 top target set in GSE30122",
  paste0("Date: ", Sys.Date()),
  paste0("Random seed: ", seed),
  paste0("Matched random gene sets: ", n_perm),
  paste0("Present NicheNet target genes in GSE30122: ", target_n, "/", nrow(targets)),
  "",
  "Primary tubules:",
  paste0(
    "- Mean z-score observed delta=", fmt(primary$observed_delta[primary$score_type == "mean_z_module_score"]),
    "; empirical right-tail P=", fmt(primary$empirical_p_right[primary$score_type == "mean_z_module_score"]),
    "; empirical two-sided |delta| P=", fmt(primary$empirical_p_two_sided_abs[primary$score_type == "mean_z_module_score"])
  ),
  paste0(
    "- Rank-based ssGSEA-like observed delta=", fmt(primary$observed_delta[primary$score_type == "rank_ssgsea_like_score"]),
    "; empirical right-tail P=", fmt(primary$empirical_p_right[primary$score_type == "rank_ssgsea_like_score"]),
    "; empirical two-sided |delta| P=", fmt(primary$empirical_p_two_sided_abs[primary$score_type == "rank_ssgsea_like_score"])
  ),
  "",
  "Glomerulus specificity context:",
  paste0(
    "- Mean z-score observed delta=", fmt(glom$observed_delta[glom$score_type == "mean_z_module_score"]),
    "; empirical right-tail P=", fmt(glom$empirical_p_right[glom$score_type == "mean_z_module_score"])
  ),
  paste0(
    "- Rank-based ssGSEA-like observed delta=", fmt(glom$observed_delta[glom$score_type == "rank_ssgsea_like_score"]),
    "; empirical right-tail P=", fmt(glom$empirical_p_right[glom$score_type == "rank_ssgsea_like_score"])
  ),
  "",
  "Interpretation guardrail:",
  "- Random sets were matched to the detected GSE30122 background by target-gene mean-expression bins.",
  "- This control asks whether the observed target-set DKD-Control effect is larger than expected for expression-matched gene sets of the same size.",
  "- It does not validate causal SPP1/CD44 signaling and should be reported only as a robustness check for the receiver-target projection."
)
writeLines(enc2utf8(summary_lines), con = file.path(out_dir, "nichenet_targetset_matched_random_control_summary.txt"), useBytes = TRUE)

writeLines(
  enc2utf8(c(
    "nichenet_targetset_matched_random_control_fast completed",
    paste0("timestamp: ", Sys.time()),
    paste0("n_perm: ", n_perm),
    paste0("seed: ", seed)
  )),
  con = file.path(out_dir, "status.txt"),
  useBytes = TRUE
)
message("Done. Outputs written to: ", out_dir)
