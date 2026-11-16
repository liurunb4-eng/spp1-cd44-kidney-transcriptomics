try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
source(file.path(dirname(normalizePath(.this, winslash = "/")), "00_config.R"))

set.seed(SEED)
scores <- readr::read_csv(file.path(OUT_DIR, "harmonized_sample_scores.csv"), show_col_types = FALSE)

binary_specs <- tribble(
  ~dataset, ~reference, ~contrast, ~orientation, ~endpoint_class, ~primary,
  "GSE30122_DKD_tubules", "Control", "DKD", 1, "human disease", TRUE,
  "GSE137570_progression", "non-progressive", "progressive", 1, "human progression", TRUE,
  "GSE216376_adenine", "control", "adenine", 1, "rat disease", TRUE,
  "GSE216376_UUO", "sham", "uuo", 1, "rat disease", TRUE,
  "GSE183841_DOCA", "control", "doca", 1, "rat disease", TRUE,
  "GSE183841_finerenone", "doca", "finerenone", -1, "rat intervention", TRUE,
  "GSE183841_spironolactone", "doca", "spironolactone", -1, "rat intervention", TRUE,
  "GSE303653_feline_cortex", "Control", "CKD3/4", 1, "feline late CKD", FALSE,
  "GSE303653_feline_medulla", "Control", "CKD3/4", 1, "feline late CKD", FALSE
)

continuous_specs <- tribble(
  ~dataset, ~outcome, ~orientation, ~endpoint_class, ~primary,
  "GSE175759_eGFR", "eGFR", -1, "human kidney function", TRUE,
  "GSE303653_feline_cortex", "ordinal CKD group", 1, "feline ordinal CKD", TRUE,
  "GSE303653_feline_medulla", "ordinal CKD group", 1, "feline ordinal CKD", TRUE
)

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(quantile(x, p, na.rm = TRUE, names = FALSE, type = 7))
}

boot_binary_delta <- function(reference, contrast, B = 1000L) {
  replicate(B, cliffs_delta(
    sample(reference, length(reference), replace = TRUE),
    sample(contrast, length(contrast), replace = TRUE)
  ))
}

boot_spearman <- function(score, outcome, B = 1000L) {
  n <- length(score)
  replicate(B, {
    ii <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(score[ii], outcome[ii], method = "spearman", use = "complete.obs"))
  })
}

analyse_binary <- function(d, reference, contrast, orientation, endpoint_class, primary) {
  x <- d$score[d$group == reference]
  y <- d$score[d$group == contrast]
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  delta <- cliffs_delta(x, y)
  bg <- boot_binary_delta(x, y)
  hg <- hedges_g(x, y)
  tibble(
    endpoint_type = "binary",
    endpoint_class = endpoint_class,
    primary_endpoint = primary,
    reference = reference,
    contrast = contrast,
    n_reference = length(x),
    n_contrast = length(y),
    raw_effect = delta,
    aligned_effect = orientation * delta,
    ci_low = safe_quantile(bg, 0.025),
    ci_high = safe_quantile(bg, 0.975),
    p_value = suppressWarnings(wilcox.test(y, x, exact = FALSE)$p.value),
    hedges_g = unname(hg[["g"]]),
    hedges_g_variance = unname(hg[["var_g"]])
  )
}

analyse_continuous <- function(d, orientation, endpoint_class, primary) {
  keep <- is.finite(d$score) & is.finite(d$continuous_outcome)
  s <- d$score[keep]; o <- d$continuous_outcome[keep]
  rho <- suppressWarnings(cor(s, o, method = "spearman"))
  br <- boot_spearman(s, o)
  tibble(
    endpoint_type = "continuous",
    endpoint_class = endpoint_class,
    primary_endpoint = primary,
    reference = NA_character_,
    contrast = unique(d$outcome_name)[1],
    n_reference = length(s),
    n_contrast = NA_integer_,
    raw_effect = rho,
    aligned_effect = orientation * rho,
    ci_low = safe_quantile(br, 0.025),
    ci_high = safe_quantile(br, 0.975),
    p_value = suppressWarnings(cor.test(s, o, method = "spearman", exact = FALSE)$p.value),
    hedges_g = NA_real_,
    hedges_g_variance = NA_real_
  )
}

effect_rows <- list()
row_id <- 0L
for (i in seq_len(nrow(binary_specs))) {
  sp <- binary_specs[i, ]
  ds <- scores |> filter(dataset == sp$dataset)
  for (key in unique(paste(ds$score_name, ds$transform, sep = "||"))) {
    bits <- strsplit(key, "\\|\\|")[[1]]
    d <- ds |> filter(score_name == bits[1], transform == bits[2])
    if (!all(c(sp$reference, sp$contrast) %in% d$group)) next
    row_id <- row_id + 1L
    effect_rows[[row_id]] <- analyse_binary(
      d, sp$reference, sp$contrast, sp$orientation, sp$endpoint_class, sp$primary
    ) |>
      mutate(dataset = sp$dataset, species = unique(d$species)[1], score_name = bits[1], transform = bits[2], .before = 1)
  }
}
for (i in seq_len(nrow(continuous_specs))) {
  sp <- continuous_specs[i, ]
  ds <- scores |> filter(dataset == sp$dataset)
  for (key in unique(paste(ds$score_name, ds$transform, sep = "||"))) {
    bits <- strsplit(key, "\\|\\|")[[1]]
    d <- ds |> filter(score_name == bits[1], transform == bits[2])
    row_id <- row_id + 1L
    effect_rows[[row_id]] <- analyse_continuous(d, sp$orientation, sp$endpoint_class, sp$primary) |>
      mutate(dataset = sp$dataset, species = unique(d$species)[1], score_name = bits[1], transform = bits[2], .before = 1)
  }
}

effects <- bind_rows(effect_rows)
write_csv_stable(effects, file.path(OUT_DIR, "benchmark_endpoint_effects.csv"))

consistency <- effects |>
  filter(transform == "z", primary_endpoint) |>
  mutate(
    feature_class = case_when(
      score_name %in% c("full_available", "common9") ~ "composite",
      score_name %in% c("injury", "remodeling") ~ "module",
      grepl("^gene_", score_name) ~ "single gene",
      TRUE ~ "other"
    )
  ) |>
  group_by(score_name, feature_class) |>
  summarise(
    n_endpoints = n(),
    median_aligned_effect = median(aligned_effect, na.rm = TRUE),
    q1 = quantile(aligned_effect, 0.25, na.rm = TRUE),
    q3 = quantile(aligned_effect, 0.75, na.rm = TRUE),
    positive_fraction = mean(aligned_effect > 0, na.rm = TRUE),
    nominal_fraction = mean(p_value < 0.05 & aligned_effect > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(median_aligned_effect), desc(positive_fraction))
write_csv_stable(consistency, file.path(OUT_DIR, "feature_consistency_summary.csv"))

# Feline scaling analysis. Re-standardizing each globally standardized gene on
# controls is algebraically equivalent to control-referenced standardization of
# the original expression, while retaining the exact samples used above.
feline_base <- scores |>
  filter(grepl("GSE303653_feline", dataset), score_name == "common9") |>
  transmute(dataset, species, sample_id, group, stage_num, scale = if_else(transform == "z", "global_z", "within_sample_rank"), score)

feline_gene_z <- scores |>
  filter(grepl("GSE303653_feline", dataset), transform == "z", grepl("^gene_", score_name)) |>
  mutate(gene = sub("^gene_", "", score_name)) |>
  filter(gene %in% GENE_SETS$common9) |>
  group_by(dataset, gene) |>
  mutate(
    control_mean = mean(score[group == "Control"], na.rm = TRUE),
    control_sd = sd(score[group == "Control"], na.rm = TRUE),
    control_z = if_else(is.finite(control_sd) & control_sd > 0, (score - control_mean) / control_sd, NA_real_)
  ) |>
  ungroup()

feline_control_ref <- feline_gene_z |>
  group_by(dataset, species, sample_id, group, stage_num) |>
  summarise(score = mean(control_z, na.rm = TRUE), n_genes = sum(is.finite(control_z)), .groups = "drop") |>
  mutate(scale = "control_referenced_z") |>
  select(dataset, species, sample_id, group, stage_num, scale, score)

feline_scales <- bind_rows(feline_base, feline_control_ref)
write_csv_stable(feline_scales, file.path(OUT_DIR, "feline_score_scaling_values.csv"))

feline_metric <- function(d) {
  rho <- suppressWarnings(cor(d$score, d$stage_num, method = "spearman", use = "complete.obs"))
  p <- suppressWarnings(cor.test(d$score, d$stage_num, method = "spearman", exact = FALSE)$p.value)
  delta <- cliffs_delta(d$score[d$group == "Control"], d$score[d$group == "CKD3/4"])
  tibble(n = nrow(d), ordinal_rho = rho, ordinal_p = p, late_vs_control_delta = delta)
}

feline_sensitivity <- feline_scales |>
  group_by(dataset, scale) |>
  group_modify(~{
    qc_id <- if (.y$dataset == "GSE303653_feline_cortex") "Cat_09" else "Cat_20"
    bind_rows(
      feline_metric(.x) |> mutate(analysis = "all_processed_samples"),
      feline_metric(.x |> filter(sample_id != qc_id)) |>
        mutate(analysis = "exclude_published_QC_flags")
    )
  }) |>
  ungroup()
write_csv_stable(feline_sensitivity, file.path(OUT_DIR, "feline_scaling_qc_sensitivity.csv"))

loo_rows <- list()
k <- 0L
for (key in unique(paste(feline_scales$dataset, feline_scales$scale, sep = "||"))) {
  bits <- strsplit(key, "\\|\\|")[[1]]
  d <- feline_scales |> filter(dataset == bits[1], scale == bits[2])
  for (sid in d$sample_id) {
    k <- k + 1L
    dd <- d |> filter(sample_id != sid)
    loo_rows[[k]] <- tibble(
      dataset = bits[1], scale = bits[2], omitted_sample = sid,
      omitted_group = d$group[d$sample_id == sid][1],
      ordinal_rho = suppressWarnings(cor(dd$score, dd$stage_num, method = "spearman", use = "complete.obs"))
    )
  }
}
feline_loo <- bind_rows(loo_rows)
feline_loo_summary <- feline_loo |>
  group_by(dataset, scale) |>
  summarise(min_rho = min(ordinal_rho, na.rm = TRUE), median_rho = median(ordinal_rho, na.rm = TRUE),
            max_rho = max(ordinal_rho, na.rm = TRUE), .groups = "drop")
write_csv_stable(feline_loo, file.path(OUT_DIR, "feline_leave_one_out.csv"))
write_csv_stable(feline_loo_summary, file.path(OUT_DIR, "feline_leave_one_out_summary.csv"))

message("Wrote benchmark effects for ", n_distinct(effects$dataset),
        " endpoints and feline robustness analyses across ", n_distinct(feline_scales$scale), " scales.")
