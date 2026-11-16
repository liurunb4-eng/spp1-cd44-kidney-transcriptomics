options(stringsAsFactors = FALSE)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
})

SEED <- 20261113L
set.seed(SEED)

locate_this_file <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) {
    if (is.null(x$ofile)) "" else as.character(x$ofile)
  }, character(1))
  hit <- tail(ofiles[nzchar(ofiles)], 1)
  if (!length(hit)) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

CONFIG_FILE <- locate_this_file()
PACKAGE_DIR <- normalizePath(file.path(dirname(CONFIG_FILE), ".."), winslash = "/", mustWork = TRUE)

normalize_optional <- function(path) {
  if (!nzchar(path)) return("")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

env_or_default <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  normalize_optional(if (nzchar(value)) value else default)
}

# The release never assumes the author's local directory hierarchy. Users may
# either populate data/input/ as described in docs/input_data_manifest.md or
# point KIDNEY_PROJECT_DIR at the legacy project root used by earlier releases.
INPUT_DIR <- env_or_default(
  "KIDNEY_COMPARATIVE_INPUT_DIR",
  file.path(PACKAGE_DIR, "data", "input")
)
PROJECT_DIR <- normalize_optional(trimws(Sys.getenv("KIDNEY_PROJECT_DIR", unset = "")))

OUT_DIR <- file.path(PACKAGE_DIR, "results")
FIG_DIR <- file.path(PACKAGE_DIR, "figures")
MANUSCRIPT_DIR <- file.path(PACKAGE_DIR, "manuscript")
TEST_DIR <- file.path(PACKAGE_DIR, "tests")
for (p in c(OUT_DIR, FIG_DIR, MANUSCRIPT_DIR, TEST_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

GENE_SETS <- list(
  full10 = c("HAVCR1", "LCN2", "VCAM1", "CLU", "COL1A1", "FN1", "ACTA2", "TGFB1", "TLR4", "NLRP3"),
  common9 = c("HAVCR1", "VCAM1", "CLU", "COL1A1", "FN1", "ACTA2", "TGFB1", "TLR4", "NLRP3"),
  injury = c("HAVCR1", "LCN2", "VCAM1", "CLU", "TLR4", "NLRP3"),
  remodeling = c("COL1A1", "FN1", "ACTA2", "TGFB1")
)

if (nzchar(PROJECT_DIR)) {
  default_paths <- list(
    gse175759_counts = file.path(PROJECT_DIR, "results", "human_bulk_gse175759_clinical_correlation", "gse175759_merged_counts_by_ensembl.csv.gz"),
    gse175759_map = file.path(PROJECT_DIR, "results", "human_bulk_gse175759_clinical_correlation", "gse175759_ensembl_to_symbol_map.csv"),
    gse175759_meta = file.path(PROJECT_DIR, "results", "human_bulk_gse175759_clinical_correlation", "gse175759_metadata_clean.csv"),
    gse30122_series = file.path(PROJECT_DIR, "data", "public_datasets", "DKD_bulk_candidate", "GSE30122", "GSE30122_series_matrix.txt.gz"),
    gse30122_gpl_dir = file.path(PROJECT_DIR, "data", "public_datasets", "DKD_bulk_candidate", "GSE30122"),
    gse137570_components = file.path(PROJECT_DIR, "manuscript_drafts", "Diabetologia_relaunch_2026-10-04", "35_gse137570_ckd_progression_screen", "GSE137570_cohort2_component_gene_scores_long.csv"),
    rat_gse216376_adenine = file.path(PROJECT_DIR, "data", "public_datasets", "GSE216376", "processed", "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz"),
    rat_gse216376_uuo = file.path(PROJECT_DIR, "data", "public_datasets", "GSE216376", "processed", "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz"),
    rat_gse183841 = file.path(PROJECT_DIR, "data", "public_datasets", "GSE183841", "processed", "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz"),
    rat_map = file.path(PROJECT_DIR, "results", "bulk_signature_projection_spp1_cd44_program", "rat_symbol_to_ensembl_mapping.csv"),
    rat_meta = file.path(PROJECT_DIR, "results", "public_data_gse183841", "sample_metadata.csv"),
    feline_xlsx = file.path(INPUT_DIR, "GSE303653", "42003_2025_9164_MOESM3_ESM.xlsx")
  )
} else {
  default_paths <- list(
    gse175759_counts = file.path(INPUT_DIR, "GSE175759", "gse175759_merged_counts_by_ensembl.csv.gz"),
    gse175759_map = file.path(PACKAGE_DIR, "data", "derived", "GSE175759", "gse175759_ensembl_to_symbol_map.csv"),
    gse175759_meta = file.path(PACKAGE_DIR, "data", "derived", "GSE175759", "gse175759_metadata_clean.csv"),
    gse30122_series = file.path(INPUT_DIR, "GSE30122", "GSE30122_series_matrix.txt.gz"),
    gse30122_gpl_dir = file.path(INPUT_DIR, "GSE30122"),
    gse137570_components = file.path(PACKAGE_DIR, "data", "derived", "GSE137570", "GSE137570_cohort2_component_gene_scores_long.csv"),
    rat_gse216376_adenine = file.path(INPUT_DIR, "GSE216376", "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz"),
    rat_gse216376_uuo = file.path(INPUT_DIR, "GSE216376", "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz"),
    rat_gse183841 = file.path(INPUT_DIR, "GSE183841", "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz"),
    rat_map = file.path(PACKAGE_DIR, "data", "derived", "mappings", "rat_symbol_to_ensembl_mapping.csv"),
    rat_meta = file.path(PACKAGE_DIR, "data", "derived", "GSE183841", "sample_metadata.csv"),
    feline_xlsx = file.path(INPUT_DIR, "GSE303653", "42003_2025_9164_MOESM3_ESM.xlsx")
  )
}

default_paths$gse137570_components <- env_or_default(
  "GSE137570_COMPONENTS_FILE", default_paths$gse137570_components
)
default_paths$feline_xlsx <- env_or_default(
  "FELINE_CKD_XLSX", default_paths$feline_xlsx
)
PATHS <- lapply(default_paths, normalize_optional)

required_paths <- function() unlist(PATHS, use.names = TRUE)

assert_data_contracts <- function() {
  missing <- required_paths()[!file.exists(required_paths())]
  if (length(missing)) {
    stop(
      "Missing required input files:\n", paste(missing, collapse = "\n"),
      "\n\nSee docs/input_data_manifest.md. You may set KIDNEY_COMPARATIVE_INPUT_DIR, ",
      "KIDNEY_PROJECT_DIR, GSE137570_COMPONENTS_FILE, or FELINE_CKD_XLSX."
    )
  }
  stopifnot(length(GENE_SETS$full10) == 10L)
  stopifnot(length(GENE_SETS$common9) == 9L)
  stopifnot(identical(setdiff(GENE_SETS$full10, GENE_SETS$common9), "LCN2"))
  stopifnot(setequal(c(GENE_SETS$injury, GENE_SETS$remodeling), GENE_SETS$full10))
  stopifnot(startsWith(normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE), PACKAGE_DIR))
  invisible(TRUE)
}

gene_z <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mu <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- NA_real_
  z <- sweep(sweep(mat, 1, mu, "-"), 1, s, "/")
  z[!is.finite(z)] <- 0
  z
}

rank_fraction <- function(mat) {
  mat <- as.matrix(mat)
  apply(mat, 2, function(x) rank(x, ties.method = "average", na.last = "keep") / sum(is.finite(x)))
}

score_gene_sets <- function(expr, dataset, species, meta, transform = c("z", "rank")) {
  transform <- match.arg(transform)
  rownames(expr) <- toupper(rownames(expr))
  score_mat <- if (transform == "z") gene_z(expr) else rank_fraction(expr)
  defs <- c(
    list(full_available = intersect(GENE_SETS$full10, rownames(score_mat))),
    list(common9 = intersect(GENE_SETS$common9, rownames(score_mat))),
    list(injury = intersect(GENE_SETS$injury, rownames(score_mat))),
    list(remodeling = intersect(GENE_SETS$remodeling, rownames(score_mat)))
  )
  for (g in intersect(GENE_SETS$full10, rownames(score_mat))) defs[[paste0("gene_", g)]] <- g
  bind_rows(lapply(names(defs), function(nm) {
    gs <- defs[[nm]]
    if (!length(gs)) return(NULL)
    tibble(
      dataset = dataset,
      species = species,
      sample_id = colnames(score_mat),
      score_name = nm,
      score = colMeans(score_mat[gs, , drop = FALSE], na.rm = TRUE),
      transform = transform,
      n_genes = length(gs),
      genes = paste(gs, collapse = ";")
    )
  })) |>
    left_join(meta, by = "sample_id")
}

hedges_g <- function(reference, contrast) {
  reference <- reference[is.finite(reference)]
  contrast <- contrast[is.finite(contrast)]
  n1 <- length(reference); n2 <- length(contrast)
  if (n1 < 2 || n2 < 2) return(c(g = NA_real_, var_g = NA_real_))
  sp2 <- ((n1 - 1) * var(reference) + (n2 - 1) * var(contrast)) / (n1 + n2 - 2)
  if (!is.finite(sp2) || sp2 <= 0) return(c(g = NA_real_, var_g = NA_real_))
  d <- (mean(contrast) - mean(reference)) / sqrt(sp2)
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  vg <- (n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2 - 2))
  c(g = g, var_g = vg)
}

cliffs_delta <- function(reference, contrast) {
  reference <- reference[is.finite(reference)]
  contrast <- contrast[is.finite(contrast)]
  if (!length(reference) || !length(contrast)) return(NA_real_)
  cmp <- outer(contrast, reference, "-")
  (sum(cmp > 0) - sum(cmp < 0)) / length(cmp)
}

bootstrap_stat <- function(df, fun, B = 2000L, seed = SEED) {
  set.seed(seed)
  out <- replicate(B, {
    idx <- unlist(lapply(split(seq_len(nrow(df)), df$group), function(ii) sample(ii, length(ii), replace = TRUE)))
    fun(df[idx, , drop = FALSE])
  })
  quantile(out[is.finite(out)], c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
}

write_csv_stable <- function(x, path) readr::write_csv(x, path, na = "")
