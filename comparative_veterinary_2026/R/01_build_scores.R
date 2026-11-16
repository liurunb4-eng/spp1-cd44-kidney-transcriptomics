try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
source(file.path(dirname(normalizePath(.this, winslash = "/")), "00_config.R"))
assert_data_contracts()

presence_rows <- list()
score_rows <- list()

record_presence <- function(dataset, species, genes_present) {
  tibble(
    dataset = dataset,
    species = species,
    gene = GENE_SETS$full10,
    present = GENE_SETS$full10 %in% toupper(genes_present)
  )
}

add_scores <- function(expr, dataset, species, meta, log_transform = FALSE, add_rank = FALSE) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- toupper(rownames(expr))
  if (log_transform) expr <- log2(expr + 1)
  presence_rows[[dataset]] <<- record_presence(dataset, species, rownames(expr))
  z <- score_gene_sets(expr, dataset, species, meta, "z")
  if (add_rank) {
    r <- score_gene_sets(expr, dataset, species, meta, "rank")
    score_rows[[dataset]] <<- bind_rows(z, r)
  } else {
    score_rows[[dataset]] <<- z
  }
}

load_gse175759 <- function() {
  message("Loading GSE175759...")
  raw_counts <- readr::read_csv(PATHS$gse175759_counts, show_col_types = FALSE)
  sample_columns <- setdiff(names(raw_counts), "ensembl_version")
  library_sizes <- colSums(as.matrix(raw_counts[, sample_columns, drop = FALSE]), na.rm = TRUE)
  map <- readr::read_csv(PATHS$gse175759_map, show_col_types = FALSE) |>
    mutate(ensembl = as.character(ensembl), symbol = toupper(symbol)) |>
    filter(symbol %in% GENE_SETS$full10)
  counts <- raw_counts |>
    mutate(ensembl = sub("\\..*$", "", ensembl_version)) |>
    inner_join(map, by = "ensembl") |>
    select(-ensembl_version, -ensembl) |>
    relocate(symbol) |>
    group_by(symbol) |>
    summarise(across(everything(), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  mat <- as.matrix(counts[, -1, drop = FALSE])
  rownames(mat) <- counts$symbol
  log_cpm <- log2(sweep(mat, 2, library_sizes[colnames(mat)], "/") * 1e6 + 1)
  meta <- readr::read_csv(PATHS$gse175759_meta, show_col_types = FALSE) |>
    transmute(
      sample_id = gsm,
      group = diagnosis,
      stage_num = NA_real_,
      continuous_outcome = egfr_ckd_epi,
      outcome_name = "eGFR",
      technical_outlier = technical_outlier,
      tissue = tissue,
      role = "human_continuous"
    ) |>
    filter(technical_outlier == "No", is.finite(continuous_outcome))
  keep <- intersect(colnames(log_cpm), meta$sample_id)
  add_scores(log_cpm[, keep, drop = FALSE], "GSE175759_eGFR", "Human", meta, FALSE, FALSE)
}

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) stop("Invalid GEO series matrix: ", path)
  sample_rows <- lines[grep("^!Sample_", lines)]
  row_names <- sub("\\t.*$", "", sample_rows)
  row_names <- sub("^!Sample_", "", row_names)
  counts <- ave(seq_along(row_names), row_names, FUN = seq_along)
  unique_names <- ifelse(duplicated(row_names) | duplicated(row_names, fromLast = TRUE), paste0(row_names, "_", counts), row_names)
  sm <- NULL
  for (i in seq_along(sample_rows)) {
    parts <- strsplit(sample_rows[[i]], "\t", fixed = TRUE)[[1]]
    values <- gsub('^"|"$', "", parts[-1])
    if (is.null(sm)) sm <- tibble(sample_index = seq_along(values))
    sm[[unique_names[[i]]]] <- values
  }
  sm <- sm |>
    rename(sample_id = geo_accession) |>
    mutate(
      disease_state = sub("^disease state: ", "", characteristics_ch1_3),
      tissue_subregion = sub("^tissue subregion: ", "", characteristics_ch1_2)
    )
  txt <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
  df <- read.delim(text = txt, check.names = FALSE, stringsAsFactors = FALSE)
  names(df)[1] <- "probe_id"
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df$probe_id
  storage.mode(mat) <- "numeric"
  colnames(mat) <- gsub('^"|"$', "", colnames(mat))
  list(meta = sm, expr = mat)
}

load_gse30122 <- function() {
  message("Loading GSE30122...")
  parsed <- parse_series_matrix(PATHS$gse30122_series)
  suppressPackageStartupMessages(library(GEOquery))
  gpl <- GEOquery::getGEO("GPL571", destdir = PATHS$gse30122_gpl_dir)
  ann <- GEOquery::Table(gpl) |>
    as_tibble() |>
    transmute(probe_id = ID, raw = `Gene Symbol`) |>
    filter(!is.na(raw), raw != "") |>
    mutate(raw = gsub(" /// ", ";", raw, fixed = TRUE)) |>
    separate_rows(raw, sep = ";") |>
    transmute(probe_id, symbol = toupper(trimws(raw))) |>
    filter(symbol %in% GENE_SETS$full10) |>
    distinct(probe_id, symbol)
  df <- as.data.frame(parsed$expr) |>
    tibble::rownames_to_column("probe_id") |>
    inner_join(ann, by = "probe_id") |>
    select(symbol, all_of(colnames(parsed$expr))) |>
    group_by(symbol) |>
    summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  mat <- as.matrix(df[, -1, drop = FALSE]); rownames(mat) <- df$symbol
  meta <- parsed$meta |>
    transmute(
      sample_id,
      group = case_when(
        disease_state == "diabetic kidney disease (DKD)" ~ "DKD",
        disease_state == "control" ~ "Control",
        TRUE ~ NA_character_
      ),
      stage_num = if_else(group == "DKD", 1, 0),
      continuous_outcome = NA_real_, outcome_name = "DKD status",
      tissue = tissue_subregion, role = "human_disease"
    ) |>
    filter(tissue == "tubules", group %in% c("Control", "DKD"))
  keep <- intersect(colnames(mat), meta$sample_id)
  add_scores(mat[, keep, drop = FALSE], "GSE30122_DKD_tubules", "Human", meta, FALSE, FALSE)
}

load_gse137570 <- function() {
  message("Loading GSE137570...")
  x <- readr::read_csv(PATHS$gse137570_components, show_col_types = FALSE)
  genes <- x |>
    filter(grepl("^gene_", score_name)) |>
    transmute(sample_id = sample, gene = sub("^gene_", "", score_name), score) |>
    pivot_wider(names_from = sample_id, values_from = score)
  mat <- as.matrix(genes[, -1, drop = FALSE])
  rownames(mat) <- genes$gene
  meta <- x |>
    distinct(sample_id = sample, progression_status, CKD_Progression) |>
    transmute(
      sample_id,
      group = progression_status,
      stage_num = CKD_Progression,
      continuous_outcome = CKD_Progression,
      outcome_name = "progression status",
      tissue = "kidney biopsy",
      role = "human_progression"
    )
  add_scores(mat, "GSE137570_progression", "Human", meta, FALSE, FALSE)
}

read_symbol_matrix <- function(path) {
  d <- readr::read_tsv(path, show_col_types = FALSE)
  names(d)[1] <- "gene"
  d |>
    mutate(gene = toupper(trimws(as.character(gene)))) |>
    filter(gene %in% GENE_SETS$full10) |>
    group_by(gene) |>
    summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
}

read_ensembl_rat <- function(path) {
  d <- readr::read_tsv(path, show_col_types = FALSE)
  names(d)[1] <- "ensembl"
  mp <- readr::read_csv(PATHS$rat_map, show_col_types = FALSE) |>
    transmute(symbol = toupper(symbol), ensembl = as.character(rat_ensembl_id)) |>
    filter(symbol %in% GENE_SETS$full10)
  d |>
    mutate(ensembl = sub("\\..*$", "", as.character(ensembl))) |>
    inner_join(mp, by = "ensembl") |>
    select(-ensembl) |>
    relocate(symbol) |>
    group_by(symbol) |>
    summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop") |>
    rename(gene = symbol)
}

add_rat_arm <- function(df, dataset, group_map, reference, contrast, role) {
  mat <- as.matrix(df[, -1, drop = FALSE]); rownames(mat) <- df$gene
  keep <- intersect(colnames(mat), names(group_map))
  keep <- keep[group_map[keep] %in% c(reference, contrast)]
  meta <- tibble(
    sample_id = keep,
    group = unname(group_map[keep]),
    stage_num = if_else(group == contrast, 1, 0),
    continuous_outcome = NA_real_, outcome_name = paste(contrast, "vs", reference),
    tissue = "whole kidney", role = role
  )
  add_scores(mat[, keep, drop = FALSE], dataset, "Rat", meta, TRUE, FALSE)
}

load_rat <- function() {
  message("Loading rat datasets...")
  ad <- read_symbol_matrix(PATHS$rat_gse216376_adenine)
  uu <- read_symbol_matrix(PATHS$rat_gse216376_uuo)
  gd_ad <- setNames(ifelse(grepl("control", names(ad)[-1], ignore.case = TRUE), "control", "adenine"), names(ad)[-1])
  gd_uu <- setNames(ifelse(grepl("sham", names(uu)[-1], ignore.case = TRUE), "sham", "uuo"), names(uu)[-1])
  add_rat_arm(ad, "GSE216376_adenine", gd_ad, "control", "adenine", "rat_disease")
  add_rat_arm(uu, "GSE216376_UUO", gd_uu, "sham", "uuo", "rat_disease")
  dc <- read_ensembl_rat(PATHS$rat_gse183841)
  md <- readr::read_csv(PATHS$rat_meta, show_col_types = FALSE)
  gd <- setNames(md$normalized_group, md$sample_id)
  add_rat_arm(dc, "GSE183841_DOCA", gd, "control", "doca", "rat_disease")
  add_rat_arm(dc, "GSE183841_finerenone", gd, "doca", "finerenone", "rat_treatment")
  add_rat_arm(dc, "GSE183841_spironolactone", gd, "doca", "spironolactone", "rat_treatment")
}

read_feline_rna <- function(sheet, tissue) {
  raw <- readxl::read_excel(PATHS$feline_xlsx, sheet = sheet, col_names = FALSE, .name_repair = "minimal")
  sample_id <- as.character(unlist(raw[2, -1], use.names = FALSE))
  meta <- tibble(
    sample_id = sample_id,
    group = as.character(unlist(raw[3, -1], use.names = FALSE)),
    sex = as.character(unlist(raw[4, -1], use.names = FALSE)),
    age = as.numeric(unlist(raw[5, -1], use.names = FALSE)),
    breed = as.character(unlist(raw[6, -1], use.names = FALSE))
  ) |>
    filter(!is.na(sample_id), nzchar(sample_id), !is.na(group)) |>
    mutate(
      stage_num = match(group, c("Control", "CKD1/2", "CKD3/4")) - 1,
      continuous_outcome = stage_num,
      outcome_name = "ordinal CKD group",
      tissue = tissue,
      role = "feline_natural_disease"
    )
  dat <- as.data.frame(raw[-(1:6), , drop = FALSE], stringsAsFactors = FALSE)
  genes <- toupper(as.character(dat[[1]]))
  mat <- apply(dat[, -1, drop = FALSE], 2, as.numeric)
  rownames(mat) <- genes; colnames(mat) <- sample_id
  mat <- mat[, meta$sample_id, drop = FALSE]
  keep <- !is.na(rownames(mat)) & nzchar(rownames(mat)) & rowSums(is.finite(mat)) > 0
  mat <- mat[keep, , drop = FALSE]
  if (anyDuplicated(rownames(mat))) {
    ug <- unique(rownames(mat))
    mat <- do.call(rbind, lapply(ug, function(g) colMeans(mat[rownames(mat) == g, , drop = FALSE], na.rm = TRUE)))
    rownames(mat) <- ug
  }
  list(expr = mat, meta = meta)
}

load_feline <- function() {
  message("Loading feline GSE303653 processed matrices...")
  cx <- read_feline_rna("S3 RNA cortex", "Cortex")
  md <- read_feline_rna("S4 RNA medulla", "Medulla")
  # The supplementary workbook already contains normalized log-scale expression,
  # including negative centered values, so no second log transformation is valid.
  add_scores(cx$expr, "GSE303653_feline_cortex", "Cat", cx$meta, FALSE, TRUE)
  add_scores(md$expr, "GSE303653_feline_medulla", "Cat", md$meta, FALSE, TRUE)
}

load_gse175759()
load_gse30122()
load_gse137570()
load_rat()
load_feline()

scores <- bind_rows(score_rows) |>
  relocate(dataset, species, sample_id, score_name, score, transform)
presence <- bind_rows(presence_rows)

write_csv_stable(scores, file.path(OUT_DIR, "harmonized_sample_scores.csv"))
write_csv_stable(presence, file.path(OUT_DIR, "gene_presence_by_dataset.csv"))

message("Wrote ", nrow(scores), " score rows from ", n_distinct(scores$dataset), " datasets/arms.")
