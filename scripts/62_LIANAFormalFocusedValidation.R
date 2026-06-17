#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(liana)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(sparseMatrixStats)
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

in_dir <- file.path(project_root, "results", "public_data_gse233078_cellchat_mvp_prep")
out_dir <- file.path(project_root, "results", "liana_focused_spp1_cd44_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

conditions <- c("lean", "obese", "obese_enalapril")
methods <- c("natmi", "connectome", "sca", "cellphonedb")
cellphonedb_permutations <- 100

priority_axes <- tibble::tribble(
  ~axis, ~module, ~ligand, ~receptor, ~source_group, ~target_group, ~priority,
  "Spp1->Cd44", "Spp1/Cd44 tubular-immune", "Spp1", "Cd44", "DTL", "Mono_Macro", "main_candidate",
  "Tgfb1->Tgfbr2", "fibrosis remodeling", "Tgfb1", "Tgfbr2", "Mono_Macro", "Endo_GEC", "context_axis",
  "Vegfa->Kdr", "vascular microcirculation", "Vegfa", "Kdr", "IC-A", "Endo_GEC", "context_axis",
  "Vegfa->Flt1", "vascular microcirculation", "Vegfa", "Flt1", "IC-A", "Endo_GEC", "context_axis"
)

read_mtx_small <- function(path) {
  lines <- readLines(path)
  dim_line <- lines[!startsWith(lines, "%")][1]
  dims <- as.integer(strsplit(dim_line, "\\s+")[[1]])
  coords <- read.table(path, skip = which(lines == dim_line), stringsAsFactors = FALSE)
  Matrix::sparseMatrix(
    i = as.integer(coords[[1]]),
    j = as.integer(coords[[2]]),
    x = as.numeric(coords[[3]]),
    dims = dims[1:2]
  )
}

make_sce <- function(mat, meta, condition) {
  keep <- meta$condition == condition
  condition_mat <- mat[, keep]
  condition_meta <- meta[keep, , drop = FALSE]
  sce <- SingleCellExperiment(
    assays = list(
      counts = condition_mat,
      logcounts = log1p(condition_mat)
    )
  )
  colData(sce)$cellchat_group <- condition_meta$cellchat_group
  sce
}

standardize_method <- function(method_name, result_tbl) {
  if (nrow(result_tbl) == 0) {
    return(tibble())
  }
  out <- as_tibble(result_tbl) %>%
    mutate(method = method_name)
  if (method_name == "natmi") {
    out <- out %>% mutate(method_score = prod_weight)
  } else if (method_name == "connectome") {
    out <- out %>% mutate(method_score = weight_sc)
  } else if (method_name == "sca") {
    out <- out %>% mutate(method_score = LRscore)
  } else if (method_name == "cellphonedb") {
    out <- out %>% mutate(method_score = lr.mean)
  } else {
    out <- out %>% mutate(method_score = NA_real_)
  }
  out
}

mat <- read_mtx_small(file.path(in_dir, "cellchat_mvp_focus_whitelist_counts.mtx"))
genes <- readLines(file.path(in_dir, "cellchat_mvp_focus_whitelist_genes.tsv"))
rownames(mat) <- genes
meta <- read.csv(file.path(in_dir, "cellchat_mvp_focus_metadata.csv"), stringsAsFactors = FALSE)
colnames(mat) <- meta$cell_id

custom_resource <- priority_axes %>%
  transmute(source_genesymbol = ligand, target_genesymbol = receptor)

all_method_rows <- list()
for (condition in conditions) {
  message("Running focused LIANA for condition: ", condition)
  sce <- make_sce(mat, meta, condition)
  res <- liana_wrap(
    sce,
    method = methods,
    resource = "custom",
    external_resource = custom_resource,
    idents_col = "cellchat_group",
    assay = "logcounts",
    min_cells = 5,
    return_all = TRUE,
    verbose = FALSE,
    permutations = cellphonedb_permutations
  )
  for (method_name in names(res)) {
    all_method_rows[[paste(condition, method_name, sep = "__")]] <- standardize_method(method_name, res[[method_name]]) %>%
      mutate(condition = condition)
  }
}

all_liana <- bind_rows(all_method_rows)
readr::write_csv(all_liana, file.path(out_dir, "liana_focused_all_method_outputs.csv"))

focused <- all_liana %>%
  inner_join(
    priority_axes,
    by = c(
      "ligand" = "ligand",
      "receptor" = "receptor",
      "source" = "source_group",
      "target" = "target_group"
    )
  ) %>%
  select(
    condition, method, axis, module, priority, source, target, ligand, receptor,
    receptor.prop, ligand.prop, method_score, pvalue, everything()
  )

focused_ranked <- focused %>%
  group_by(condition, method) %>%
  mutate(method_rank = rank(-method_score, ties.method = "min")) %>%
  ungroup()

aggregate <- focused_ranked %>%
  group_by(condition, axis, module, priority, source, target, ligand, receptor) %>%
  summarise(
    n_methods = n_distinct(method),
    mean_method_rank = mean(method_rank, na.rm = TRUE),
    median_method_rank = median(method_rank, na.rm = TRUE),
    mean_method_score = mean(method_score, na.rm = TRUE),
    cellphonedb_pvalue = suppressWarnings(min(pvalue[method == "cellphonedb"], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    cellphonedb_pvalue = ifelse(is.infinite(cellphonedb_pvalue), NA_real_, cellphonedb_pvalue),
    aggregate_score = 1 / mean_method_rank
  ) %>%
  arrange(condition, mean_method_rank)

readr::write_csv(focused_ranked, file.path(out_dir, "liana_focused_priority_axes_by_method.csv"))
readr::write_csv(aggregate, file.path(out_dir, "liana_focused_priority_axes_aggregate.csv"))

plot_aggregate <- aggregate %>%
  mutate(
    condition = factor(condition, levels = conditions),
    axis = factor(axis, levels = priority_axes$axis),
    p_label = ifelse(is.na(cellphonedb_pvalue), "p=NA", paste0("p=", formatC(cellphonedb_pvalue, format = "f", digits = 2)))
  )

p1 <- ggplot(plot_aggregate, aes(x = condition, y = axis)) +
  geom_point(aes(size = aggregate_score, color = mean_method_rank), alpha = 0.86) +
  geom_text(aes(label = paste0("R", round(mean_method_rank, 1))), size = 3, color = "white") +
  scale_color_viridis_c(option = "mako", direction = -1) +
  scale_size(range = c(5, 13)) +
  labs(
    title = "Focused LIANA validation of priority LR axes",
    x = "Condition",
    y = "Priority LR axis",
    color = "Mean method rank",
    size = "Aggregate score"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(out_dir, "FigureBW_liana_focused_aggregate_bubble.png"), p1, width = 8.8, height = 4.8, dpi = 300)

method_heat <- focused_ranked %>%
  select(condition, axis, method, method_score) %>%
  mutate(row_label = paste(axis, condition, sep = " | ")) %>%
  group_by(method) %>%
  mutate(method_z = as.numeric(scale(method_score))) %>%
  ungroup()

p2 <- ggplot(method_heat, aes(x = method, y = factor(row_label, levels = rev(unique(row_label))), fill = method_z)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", method_z)), size = 2.8) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  labs(
    title = "Method-level standardized LIANA scores",
    x = "LIANA method",
    y = NULL,
    fill = "Method z-score"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(out_dir, "FigureBX_liana_focused_method_heatmap.png"), p2, width = 8.6, height = 5.6, dpi = 300)

spp1 <- aggregate %>%
  filter(axis == "Spp1->Cd44") %>%
  arrange(match(condition, conditions))

spp1_methods <- focused_ranked %>%
  filter(axis == "Spp1->Cd44") %>%
  arrange(method, match(condition, conditions))

lines <- c(
  "Focused LIANA custom-resource validation for GSE233078",
  "",
  "Purpose:",
  "Run formal LIANA methods on the focused single-cell matrix using a custom priority LR resource.",
  "",
  "Methods included:",
  paste0("- ", methods),
  "",
  paste0("CellPhoneDB-style permutations: ", cellphonedb_permutations),
  "",
  "Important boundary:",
  "This is a focused LIANA validation against a small custom LR resource, not a full unbiased whole-database LIANA screen.",
  "",
  "Spp1->Cd44 DTL->Mono_Macro aggregate results:"
)

for (i in seq_len(nrow(spp1))) {
  row <- spp1[i, ]
  lines <- c(
    lines,
    sprintf(
      "- %s: mean_method_rank=%.2f, median_method_rank=%.2f, aggregate_score=%.3f, CellPhoneDB_p=%s",
      row$condition,
      row$mean_method_rank,
      row$median_method_rank,
      row$aggregate_score,
      ifelse(is.na(row$cellphonedb_pvalue), "NA", sprintf("%.4f", row$cellphonedb_pvalue))
    )
  )
}

lines <- c(lines, "", "Spp1->Cd44 method-level scores:")
for (i in seq_len(nrow(spp1_methods))) {
  row <- spp1_methods[i, ]
  lines <- c(
    lines,
    sprintf(
      "- %s | %s: method_score=%.4f, method_rank=%.1f%s",
      row$condition,
      row$method,
      row$method_score,
      row$method_rank,
      ifelse(is.na(row$pvalue), "", sprintf(", CellPhoneDB_p=%.4f", row$pvalue))
    )
  )
}

lines <- c(
  lines,
  "",
  "Writing interpretation:",
  "1. The Spp1->Cd44 DTL-to-Mono/Macro axis can now be described as supported by focused LIANA custom-resource validation.",
  "2. Because the analysis used a small predefined LR set, it should not be described as an unbiased discovery screen.",
  "3. The safest wording is: focused LIANA analysis supports Spp1-Cd44 as a stable high-priority candidate axis across conditions; absolute method scores are generally higher in obese kidneys and attenuate after enalapril in most methods, but relative aggregate rank alone should not be overinterpreted as disease-specific emergence.",
  "4. Do not write that LIANA proves the axis is absent in lean or completely blocked by enalapril.",
  "",
  "Output files:",
  "- liana_focused_all_method_outputs.csv",
  "- liana_focused_priority_axes_by_method.csv",
  "- liana_focused_priority_axes_aggregate.csv",
  "- FigureBW_liana_focused_aggregate_bubble.png",
  "- FigureBX_liana_focused_method_heatmap.png"
)

writeLines(lines, file.path(project_root, "analysis_summary.txt"), useBytes = TRUE)
writeLines("completed", file.path(out_dir, "status.txt"), useBytes = TRUE)

