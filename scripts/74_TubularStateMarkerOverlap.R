# GPB defensive supplement:
# Cross-species tubular injury / failed-repair marker-panel overlap.
#
# Scope:
# - This is a bounded marker-panel overlap check, not a new discovery branch.
# - It tests whether the pre-specified rat DTL-associated injured-tubular
#   marker panel is enriched for published human/AKI-CKD maladaptive tubular
#   marker panels.
# - It should be reported as supplementary state-conservation evidence only.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
  project_dir <- dirname(script_path)
} else {
  project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "gpb_tubular_state_marker_overlap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_safe <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

split_genes <- function(x) {
  unique(trimws(unlist(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE))))
}

to_human_symbol <- function(x) {
  unique(toupper(trimws(x[nzchar(trimws(x))])))
}

fisher_overlap <- function(query_genes, reference_genes, universe_genes,
                           query_name, reference_name, query_note, reference_note) {
  universe <- unique(universe_genes)
  query <- intersect(to_human_symbol(query_genes), universe)
  reference <- intersect(to_human_symbol(reference_genes), universe)
  overlap <- intersect(query, reference)

  k <- length(overlap)
  a <- length(query)
  b <- length(reference)
  n <- length(universe)
  neither <- n - a - b + k
  if (neither < 0) stop("Invalid overlap table for ", query_name, " vs ", reference_name)

  mat <- matrix(c(k, a - k, b - k, neither), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "greater")

  data.frame(
    query_set = query_name,
    reference_set = reference_name,
    query_n_in_universe = a,
    reference_n_in_universe = b,
    overlap_n = k,
    universe_n = n,
    overlap_fraction_query = ifelse(a > 0, k / a, NA_real_),
    fisher_odds_ratio = unname(ft$estimate),
    fisher_p = ft$p.value,
    overlap_genes = paste(sort(overlap), collapse = ";"),
    query_note = query_note,
    reference_note = reference_note,
    stringsAsFactors = FALSE
  )
}

marker_presence_path <- file.path(
  project_dir,
  "results", "public_data_gse233078_dtl_annotation_audit",
  "dtl_annotation_marker_gene_presence.csv"
)
signature_path <- file.path(
  project_dir,
  "results", "bulk_signature_projection_spp1_cd44_program",
  "signature_definitions.csv"
)
universe_path <- file.path(
  project_dir,
  "data", "public_datasets", "GSE233078", "processed",
  "GSE233078_detected_genes_universe.txt"
)

marker_presence <- read_csv_safe(marker_presence_path)
signature_defs <- read_csv_safe(signature_path)
if (!file.exists(universe_path)) {
  stop("Gene universe file missing. Expected: ", universe_path,
       "\nCreate it from rownames of GSE233078_EXPORT_GEO_counts_postfilter.rds.gz.")
}
universe_genes <- to_human_symbol(readLines(universe_path, warn = FALSE))

injured_pt_genes <- marker_presence$gene[
  marker_presence$marker_set == "injured_pt" & marker_presence$present
]
dtl_candidate_genes <- marker_presence$gene[
  marker_presence$marker_set == "dtl_thin_limb_candidate" & marker_presence$present
]
spp1_anchor <- marker_presence$gene[
  marker_presence$marker_set == "spp1_cd44_anchor" &
    toupper(marker_presence$gene) == "SPP1" &
    marker_presence$present
]

dtl_signature_genes <- split_genes(signature_defs$genes[
  signature_defs$signature == "DTL_obese_up_enalapril_down_signature"
])
composite_genes <- split_genes(signature_defs$genes[
  signature_defs$signature == "Spp1_Cd44_tubular_immune_program"
])

query_sets <- list(
  Rat_DTL_associated_injured_tubular_anchor_panel = unique(c(injured_pt_genes, spp1_anchor)),
  Rat_DTL_candidate_marker_panel = dtl_candidate_genes,
  Rat_DTL_obese_up_enalapril_down_signature = dtl_signature_genes,
  Rat_SPP1_CD44_composite_program = composite_genes
)
query_notes <- c(
  Rat_DTL_associated_injured_tubular_anchor_panel =
    "Pre-specified GSE233078 DTL-annotation audit injured-tubular markers plus Spp1 anchor; used as state-conservation check.",
  Rat_DTL_candidate_marker_panel =
    "Pre-specified DTL candidate markers from the GSE233078 annotation audit; included as a negative/specificity contrast.",
  Rat_DTL_obese_up_enalapril_down_signature =
    "Rat DTL-annotated obese-up/enalapril-down disease-reversal signature from the existing manuscript analysis.",
  Rat_SPP1_CD44_composite_program =
    "Composite rat tubular-immune Spp1/Cd44 program; included for transparency but not used as the primary tubular-state overlap claim."
)

reference_sets <- list(
  Kirita_FR_PTC_core = c("VCAM1", "SEMA5A", "DCDC2", "CCL2"),
  Human_inflammatory_PT_NatCommun2025 = c(
    "HAVCR1", "VCAM1", "TPM1", "VIM", "SPP1", "ITGB6", "ITGB8",
    "CCL2", "CXCL1", "ICAM1", "MMP7"
  ),
  Conservative_maladaptive_tubular_injury_core = c(
    "HAVCR1", "LCN2", "VCAM1", "KRT8", "KRT18", "SOX9",
    "CLU", "SPP1", "VIM", "PROM1", "DCDC2"
  )
)
reference_notes <- c(
  Kirita_FR_PTC_core =
    "Failed-repair proximal tubule core markers reported in the AKI repair literature; compact panel for specificity testing.",
  Human_inflammatory_PT_NatCommun2025 =
    "Human kidney disease inflammatory PT markers reported for HAVCR1+/VCAM1+ PT cells localized to fibrotic niches.",
  Conservative_maladaptive_tubular_injury_core =
    "Conservative cross-study tubular injury/maladaptive marker panel assembled from commonly used failed-repair/injured PT markers; partly overlaps by design with the annotation audit and should be treated as a sanity check rather than independent validation."
)

rows <- list()
for (q in names(query_sets)) {
  for (r in names(reference_sets)) {
    rows[[paste(q, r, sep = "__")]] <- fisher_overlap(
      query_genes = query_sets[[q]],
      reference_genes = reference_sets[[r]],
      universe_genes = universe_genes,
      query_name = q,
      reference_name = r,
      query_note = query_notes[[q]],
      reference_note = reference_notes[[r]]
    )
  }
}
overlap_df <- do.call(rbind, rows)
overlap_df$fdr <- p.adjust(overlap_df$fisher_p, method = "BH")
overlap_df$minus_log10_fdr <- -log10(pmax(overlap_df$fdr, .Machine$double.xmin))
overlap_df <- overlap_df[order(overlap_df$fdr, -overlap_df$overlap_fraction_query), ]

write.csv(
  overlap_df,
  file.path(out_dir, "gpb_tubular_state_marker_overlap_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    set_type = c(rep("query", length(query_sets)), rep("reference", length(reference_sets))),
    set_name = c(names(query_sets), names(reference_sets)),
    n_genes = c(vapply(query_sets, function(x) length(to_human_symbol(x)), integer(1)),
                vapply(reference_sets, function(x) length(to_human_symbol(x)), integer(1))),
    genes = c(vapply(query_sets, function(x) paste(to_human_symbol(x), collapse = ";"), character(1)),
              vapply(reference_sets, function(x) paste(to_human_symbol(x), collapse = ";"), character(1))),
    note = c(query_notes[names(query_sets)], reference_notes[names(reference_sets)]),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "gpb_tubular_state_marker_sets.csv"),
  row.names = FALSE
)

top_primary <- overlap_df[
  overlap_df$query_set == "Rat_DTL_associated_injured_tubular_anchor_panel",
]
top_primary <- top_primary[order(top_primary$fdr), ]

summary_lines <- c(
  "GPB tubular-state marker overlap supplement",
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Detected-gene universe from GSE233078: ", length(universe_genes), " genes"),
  "",
  "Primary interpretation:",
  "- This analysis is a bounded marker-panel overlap check. It is not a new network analysis and should not be used to claim causal validation.",
  "- The primary query set is the pre-specified DTL-associated injured-tubular audit panel plus the Spp1 anchor.",
  "- DTL candidate markers and the disease-reversal signature are retained as specificity/transparency contrasts.",
  "- The conservative maladaptive tubular panel partly shares canonical injury markers with the annotation audit by design; the cleaner external signal is the overlap with compact failed-repair/inflammatory PT panels.",
  "",
  "Primary overlap rows:",
  apply(top_primary, 1, function(x) {
    paste0(
      "- ", x[["reference_set"]],
      ": overlap ", x[["overlap_n"]], "/", x[["query_n_in_universe"]],
      "; Fisher P=", signif(as.numeric(x[["fisher_p"]]), 3),
      "; BH-FDR=", signif(as.numeric(x[["fdr"]]), 3),
      "; genes=", x[["overlap_genes"]]
    )
  }),
  "",
  "Suggested Results text:",
  paste(
    "As a narrowly bounded cross-species state-conservation check, we compared the",
    "pre-specified rat DTL-associated injured-tubular marker panel with published",
    "failed-repair or inflammatory proximal-tubular marker panels. Using the detected",
    "GSE233078 gene universe as background, the rat injured-tubular anchor panel showed",
    "expected overlap with a conservative maladaptive tubular injury panel and, more",
    "informatively, shared HAVCR1/SPP1/VCAM1 with a human inflammatory PT panel and",
    "VCAM1 with a compact failed-repair PT core. These results support the cautious",
    "interpretation that the Spp1-high rat DTL-annotated compartment reflects an injured",
    "tubular state with cross-study marker concordance rather than a purely homeostatic",
    "DTL identity."
  ),
  "",
  "Suggested limitation sentence:",
  paste(
    "This marker-panel overlap is intentionally restricted to anchor-state concordance;",
    "it does not establish a full one-to-one equivalence between rat DTL-annotated cells",
    "and human failed-repair proximal tubules, nor does it validate ligand-receptor causality."
  )
)
writeLines(summary_lines, file.path(out_dir, "gpb_tubular_state_marker_overlap_summary.txt"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  plot_df <- overlap_df
  plot_df$query_set <- factor(plot_df$query_set, levels = rev(unique(plot_df$query_set)))
  plot_df$reference_set <- factor(plot_df$reference_set, levels = unique(plot_df$reference_set))

  p <- ggplot(plot_df, aes(x = reference_set, y = query_set)) +
    geom_point(aes(size = overlap_fraction_query, color = minus_log10_fdr)) +
    scale_color_gradient(low = "#8a8f98", high = "#c23b22", name = "-log10(FDR)") +
    scale_size_continuous(range = c(2.5, 8), name = "Overlap fraction") +
    labs(
      x = "Published tubular injury / failed-repair marker panel",
      y = "Rat query set",
      title = "Cross-species tubular-state marker overlap",
      subtitle = "Fisher enrichment against the detected GSE233078 gene universe"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(out_dir, "Supplementary_Figure_S6_tubular_state_marker_overlap.png"),
    p, width = 8.8, height = 4.8, dpi = 300
  )
  ggsave(
    file.path(out_dir, "Supplementary_Figure_S6_tubular_state_marker_overlap.svg"),
    p, width = 8.8, height = 4.8
  )
}

message("Done. Outputs written to: ", out_dir)
