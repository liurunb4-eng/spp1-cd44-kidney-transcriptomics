#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "public_data_gse233078_dtl_annotation_audit"
RSCRIPT = Path(r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe")


R_CODE = r'''
suppressPackageStartupMessages(library(Matrix))

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "results/public_data_gse233078_dtl_annotation_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_meta.data.csv.gz")
mat_path <- file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_counts_postfilter.rds.gz")

meta <- read.csv(meta_path, check.names = FALSE)
colnames(meta)[colnames(meta) == "Unnamed: 0"] <- "cell_id"
colnames(meta)[colnames(meta) == ""] <- "cell_id"
meta$condition <- ifelse(meta$exp.cond == "Lean", "lean", ifelse(meta$exp.cond == "Obese", "obese", "obese_enalapril"))
rownames(meta) <- meta$cell_id

con <- gzcon(file(mat_path, "rb"))
raw1 <- readBin(con, what = "raw", n = 300000000)
close(con)
raw2 <- memDecompress(raw1, type = "gzip")
mat <- unserialize(raw2)

cells <- intersect(meta$cell_id, colnames(mat))
meta <- meta[cells, , drop = FALSE]
mat <- mat[, cells, drop = FALSE]
lib <- Matrix::colSums(mat)

marker_sets <- list(
  dtl_thin_limb_candidate = c("Aqp1", "Slc14a2", "Cldn10", "Cldn16", "Pdpn"),
  proximal_tubule = c("Lrp2", "Cubn", "Slc5a2", "Slc34a1", "Slc13a3", "Acy3", "Kap"),
  injured_pt = c("Havcr1", "Vcam1", "Krt8", "Krt18", "Lcn2", "Sox9", "Clu"),
  tal = c("Slc12a1", "Umod", "Pvalb", "Kcnj1"),
  dct_cnt = c("Slc12a3", "Calb1", "Trpm6", "Pvalb"),
  collecting_duct = c("Aqp2", "Aqp3", "Avpr2", "Scnn1g", "Fxyd4"),
  mono_macro = c("Lyz2", "Csf1r", "Adgre1", "Cd68", "C1qa", "C1qb", "C1qc", "Tyrobp", "Apoe", "Trem2"),
  spp1_cd44_anchor = c("Spp1", "Cd44")
)
genes <- unique(unlist(marker_sets))
present_genes <- intersect(genes, rownames(mat))
expr <- as.matrix(mat[present_genes, , drop = FALSE])
norm_expr <- log1p(t(t(expr) / lib * 10000))

presence <- data.frame(
  marker_set = rep(names(marker_sets), lengths(marker_sets)),
  gene = unlist(marker_sets, use.names = FALSE)
)
presence$present <- presence$gene %in% rownames(mat)
write.csv(presence, file.path(out_dir, "dtl_annotation_marker_gene_presence.csv"), row.names = FALSE)

cluster_levels <- sort(unique(meta$clusters3))
rows <- list()
idx <- 1
for (cluster in cluster_levels) {
  cluster_cells <- rownames(meta)[meta$clusters3 == cluster]
  n_cells <- length(cluster_cells)
  for (gene in genes) {
    present <- gene %in% rownames(norm_expr)
    if (present) {
      vals <- as.numeric(norm_expr[gene, cluster_cells])
      raw_counts <- mat[gene, cluster_cells]
      mean_expr <- mean(vals)
      pct_detected <- mean(raw_counts > 0)
    } else {
      mean_expr <- 0
      pct_detected <- 0
    }
    rows[[idx]] <- data.frame(
      cluster = cluster,
      n_cells = n_cells,
      gene = gene,
      present = present,
      mean_lognorm = mean_expr,
      pct_detected = pct_detected
    )
    idx <- idx + 1
  }
}
gene_cluster <- do.call(rbind, rows)
write.csv(gene_cluster, file.path(out_dir, "dtl_annotation_marker_by_cluster.csv"), row.names = FALSE)

score_df <- data.frame(
  cell_id = rownames(meta),
  sample_id = meta$orig.ident,
  condition = meta$condition,
  cluster = meta$clusters3
)
for (set_name in names(marker_sets)) {
  set_genes <- intersect(marker_sets[[set_name]], rownames(norm_expr))
  score_col <- paste0(set_name, "_score")
  if (length(set_genes) == 0) {
    score_df[[score_col]] <- 0
  } else {
    score_df[[score_col]] <- Matrix::colMeans(norm_expr[set_genes, , drop = FALSE])
  }
}
write.csv(score_df, file.path(out_dir, "dtl_annotation_cell_marker_scores.csv"), row.names = FALSE)

cluster_score_rows <- list()
idx <- 1
for (cluster in cluster_levels) {
  sub <- score_df[score_df$cluster == cluster, ]
  out <- data.frame(cluster = cluster, n_cells = nrow(sub))
  for (set_name in names(marker_sets)) {
    score_col <- paste0(set_name, "_score")
    out[[paste0(score_col, "_mean")]] <- mean(sub[[score_col]])
    out[[paste0(score_col, "_median")]] <- median(sub[[score_col]])
  }
  cluster_score_rows[[idx]] <- out
  idx <- idx + 1
}
cluster_scores <- do.call(rbind, cluster_score_rows)
write.csv(cluster_scores, file.path(out_dir, "dtl_annotation_marker_scores_by_cluster.csv"), row.names = FALSE)

condition_rows <- list()
idx <- 1
focus <- score_df[score_df$cluster %in% c("DTL", "PCT/PST", "Mono", "Macro", "Mono/Macro"), ]
for (cluster in sort(unique(focus$cluster))) {
  for (condition in c("lean", "obese", "obese_enalapril")) {
    sub <- focus[focus$cluster == cluster & focus$condition == condition, ]
    out <- data.frame(cluster = cluster, condition = condition, n_cells = nrow(sub))
    for (set_name in names(marker_sets)) {
      score_col <- paste0(set_name, "_score")
      out[[paste0(score_col, "_mean")]] <- ifelse(nrow(sub) > 0, mean(sub[[score_col]]), NA)
    }
    condition_rows[[idx]] <- out
    idx <- idx + 1
  }
}
condition_scores <- do.call(rbind, condition_rows)
write.csv(condition_scores, file.path(out_dir, "dtl_annotation_marker_scores_by_condition.csv"), row.names = FALSE)

selected_genes <- c("Spp1", "Cd44", "Aqp1", "Slc14a2", "Cldn10", "Lrp2", "Slc5a2", "Slc34a1", "Havcr1", "Vcam1", "Krt8", "Lcn2", "Csf1r", "Adgre1", "C1qa")
selected <- gene_cluster[gene_cluster$gene %in% selected_genes & gene_cluster$cluster %in% c("DTL", "PCT/PST", "Mono", "Macro", "Mono/Macro", "TAL", "DCT"), ]
write.csv(selected, file.path(out_dir, "dtl_annotation_selected_marker_summary.csv"), row.names = FALSE)

dtl_scores <- cluster_scores[cluster_scores$cluster == "DTL", ]
pct_scores <- cluster_scores[cluster_scores$cluster == "PCT/PST", ]
mono_scores <- cluster_scores[cluster_scores$cluster %in% c("Mono", "Macro", "Mono/Macro"), ]

safe_num <- function(x) { if (length(x) == 0 || is.na(x)) return(0); as.numeric(x) }
dtl_n <- safe_num(dtl_scores$n_cells[1])
pct_n <- safe_num(pct_scores$n_cells[1])
dtl_dtl <- safe_num(dtl_scores$dtl_thin_limb_candidate_score_mean[1])
dtl_pt <- safe_num(dtl_scores$proximal_tubule_score_mean[1])
dtl_inj <- safe_num(dtl_scores$injured_pt_score_mean[1])
pct_pt <- safe_num(pct_scores$proximal_tubule_score_mean[1])
pct_inj <- safe_num(pct_scores$injured_pt_score_mean[1])
dtl_spp1 <- safe_num(gene_cluster$mean_lognorm[gene_cluster$cluster == "DTL" & gene_cluster$gene == "Spp1"][1])
dtl_spp1_pct <- safe_num(gene_cluster$pct_detected[gene_cluster$cluster == "DTL" & gene_cluster$gene == "Spp1"][1])
dtl_cd44 <- safe_num(gene_cluster$mean_lognorm[gene_cluster$cluster == "DTL" & gene_cluster$gene == "Cd44"][1])
mono_cd44_mean <- mean(gene_cluster$mean_lognorm[gene_cluster$cluster %in% c("Mono", "Macro", "Mono/Macro") & gene_cluster$gene == "Cd44"], na.rm = TRUE)

summary <- c(
  "GSE233078 DTL annotation audit",
  paste0("Cells analyzed: ", nrow(meta), "; DTL cells: ", dtl_n, "; PCT/PST cells: ", pct_n, "."),
  paste0("DTL marker-set mean score in DTL: ", sprintf("%.3f", dtl_dtl), "."),
  paste0("Proximal-tubule marker-set mean score in DTL: ", sprintf("%.3f", dtl_pt), "; in PCT/PST: ", sprintf("%.3f", pct_pt), "."),
  paste0("Injured-PT marker-set mean score in DTL: ", sprintf("%.3f", dtl_inj), "; in PCT/PST: ", sprintf("%.3f", pct_inj), "."),
  paste0("DTL Spp1 mean lognorm: ", sprintf("%.3f", dtl_spp1), "; DTL Spp1 detection fraction: ", sprintf("%.3f", dtl_spp1_pct), "."),
  paste0("DTL Cd44 mean lognorm: ", sprintf("%.3f", dtl_cd44), "; Mono/Macro-family Cd44 mean lognorm: ", sprintf("%.3f", mono_cd44_mean), "."),
  "Interpretation rule: support for DTL identity is stronger when DTL candidate markers exceed proximal-tubule and injured-PT marker scores within the DTL cluster, while proximal-tubule markers are highest in PCT/PST. If DTL shows high injured-PT or proximal-tubule scores, the Spp1 source should be framed as DTL-associated/injured tubular rather than definitive DTL."
)
writeLines(summary, file.path(out_dir, "dtl_annotation_audit_summary.txt"))
writeLines("OK", file.path(out_dir, "status.txt"))
'''


def main() -> None:
    if not RSCRIPT.exists():
        raise FileNotFoundError(f"Rscript not found: {RSCRIPT}")
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(R_CODE)
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], cwd=ROOT, check=True, timeout=1800)
    finally:
        script_path.unlink(missing_ok=True)
    print(f"Wrote DTL annotation audit to {OUT}")


if __name__ == "__main__":
    main()

