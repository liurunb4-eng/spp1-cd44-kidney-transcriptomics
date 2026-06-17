#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import math
import shutil
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gseapy as gp


RSCRIPT = Path("C:/Program Files/R/R-4.5.2/bin/Rscript.exe")
FOCUS_GROUPS = ["DTL", "Mono_Macro"]
CONTRASTS = [
    ("obese_vs_lean", "obese", "lean", "up"),
    ("enalapril_vs_obese", "obese_enalapril", "obese", "down"),
]
ENRICHR_LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2019_Mouse", "Reactome_2022"]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_r_pseudobulk_script(project_root: Path, out_dir: Path) -> str:
    return f"""
suppressPackageStartupMessages(library(Matrix))
project_root <- "{project_root.as_posix()}"
out_dir <- "{out_dir.as_posix()}"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_meta.data.csv.gz"), check.names = FALSE)
colnames(meta)[colnames(meta) == "Unnamed: 0"] <- "cell_id"
colnames(meta)[colnames(meta) == ""] <- "cell_id"
meta$condition <- ifelse(meta$exp.cond == "Lean", "lean", ifelse(meta$exp.cond == "Obese", "obese", "obese_enalapril"))
meta$focus_group <- ifelse(meta$clusters3 == "DTL", "DTL", ifelse(meta$clusters3 %in% c("Mono", "Macro", "Mono/Macro"), "Mono_Macro", NA))
meta <- meta[!is.na(meta$focus_group), ]
rownames(meta) <- meta$cell_id

con <- gzcon(file(file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_counts_postfilter.rds.gz"), "rb"))
raw1 <- readBin(con, what="raw", n=300000000)
close(con)
raw2 <- memDecompress(raw1, type="gzip")
mat <- unserialize(raw2)

focus_cells <- intersect(meta$cell_id, colnames(mat))
meta <- meta[focus_cells, , drop = FALSE]
mat <- mat[, focus_cells, drop = FALSE]
meta$group_id <- paste(meta$orig.ident, meta$condition, meta$focus_group, sep = "__")
group_ids <- unique(meta$group_id)

pseudobulk <- matrix(0, nrow = nrow(mat), ncol = length(group_ids))
rownames(pseudobulk) <- rownames(mat)
colnames(pseudobulk) <- group_ids
group_meta <- data.frame()

for (group_id in group_ids) {{
  cells <- rownames(meta)[meta$group_id == group_id]
  pseudobulk[, group_id] <- Matrix::rowSums(mat[, cells, drop = FALSE])
  parts <- strsplit(group_id, "__", fixed = TRUE)[[1]]
  group_meta <- rbind(
    group_meta,
    data.frame(
      group_id = group_id,
      sample_id = parts[1],
      condition = parts[2],
      focus_group = parts[3],
      n_cells = length(cells),
      library_size = sum(pseudobulk[, group_id])
    )
  )
}}

write.csv(data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE), file.path(out_dir, "pseudobulk_counts_focus_groups.csv"), row.names = FALSE)
write.csv(group_meta, file.path(out_dir, "pseudobulk_group_metadata.csv"), row.names = FALSE)
"""


def run_r_export(project_root: Path, out_dir: Path) -> None:
    if not RSCRIPT.exists():
        raise FileNotFoundError(f"Rscript not found: {RSCRIPT}")
    counts_path = out_dir / "pseudobulk_counts_focus_groups.csv"
    meta_path = out_dir / "pseudobulk_group_metadata.csv"
    if counts_path.exists() and meta_path.exists():
        return
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_pseudobulk_script(project_root, out_dir))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=1800)
    finally:
        script_path.unlink(missing_ok=True)


def cohen_d(group_a: np.ndarray, group_b: np.ndarray) -> np.ndarray:
    mean_a = np.nanmean(group_a, axis=1)
    mean_b = np.nanmean(group_b, axis=1)
    var_a = np.nanvar(group_a, axis=1, ddof=1)
    var_b = np.nanvar(group_b, axis=1, ddof=1)
    n_a = group_a.shape[1]
    n_b = group_b.shape[1]
    pooled = (((n_a - 1) * var_a) + ((n_b - 1) * var_b)) / max(n_a + n_b - 2, 1)
    pooled[pooled <= 0] = np.nan
    out = (mean_a - mean_b) / np.sqrt(pooled)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def build_logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    numeric = counts.set_index("gene")
    lib_sizes = numeric.sum(axis=0)
    cpm = numeric.divide(lib_sizes, axis=1) * 1_000_000
    return np.log2(cpm + 0.5)


def run_deg_for_contrast(
    logcpm: pd.DataFrame,
    metadata: pd.DataFrame,
    focus_group: str,
    contrast_name: str,
    condition_a: str,
    condition_b: str,
) -> pd.DataFrame:
    group_a_ids = metadata[(metadata["focus_group"] == focus_group) & (metadata["condition"] == condition_a)]["group_id"].tolist()
    group_b_ids = metadata[(metadata["focus_group"] == focus_group) & (metadata["condition"] == condition_b)]["group_id"].tolist()
    a = logcpm[group_a_ids].to_numpy()
    b = logcpm[group_b_ids].to_numpy()
    mean_a = np.nanmean(a, axis=1)
    mean_b = np.nanmean(b, axis=1)
    logfc = mean_a - mean_b
    stat = ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    p_values = np.nan_to_num(stat.pvalue, nan=1.0, posinf=1.0, neginf=1.0)
    fdr = multipletests(p_values, method="fdr_bh")[1]
    effect = cohen_d(a, b)
    out = pd.DataFrame(
        {
            "gene": logcpm.index,
            "focus_group": focus_group,
            "contrast": contrast_name,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "n_a": len(group_a_ids),
            "n_b": len(group_b_ids),
            "mean_logcpm_a": mean_a,
            "mean_logcpm_b": mean_b,
            "logFC": logfc,
            "cohen_d": effect,
            "welch_p": p_values,
            "fdr": fdr,
        }
    )
    out["signed_score"] = np.sign(out["logFC"]) * -np.log10(out["welch_p"].clip(lower=1e-300))
    out["is_anchor_gene"] = out["gene"].isin(["Spp1", "Cd44"])
    return out.sort_values(["welch_p", "fdr"], ascending=True)


def select_candidate_genes(deg: pd.DataFrame, desired_direction: str, min_genes: int = 20, fallback_n: int = 250) -> tuple[list[str], str]:
    if desired_direction == "up":
        directed = deg[deg["logFC"] > 0].copy()
        directed["rank_score"] = directed["signed_score"]
    else:
        directed = deg[deg["logFC"] < 0].copy()
        directed["rank_score"] = -directed["signed_score"]
    nominal = directed[(directed["welch_p"] < 0.1) & (directed["rank_score"] > 0)].sort_values("rank_score", ascending=False)
    if len(nominal) >= min_genes:
        return nominal["gene"].head(500).tolist(), "nominal_p_lt_0.1"
    fallback = directed.sort_values("rank_score", ascending=False).head(fallback_n)
    return fallback["gene"].tolist(), f"fallback_top_{fallback_n}_ranked"


def run_enrichr(gene_lists: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    frames = []
    for _, row in gene_lists.iterrows():
        genes = [gene for gene in str(row["genes"]).split(";") if gene]
        if len(genes) < 5:
            continue
        for library in ENRICHR_LIBRARIES:
            try:
                # Query one library at a time so rate-limit failures cannot mix library labels.
                time.sleep(1.2)
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=library,
                    organism="Mouse",
                    outdir=None,
                    cutoff=1.0,
                )
                res = enr.res2d.copy()
                if res.empty:
                    continue
                res["requested_library"] = library
                res["list_name"] = row["list_name"]
                res["focus_group"] = row["focus_group"]
                res["contrast"] = row["contrast"]
                res["direction"] = row["direction"]
                res["list_basis"] = row["list_basis"]
                frames.append(res)
            except Exception as exc:
                safe_name = f"{row['list_name']}_{library}".replace("/", "_")
                write_text(out_dir / f"enrichr_error_{safe_name}.txt", f"{type(exc).__name__}: {exc}\n")
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_gene_lists(deg_all: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for focus_group in FOCUS_GROUPS:
        for contrast_name, condition_a, condition_b, direction in CONTRASTS:
            subset = deg_all[(deg_all["focus_group"] == focus_group) & (deg_all["contrast"] == contrast_name)]
            genes, basis = select_candidate_genes(subset, direction)
            rows.append(
                {
                    "list_name": f"{focus_group}_{contrast_name}_{direction}",
                    "focus_group": focus_group,
                    "contrast": contrast_name,
                    "direction": direction,
                    "list_basis": basis,
                    "n_genes": len(genes),
                    "genes": ";".join(genes),
                }
            )
    return pd.DataFrame(rows)


def save_volcano(deg_all: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    panel_order = [
        ("DTL", "obese_vs_lean"),
        ("DTL", "enalapril_vs_obese"),
        ("Mono_Macro", "obese_vs_lean"),
        ("Mono_Macro", "enalapril_vs_obese"),
    ]
    for ax, (focus_group, contrast) in zip(axes.flat, panel_order):
        subset = deg_all[(deg_all["focus_group"] == focus_group) & (deg_all["contrast"] == contrast)].copy()
        y = -np.log10(subset["welch_p"].clip(lower=1e-300))
        sig = (subset["welch_p"] < 0.05) & (subset["logFC"].abs() > 0.5)
        ax.scatter(subset["logFC"], y, s=3, color="#94a3b8", alpha=0.45, linewidths=0)
        ax.scatter(subset.loc[sig, "logFC"], y.loc[sig], s=6, color="#dc2626", alpha=0.65, linewidths=0)
        for gene in ["Spp1", "Cd44"]:
            row = subset[subset["gene"] == gene]
            if not row.empty:
                x_val = float(row["logFC"].iloc[0])
                y_val = float(-np.log10(max(float(row["welch_p"].iloc[0]), 1e-300)))
                ax.scatter([x_val], [y_val], s=36, color="#2563eb", edgecolor="white", linewidth=0.8)
                ax.text(x_val, y_val, gene, fontsize=8, ha="left", va="bottom")
        ax.axvline(0, color="#111827", linewidth=0.8)
        ax.set_title(f"{focus_group} | {contrast}", fontweight="bold")
        ax.set_xlabel("logFC")
        ax.set_ylabel("-log10(p)")
        ax.grid(color="#e5e7eb", linewidth=0.5)
    fig.suptitle("Spp1/Cd44-axis pseudobulk DEG volcano plots", fontsize=14, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def top_enrichment(enrichment: pd.DataFrame, n_per_list: int = 6) -> pd.DataFrame:
    if enrichment.empty:
        return enrichment
    p_col = "Adjusted P-value" if "Adjusted P-value" in enrichment.columns else "P-value"
    return (
        enrichment.sort_values(p_col, ascending=True)
        .groupby(["list_name", "requested_library"], as_index=False)
        .head(n_per_list)
        .reset_index(drop=True)
    )


def save_enrichment_dotplot(enrichment: pd.DataFrame, out_path: Path) -> None:
    top = top_enrichment(enrichment, n_per_list=4)
    if top.empty:
        write_text(out_path.with_suffix(".txt"), "No enrichment results available.\n")
        return
    p_col = "Adjusted P-value" if "Adjusted P-value" in top.columns else "P-value"
    top = top.copy()
    top["minus_log10_fdr"] = -np.log10(top[p_col].astype(float).clip(lower=1e-300))
    top["display_term"] = top["Term"].astype(str).str.replace(r" \\(GO:\\d+\\)$", "", regex=True).str.slice(0, 56)
    top = top.sort_values(["list_name", "minus_log10_fdr"], ascending=[True, False]).head(40)
    y_labels = [f"{row['list_name']} | {row['display_term']}" for _, row in top.iterrows()]
    y = np.arange(len(top))
    color_map = {"GO_Biological_Process_2023": "#2563eb", "KEGG_2019_Mouse": "#dc2626", "Reactome_2022": "#16a34a"}
    colors = [color_map.get(gs, "#64748b") for gs in top["requested_library"]]
    fig, ax = plt.subplots(figsize=(11, max(6, len(top) * 0.26)), constrained_layout=True)
    ax.scatter(top["minus_log10_fdr"], y, s=np.sqrt(top["Combined Score"].astype(float).clip(lower=0)) * 10 + 25, c=colors, alpha=0.82)
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_title("GO/KEGG/Reactome enrichment of Spp1/Cd44-axis DE candidates", fontweight="bold")
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, label=label, markersize=8)
        for label, color in color_map.items()
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower right")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def build_summary(deg_all: pd.DataFrame, gene_lists: pd.DataFrame, enrichment: pd.DataFrame) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def build_readable_result(summary: str) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_to_quick(project_root: Path, out_dir: Path) -> None:
    return None

def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results/public_data_gse233078_spp1_cd44_deg_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_r_export(project_root, out_dir)
    counts = pd.read_csv(out_dir / "pseudobulk_counts_focus_groups.csv")
    metadata = pd.read_csv(out_dir / "pseudobulk_group_metadata.csv")
    logcpm = build_logcpm(counts)
    deg_frames = []
    for focus_group in FOCUS_GROUPS:
        for contrast_name, condition_a, condition_b, _direction in CONTRASTS:
            deg_frames.append(run_deg_for_contrast(logcpm, metadata, focus_group, contrast_name, condition_a, condition_b))
    deg_all = pd.concat(deg_frames, ignore_index=True)
    deg_all.to_csv(out_dir / "spp1_cd44_axis_deg_all.csv", index=False, encoding="utf-8-sig")
    priority = deg_all[(deg_all["is_anchor_gene"]) | ((deg_all["welch_p"] < 0.05) & (deg_all["logFC"].abs() >= 0.5))]
    priority.to_csv(out_dir / "spp1_cd44_axis_deg_priority.csv", index=False, encoding="utf-8-sig")
    gene_lists = build_gene_lists(deg_all)
    gene_lists.to_csv(out_dir / "spp1_cd44_axis_enrichment_gene_lists.csv", index=False, encoding="utf-8-sig")
    enrichment = run_enrichr(gene_lists, out_dir)
    enrichment.to_csv(out_dir / "spp1_cd44_axis_enrichr_results.csv", index=False, encoding="utf-8-sig")
    top_enrichment(enrichment).to_csv(out_dir / "spp1_cd44_axis_enrichr_top_terms.csv", index=False, encoding="utf-8-sig")
    save_volcano(deg_all, out_dir / "FigureA1_gse233078_spp1_cd44_deg_volcano.svg")
    save_enrichment_dotplot(enrichment, out_dir / "FigureA2_gse233078_spp1_cd44_enrichment_dotplot.svg")
    summary = build_summary(deg_all, gene_lists, enrichment)
    write_text(out_dir / "spp1_cd44_deg_enrichment_summary.txt", summary)
    write_text(out_dir / "status.txt", "Spp1/Cd44 DEG and enrichment completed.\n")
    print(f"Spp1/Cd44 DEG/enrichment outputs written to: {out_dir}")


if __name__ == "__main__":
    main()


