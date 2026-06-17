#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import math
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


RSCRIPT = Path("C:/Program Files/R/R-4.5.2/bin/Rscript.exe")
CONDITIONS = ["lean", "obese", "obese_enalapril"]
FOCUS_GROUPS = ["DTL", "Mono_Macro"]


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def build_module_gene_sets() -> dict[str, list[str]]:
    return {
        "Spp1_Cd44_axis": ["Spp1", "Cd44"],
        "adhesion_migration": ["Itgam", "Itgb1", "Vcam1", "Icam1", "Fn1", "Ccl2", "Ccr2", "Cxcl10", "Mmp9"],
        "macrophage_activation": ["Lgals3", "Lyz2", "Csf1r", "Adgre1", "Trem2", "Apoe", "Tyrobp", "Cd68"],
        "fibro_inflammatory_context": ["Tgfb1", "Fn1", "Col1a1", "Col3a1", "Vim", "Timp1", "Lox"],
        "tubular_injury": ["Havcr1", "Lcn2", "Sox9", "Vcam1", "Krt8", "Krt18", "Clu"],
        "oxidative_stress": ["Hmox1", "Nfe2l2", "Gpx1", "Sod2", "Txnrd1", "Nqo1"],
        "innate_complement": ["C3", "C1qa", "C1qb", "C1qc", "Tlr4", "Nlrp3", "Il1b"],
        "chemokine_recruitment": ["Ccl2", "Ccl3", "Ccl4", "Cxcl10", "Cxcl2", "Ccr2", "Ccr5"],
        "ecm_remodeling": ["Mmp2", "Mmp9", "Timp1", "Col1a1", "Col3a1", "Fn1", "Lox"],
    }


def build_r_export_script(project_root: Path, out_dir: Path) -> str:
    gene_sets = build_module_gene_sets()
    all_genes = sorted({gene for genes in gene_sets.values() for gene in genes})
    gene_vector = ", ".join(f'"{gene}"' for gene in all_genes)
    set_lines = []
    for name, genes in gene_sets.items():
        quoted = ", ".join(f'"{gene}"' for gene in genes)
        set_lines.append(f'gene_sets[["{name}"]] <- c({quoted})')
    set_block = "\n".join(set_lines)
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
genes <- c({gene_vector})
present_genes <- intersect(genes, rownames(mat))
lib <- Matrix::colSums(mat)
expr <- as.matrix(mat[present_genes, , drop = FALSE])
norm_expr <- log1p(t(t(expr) / lib * 10000))

gene_sets <- list()
{set_block}

score_df <- data.frame(
  cell_id = focus_cells,
  sample_id = meta[focus_cells, "orig.ident"],
  condition = meta[focus_cells, "condition"],
  original_cluster = meta[focus_cells, "clusters3"],
  focus_group = meta[focus_cells, "focus_group"]
)

for (gene in genes) {{
  if (gene %in% rownames(norm_expr)) {{
    score_df[[paste0(gene, "_lognorm")]] <- as.numeric(norm_expr[gene, ])
  }} else {{
    score_df[[paste0(gene, "_lognorm")]] <- 0
  }}
}}

for (set_name in names(gene_sets)) {{
  set_genes <- intersect(gene_sets[[set_name]], rownames(norm_expr))
  if (length(set_genes) == 0) {{
    score_df[[paste0(set_name, "_score")]] <- 0
  }} else {{
    score_df[[paste0(set_name, "_score")]] <- Matrix::colMeans(norm_expr[set_genes, , drop = FALSE])
  }}
}}

write.csv(score_df, file.path(out_dir, "spp1_cd44_module_cell_scores.csv"), row.names = FALSE)

presence <- data.frame(
  module = rep(names(gene_sets), sapply(gene_sets, length)),
  requested_gene = unlist(gene_sets, use.names = FALSE)
)
presence$present <- presence$requested_gene %in% rownames(mat)
write.csv(presence, file.path(out_dir, "spp1_cd44_module_gene_presence.csv"), row.names = FALSE)
"""


def run_r_export(project_root: Path, out_dir: Path) -> None:
    if not RSCRIPT.exists():
        raise FileNotFoundError(f"Rscript not found: {RSCRIPT}")
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_export_script(project_root, out_dir))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=1800)
    finally:
        script_path.unlink(missing_ok=True)


def zscore(values: pd.Series) -> pd.Series:
    mu = float(values.mean())
    sd = float(values.std(ddof=0))
    if sd == 0 or math.isnan(sd):
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - mu) / sd


def add_anchor_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["Spp1_z"] = out.groupby("focus_group")["Spp1_lognorm"].transform(zscore)
    out["Cd44_z"] = out.groupby("focus_group")["Cd44_lognorm"].transform(zscore)
    out["anchor_score"] = np.where(out["focus_group"] == "DTL", out["Spp1_z"], out["Cd44_z"])
    return out


def spearman_corr(x: pd.Series, y: pd.Series) -> float:
    if len(x) < 3:
        return 0.0
    xr = x.rank(method="average")
    yr = y.rank(method="average")
    if float(xr.std(ddof=0)) == 0 or float(yr.std(ddof=0)) == 0:
        return 0.0
    return float(np.corrcoef(xr, yr)[0, 1])


def classify_pattern(lean: float, obese: float, enalapril: float) -> str:
    if obese > lean and enalapril < obese:
        if enalapril <= lean:
            return "obese_up_full_return"
        return "obese_up_partial_return"
    if obese > lean:
        return "obese_up_no_return"
    if obese < lean and enalapril > obese:
        return "obese_down_recovery"
    return "mixed_or_flat"


def build_sample_rows(cell_df: pd.DataFrame) -> pd.DataFrame:
    score_cols = ["anchor_score"] + [f"{name}_score" for name in build_module_gene_sets()]
    sample_df = (
        cell_df.groupby(["sample_id", "condition", "focus_group"], as_index=False)[score_cols]
        .mean()
        .sort_values(["focus_group", "condition", "sample_id"])
    )
    counts = cell_df.groupby(["sample_id", "condition", "focus_group"], as_index=False)["cell_id"].count()
    counts = counts.rename(columns={"cell_id": "n_cells"})
    return sample_df.merge(counts, on=["sample_id", "condition", "focus_group"], how="left")


def build_condition_rows(sample_df: pd.DataFrame) -> pd.DataFrame:
    module_names = list(build_module_gene_sets())
    rows = []
    for focus_group in FOCUS_GROUPS:
        for condition in CONDITIONS:
            subset = sample_df[(sample_df["focus_group"] == focus_group) & (sample_df["condition"] == condition)]
            out = {
                "focus_group": focus_group,
                "condition": condition,
                "n_samples": int(len(subset)),
                "total_cells": int(subset["n_cells"].sum()) if len(subset) else 0,
                "anchor_score_mean": float(subset["anchor_score"].mean()) if len(subset) else 0.0,
            }
            for module in module_names:
                out[f"{module}_score_mean"] = float(subset[f"{module}_score"].mean()) if len(subset) else 0.0
            rows.append(out)
    return pd.DataFrame(rows)


def build_module_support_rows(cell_df: pd.DataFrame, condition_df: pd.DataFrame) -> list[dict[str, str]]:
    module_names = list(build_module_gene_sets())
    rows: list[dict[str, str]] = []
    for focus_group in FOCUS_GROUPS:
        focus_cells = cell_df[cell_df["focus_group"] == focus_group]
        cond_values = condition_df[condition_df["focus_group"] == focus_group]
        for module in module_names:
            score_col = f"{module}_score"
            rho = spearman_corr(focus_cells["anchor_score"], focus_cells[score_col])
            value_map = {
                row["condition"]: float(row[f"{module}_score_mean"])
                for _, row in cond_values.iterrows()
            }
            lean = value_map.get("lean", 0.0)
            obese = value_map.get("obese", 0.0)
            enalapril = value_map.get("obese_enalapril", 0.0)
            rows.append(
                {
                    "focus_group": focus_group,
                    "module": module,
                    "anchor_correlation_spearman": f"{rho:.6f}",
                    "lean_score": f"{lean:.6f}",
                    "obese_score": f"{obese:.6f}",
                    "obese_enalapril_score": f"{enalapril:.6f}",
                    "obese_vs_lean_delta": f"{(obese - lean):.6f}",
                    "enalapril_vs_obese_delta": f"{(enalapril - obese):.6f}",
                    "condition_pattern": classify_pattern(lean, obese, enalapril),
                }
            )
    return rows


def save_support_heatmap(rows: list[dict[str, str]], out_path: Path) -> None:
    module_names = list(build_module_gene_sets())
    matrix = np.zeros((len(module_names), len(FOCUS_GROUPS)))
    for row in rows:
        matrix[module_names.index(row["module"]), FOCUS_GROUPS.index(row["focus_group"])] = float(row["anchor_correlation_spearman"])
    fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(FOCUS_GROUPS)))
    ax.set_xticklabels(FOCUS_GROUPS)
    ax.set_yticks(range(len(module_names)))
    ax.set_yticklabels(module_names)
    for i in range(len(module_names)):
        for j in range(len(FOCUS_GROUPS)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8, color="#111827")
    ax.set_title("Spp1/Cd44 anchor-associated module support", fontweight="bold")
    fig.colorbar(im, ax=ax, label="Spearman correlation with anchor score")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def save_pattern_bubble(rows: list[dict[str, str]], out_path: Path) -> None:
    module_names = list(build_module_gene_sets())
    fig, ax = plt.subplots(figsize=(9, 6.5), constrained_layout=True)
    color_map = {
        "obese_up_partial_return": "#2563eb",
        "obese_up_full_return": "#16a34a",
        "obese_up_no_return": "#dc2626",
        "obese_down_recovery": "#7c3aed",
        "mixed_or_flat": "#94a3b8",
    }
    for row in rows:
        x = FOCUS_GROUPS.index(row["focus_group"])
        y = module_names.index(row["module"])
        rho = abs(float(row["anchor_correlation_spearman"]))
        ax.scatter(
            x,
            y,
            s=1200 * rho + 35,
            color=color_map.get(row["condition_pattern"], "#94a3b8"),
            alpha=0.82,
            edgecolor="white",
            linewidth=0.8,
        )
    ax.set_xticks(range(len(FOCUS_GROUPS)))
    ax.set_xticklabels(FOCUS_GROUPS)
    ax.set_yticks(range(len(module_names)))
    ax.set_yticklabels(module_names)
    ax.set_title("Condition pattern of Spp1/Cd44-associated modules", fontweight="bold")
    ax.grid(color="#e5e7eb", linewidth=0.8)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, label=label, markersize=9)
        for label, color in color_map.items()
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.35, 1.0))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def build_summary(rows: list[dict[str, str]], presence: pd.DataFrame) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def build_readable_result(summary: str) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_to_quick(project_root: Path, out_dir: Path) -> None:
    return None

def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results/public_data_gse233078_spp1_cd44_module_support"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_r_export(project_root, out_dir)
    cell_df = pd.read_csv(out_dir / "spp1_cd44_module_cell_scores.csv")
    cell_df = add_anchor_scores(cell_df)
    cell_df.to_csv(out_dir / "spp1_cd44_module_cell_scores_with_anchor.csv", index=False, encoding="utf-8-sig")
    sample_df = build_sample_rows(cell_df)
    sample_df.to_csv(out_dir / "spp1_cd44_module_sample_scores.csv", index=False, encoding="utf-8-sig")
    condition_df = build_condition_rows(sample_df)
    condition_df.to_csv(out_dir / "spp1_cd44_module_condition_scores.csv", index=False, encoding="utf-8-sig")

    support_rows = build_module_support_rows(cell_df, condition_df)
    write_csv(
        out_dir / "spp1_cd44_module_support_rows.csv",
        support_rows,
        [
            "focus_group",
            "module",
            "anchor_correlation_spearman",
            "lean_score",
            "obese_score",
            "obese_enalapril_score",
            "obese_vs_lean_delta",
            "enalapril_vs_obese_delta",
            "condition_pattern",
        ],
    )
    save_support_heatmap(support_rows, out_dir / "FigureAF_gse233078_spp1_cd44_module_support_heatmap.svg")
    save_pattern_bubble(support_rows, out_dir / "FigureAG_gse233078_spp1_cd44_module_support_bubble.svg")
    presence = pd.read_csv(out_dir / "spp1_cd44_module_gene_presence.csv")
    summary = build_summary(support_rows, presence)
    write_text(out_dir / "spp1_cd44_module_support_summary.txt", summary)
    write_text(out_dir / "status.txt", "Spp1/Cd44 module support completed.\n")
    print(f"Spp1/Cd44 module support outputs written to: {out_dir}")


if __name__ == "__main__":
    main()


