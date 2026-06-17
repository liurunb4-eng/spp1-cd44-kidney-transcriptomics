#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import html
import shutil
import subprocess
import tempfile


RSCRIPT = Path(r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def condition_order() -> list[str]:
    return ["lean", "obese", "obese_enalapril"]


def build_gene_sets() -> dict[str, list[str]]:
    return {
        "spp1_cd44_axis": ["Spp1", "Cd44"],
        "tubular_stress": ["Spp1", "Havcr1", "Lcn2", "Vcam1", "Sox9"],
        "myeloid_activation": ["Spp1", "Cd44", "Tlr4", "Nlrp3", "Lgals3", "Ccl2"],
        "fibro_inflammatory_context": ["Tgfb1", "Fn1", "Col1a1", "Col3a1", "Tlr4", "Nlrp3"],
    }


def focused_cluster_map() -> dict[str, str]:
    return {
        "DTL": "DTL",
        "Mono": "Mono_Macro",
        "Macro": "Mono_Macro",
        "Mono/Macro": "Mono_Macro",
    }


def build_r_script(project_root: Path, out_dir: Path) -> str:
    gene_sets = build_gene_sets()
    all_genes = sorted({gene for genes in gene_sets.values() for gene in genes})
    gene_vector = ", ".join(f'"{gene}"' for gene in all_genes)
    set_lines = []
    for name, genes in gene_sets.items():
        quoted_genes = ", ".join(f'"{gene}"' for gene in genes)
        set_lines.append(f'gene_sets[["{name}"]] <- c({quoted_genes})')
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
expr <- as.matrix(mat[present_genes, , drop = FALSE])
lib <- Matrix::colSums(mat[, focus_cells, drop = FALSE])
norm_expr <- log1p(t(t(expr) / lib * 10000))

gene_sets <- list()
{set_block}
score_df <- data.frame(
  cell_id = focus_cells,
  sample_id = meta[focus_cells, "orig.ident"],
  condition = meta[focus_cells, "condition"],
  original_cluster = meta[focus_cells, "clusters3"],
  focus_group = meta[focus_cells, "focus_group"],
  UMAP_1 = meta[focus_cells, "UMAP_1"],
  UMAP_2 = meta[focus_cells, "UMAP_2"],
  nCount_RNA = meta[focus_cells, "nCount_RNA"],
  nFeature_RNA = meta[focus_cells, "nFeature_RNA"]
)
for (set_name in names(gene_sets)) {{
  set_genes <- intersect(gene_sets[[set_name]], rownames(norm_expr))
  if (length(set_genes) == 0) {{
    score_df[[set_name]] <- 0
  }} else {{
    score_df[[set_name]] <- Matrix::colMeans(norm_expr[set_genes, , drop = FALSE])
  }}
}}
write.csv(score_df, file.path(out_dir, "dtl_monomacro_cell_state_scores.csv"), row.names = FALSE)

present_rows <- data.frame(
  gene_set = rep(names(gene_sets), sapply(gene_sets, length)),
  requested_gene = unlist(gene_sets, use.names = FALSE)
)
present_rows$present <- present_rows$requested_gene %in% present_genes
write.csv(present_rows, file.path(out_dir, "dtl_monomacro_state_gene_presence.csv"), row.names = FALSE)
"""


def run_r_export(project_root: Path, out_dir: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_script(project_root, out_dir))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=1800)
    finally:
        script_path.unlink(missing_ok=True)


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def build_sample_state_rows(cell_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in cell_rows:
        grouped.setdefault((row["sample_id"], row["condition"], row["focus_group"]), []).append(row)
    out_rows: list[dict[str, str]] = []
    score_cols = list(build_gene_sets().keys())
    for (sample_id, condition, focus_group), rows in sorted(grouped.items()):
        out = {
            "sample_id": sample_id,
            "condition": condition,
            "focus_group": focus_group,
            "n_cells": str(len(rows)),
        }
        for col in score_cols:
            out[col] = f"{mean([float(row[col]) for row in rows]):.6f}"
        out_rows.append(out)
    return out_rows


def zscore_by_group(rows: list[dict[str, str]], value_col: str) -> dict[int, float]:
    values_by_group: dict[str, list[float]] = {}
    for row in rows:
        values_by_group.setdefault(row["focus_group"], []).append(float(row[value_col]))
    params: dict[str, tuple[float, float]] = {}
    for group, values in values_by_group.items():
        mu = mean(values)
        variance = mean([(value - mu) ** 2 for value in values])
        sd = variance ** 0.5 if variance > 0 else 1.0
        params[group] = (mu, sd)
    out: dict[int, float] = {}
    for idx, row in enumerate(rows):
        mu, sd = params[row["focus_group"]]
        out[idx] = (float(row[value_col]) - mu) / sd
    return out


def build_condition_state_rows(sample_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    state_cols = ["spp1_cd44_axis", "tubular_stress", "myeloid_activation", "fibro_inflammatory_context"]
    z_by_col = {col: zscore_by_group(sample_rows, col) for col in state_cols}
    enriched_rows = []
    for idx, row in enumerate(sample_rows):
        enriched = dict(row)
        if row["focus_group"] == "DTL":
            state_axis = mean([z_by_col["spp1_cd44_axis"][idx], z_by_col["tubular_stress"][idx], z_by_col["fibro_inflammatory_context"][idx]])
        else:
            state_axis = mean([z_by_col["spp1_cd44_axis"][idx], z_by_col["myeloid_activation"][idx], z_by_col["fibro_inflammatory_context"][idx]])
        enriched["state_axis_z"] = f"{state_axis:.6f}"
        enriched_rows.append(enriched)

    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in enriched_rows:
        grouped.setdefault((row["focus_group"], row["condition"]), []).append(row)
    out_rows: list[dict[str, str]] = []
    for focus_group in ["DTL", "Mono_Macro"]:
        for condition in condition_order():
            rows = grouped.get((focus_group, condition), [])
            out = {
                "focus_group": focus_group,
                "condition": condition,
                "n_samples": str(len(rows)),
                "total_cells": str(sum(int(row["n_cells"]) for row in rows)),
                "state_axis_z_mean": f"{mean([float(row['state_axis_z']) for row in rows]):.6f}",
            }
            for col in state_cols:
                out[f"{col}_mean"] = f"{mean([float(row[col]) for row in rows]):.6f}"
            out_rows.append(out)
    return out_rows


def classify_state_shift(lean: float, obese: float, enalapril: float) -> str:
    if obese > lean and enalapril < obese:
        if enalapril <= lean:
            return "obese_up_with_full_reversal"
        return "obese_up_with_partial_reversal"
    if obese > lean:
        return "obese_up_without_reversal"
    if obese < lean and enalapril > obese:
        return "obese_down_with_recovery"
    return "mixed_or_flat"


def build_state_shift_rows(condition_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    value_map = {
        (row["focus_group"], row["condition"]): float(row["state_axis_z_mean"])
        for row in condition_rows
    }
    out_rows = []
    for focus_group in ["DTL", "Mono_Macro"]:
        lean = value_map.get((focus_group, "lean"), 0.0)
        obese = value_map.get((focus_group, "obese"), 0.0)
        enalapril = value_map.get((focus_group, "obese_enalapril"), 0.0)
        out_rows.append(
            {
                "focus_group": focus_group,
                "lean_state_axis_z": f"{lean:.6f}",
                "obese_state_axis_z": f"{obese:.6f}",
                "obese_enalapril_state_axis_z": f"{enalapril:.6f}",
                "obese_vs_lean_delta": f"{(obese - lean):.6f}",
                "enalapril_vs_obese_delta": f"{(enalapril - obese):.6f}",
                "state_shift_call": classify_state_shift(lean, obese, enalapril),
            }
        )
    return out_rows


def esc(text: str) -> str:
    return html.escape(str(text), quote=True)


def svg_text(svg: list[str], x: float, y: float, text: str, size: int = 18, weight: str = "400", fill: str = "#0f172a", anchor: str = "start") -> None:
    svg.append(
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}">{esc(text)}</text>'
    )


def scale(value: float, in_min: float, in_max: float, out_min: float, out_max: float) -> float:
    if in_max <= in_min:
        return (out_min + out_max) / 2
    return out_min + (value - in_min) / (in_max - in_min) * (out_max - out_min)


def build_umap_svg(cell_rows: list[dict[str, str]], out_path: Path) -> None:
    width, height = 1700, 980
    xs = [float(row["UMAP_1"]) for row in cell_rows]
    ys = [float(row["UMAP_2"]) for row in cell_rows]
    colors = {"lean": "#94a3b8", "obese": "#b91c1c", "obese_enalapril": "#0f766e"}
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    svg_text(svg, 60, 78, "Figure AB. DTL and Mono/Macro state-shift map", 38, "700")
    svg_text(svg, 60, 118, "UMAP positions from GSE233078 metadata; highlighted cells are DTL and Mono/Macro only.", 20, "400", "#475569")
    panel = {"DTL": (70, 175, 730, 660), "Mono_Macro": (900, 175, 730, 660)}
    for group, (x0, y0, w, h) in panel.items():
        svg.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" rx="24" fill="white" stroke="#cbd5e1" stroke-width="2"/>')
        svg_text(svg, x0 + 28, y0 + 44, group, 26, "700")
        for row in cell_rows:
            if row["focus_group"] != group:
                continue
            x = scale(float(row["UMAP_1"]), min(xs), max(xs), x0 + 45, x0 + w - 45)
            y = scale(float(row["UMAP_2"]), min(ys), max(ys), y0 + h - 55, y0 + 75)
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.1" fill="{colors[row["condition"]]}" fill-opacity="0.46"/>')
    legend_x = 80
    for idx, condition in enumerate(condition_order()):
        y = 900 + idx * 28
        svg.append(f'<circle cx="{legend_x}" cy="{y}" r="8" fill="{colors[condition]}" fill-opacity="0.85"/>')
        svg_text(svg, legend_x + 20, y + 6, condition, 18, "600", "#334155")
    svg_text(svg, 900, 904, "Reading rule: this is a condition/state distribution map, not causal pseudotime.", 18, "700", "#475569")
    svg.append("</svg>")
    out_path.write_text("\n".join(svg), encoding="utf-8")


def build_state_shift_svg(condition_rows: list[dict[str, str]], shift_rows: list[dict[str, str]], out_path: Path) -> None:
    width, height = 1700, 980
    value_map = {(row["focus_group"], row["condition"]): float(row["state_axis_z_mean"]) for row in condition_rows}
    values = list(value_map.values())
    y_min, y_max = min(values) - 0.25, max(values) + 0.25
    x_map = {"lean": 260, "obese": 520, "obese_enalapril": 780}
    colors = {"DTL": "#ea580c", "Mono_Macro": "#2563eb"}
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    svg_text(svg, 60, 78, "Figure AC. Sample-aware DTL and Mono/Macro state-shift axis", 38, "700")
    svg_text(svg, 60, 118, "State axis is z-scored within focus group and summarized by sample-aware condition means.", 20, "400", "#475569")
    panels = {"DTL": (80, 190, 780, 620), "Mono_Macro": (910, 190, 700, 620)}
    for group, (x0, y0, w, h) in panels.items():
        svg.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" rx="24" fill="white" stroke="#cbd5e1" stroke-width="2"/>')
        svg_text(svg, x0 + 24, y0 + 46, group, 26, "700")
        for condition, dx in x_map.items():
            x = x0 + dx - 120
            svg_text(svg, x, y0 + h - 22, condition, 16, "700", "#334155", "middle")
        points: list[tuple[float, float, str]] = []
        for condition in condition_order():
            value = value_map.get((group, condition), 0.0)
            x = x0 + x_map[condition] - 120
            y = scale(value, y_min, y_max, y0 + h - 85, y0 + 95)
            points.append((x, y, condition))
        for (x1, y1, _), (x2, y2, _) in zip(points, points[1:]):
            svg.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{colors[group]}" stroke-width="5" stroke-linecap="round" opacity="0.72"/>')
        for x, y, condition in points:
            value = value_map.get((group, condition), 0.0)
            svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="18" fill="{colors[group]}" fill-opacity="0.88" stroke="white" stroke-width="3"/>')
            svg_text(svg, x, y - 30, f"{value:.2f}", 16, "700", "#0f172a", "middle")
        call = [row["state_shift_call"] for row in shift_rows if row["focus_group"] == group][0]
        svg_text(svg, x0 + 24, y0 + h - 68, f"Call: {call}", 18, "700", "#475569")
    svg_text(svg, 80, 905, "Interpretation: this supports condition-associated state shift, not lineage causality.", 19, "700", "#475569")
    svg.append("</svg>")
    out_path.write_text("\n".join(svg), encoding="utf-8")


def build_summary(condition_rows: list[dict[str, str]], shift_rows: list[dict[str, str]], gene_presence_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results" / "public_data_gse233078_dtl_monomacro_state_shift"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_r_export(project_root, out_dir)
    cell_rows = read_csv(out_dir / "dtl_monomacro_cell_state_scores.csv")
    gene_presence_rows = read_csv(out_dir / "dtl_monomacro_state_gene_presence.csv")
    sample_rows = build_sample_state_rows(cell_rows)
    condition_rows = build_condition_state_rows(sample_rows)
    shift_rows = build_state_shift_rows(condition_rows)

    write_csv(out_dir / "dtl_monomacro_sample_state_scores.csv", sample_rows, list(sample_rows[0].keys()))
    write_csv(out_dir / "dtl_monomacro_condition_state_summary.csv", condition_rows, list(condition_rows[0].keys()))
    write_csv(out_dir / "dtl_monomacro_state_shift_calls.csv", shift_rows, list(shift_rows[0].keys()))
    build_umap_svg(cell_rows, out_dir / "FigureAB_gse233078_dtl_monomacro_state_map.svg")
    build_state_shift_svg(condition_rows, shift_rows, out_dir / "FigureAC_gse233078_dtl_monomacro_state_shift.svg")
    write_text(out_dir / "gse233078_dtl_monomacro_state_shift_summary.txt", build_summary(condition_rows, shift_rows, gene_presence_rows))
    write_text(out_dir / "status.txt", "GSE233078 DTL/Mono_Macro state-shift completed.\n")

    print(f"GSE233078 DTL/Mono_Macro state-shift written to: {out_dir}")


if __name__ == "__main__":
    main()


