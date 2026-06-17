#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import html
import shutil
import subprocess
import tempfile


W = 1780
H = 960
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


def map_cluster_to_cellchat_group(cluster: str) -> str | None:
    if cluster == "IC-A":
        return "IC-A"
    if cluster in {"Endo", "GEC"}:
        return "Endo_GEC"
    if cluster in {"Mono", "Macro", "Mono/Macro"}:
        return "Mono_Macro"
    if cluster == "DTL":
        return "DTL"
    return None


def build_priority_lr_whitelist() -> list[dict[str, str]]:
    return [
        {
            "module": "vascular_microcirculation",
            "ligand": "Vegfa",
            "receptor": "Kdr",
            "source_group": "IC-A",
            "target_group": "Endo_GEC",
            "priority_tier": "core",
            "note": "Best current vascular axis for CellChat MVP.",
        },
        {
            "module": "vascular_microcirculation",
            "ligand": "Vegfa",
            "receptor": "Flt1",
            "source_group": "IC-A",
            "target_group": "Endo_GEC",
            "priority_tier": "core",
            "note": "Companion VEGF receptor axis supporting the same vascular story.",
        },
        {
            "module": "fibrosis_ecm",
            "ligand": "Tgfb1",
            "receptor": "Tgfbr2",
            "source_group": "Mono_Macro",
            "target_group": "Endo_GEC",
            "priority_tier": "core",
            "note": "Best current fibrosis-remodeling axis for CellChat MVP.",
        },
        {
            "module": "fibrosis_ecm",
            "ligand": "Spp1",
            "receptor": "Cd44",
            "source_group": "DTL",
            "target_group": "Mono_Macro",
            "priority_tier": "hypothesis",
            "note": "Strong signal but should remain hypothesis-level until standard CellChat confirms it.",
        },
    ]


def build_focus_metadata_rows(metadata_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for row in metadata_rows:
        group = map_cluster_to_cellchat_group(row["cluster"])
        if group is None:
            continue
        out_rows.append(
            {
                "cell_id": row["cell_id"],
                "sample_id": row["sample_id"],
                "condition": row["normalized_condition"],
                "original_cluster": row["cluster"],
                "cellchat_group": group,
            }
        )
    return out_rows


def build_group_condition_counts(focus_metadata_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], int] = {}
    for row in focus_metadata_rows:
        key = (row["cellchat_group"], row["condition"])
        grouped[key] = grouped.get(key, 0) + 1
    out_rows: list[dict[str, str]] = []
    for (group, condition), count in sorted(grouped.items()):
        out_rows.append(
            {
                "cellchat_group": group,
                "condition": condition,
                "cell_count": str(count),
            }
        )
    return out_rows


def build_cluster_condition_counts(metadata_rows: list[dict[str, str]]) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = {}
    for row in metadata_rows:
        key = (row["cluster"], row["normalized_condition"])
        counts[key] = counts.get(key, 0) + 1
    return counts


def aggregate_expression_to_cellchat_groups(
    expression_rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    cluster_condition_counts = build_cluster_condition_counts(metadata_rows)
    grouped: dict[tuple[str, str, str], dict[str, float]] = {}
    for row in expression_rows:
        group = map_cluster_to_cellchat_group(row["cluster"])
        if group is None:
            continue
        weight = cluster_condition_counts.get((row["cluster"], row["condition"]), 0)
        if weight == 0:
            continue
        key = (row["gene"], group, row["condition"])
        holder = grouped.setdefault(
            key,
            {"weighted_mean_sum": 0.0, "weighted_pct_sum": 0.0, "cell_count": 0.0},
        )
        holder["weighted_mean_sum"] += float(row["mean_expr"]) * weight
        holder["weighted_pct_sum"] += float(row["pct_expr"]) * weight
        holder["cell_count"] += weight

    out_rows: list[dict[str, str]] = []
    for (gene, group, condition), holder in sorted(grouped.items()):
        total = holder["cell_count"]
        mean_expr = holder["weighted_mean_sum"] / total if total else 0.0
        pct_expr = holder["weighted_pct_sum"] / total if total else 0.0
        out_rows.append(
            {
                "gene": gene,
                "cellchat_group": group,
                "condition": condition,
                "mean_expr": f"{mean_expr:.6f}",
                "pct_expr": f"{pct_expr:.6f}",
                "cell_count": str(int(total)),
            }
        )
    return out_rows


def interaction_score(ligand_mean: float, ligand_pct: float, receptor_mean: float, receptor_pct: float) -> float:
    return (ligand_mean * ligand_pct) * (receptor_mean * receptor_pct)


def classify_interaction_pattern(lean: float, obese: float, obese_enalapril: float) -> str:
    if obese >= lean * 1.2 and obese >= obese_enalapril * 1.15 and (obese - lean) >= 0.01:
        return "obese_up_with_reversal"
    if obese >= lean * 1.2 and (obese - lean) >= 0.01:
        return "obese_up_without_clear_reversal"
    if lean >= obese * 1.2 and obese_enalapril >= obese * 1.15 and (lean - obese) >= 0.01:
        return "obese_down_with_recovery"
    return "stable_or_mixed"


def build_priority_edge_scores(
    aggregated_expression_rows: list[dict[str, str]],
    whitelist_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    value_map = {
        (row["gene"], row["cellchat_group"], row["condition"]): row
        for row in aggregated_expression_rows
    }
    out_rows: list[dict[str, str]] = []
    for pair in whitelist_rows:
        score_map: dict[str, float] = {}
        for condition in ["lean", "obese", "obese_enalapril"]:
            ligand_row = value_map.get((pair["ligand"], pair["source_group"], condition))
            receptor_row = value_map.get((pair["receptor"], pair["target_group"], condition))
            if ligand_row is None or receptor_row is None:
                score_map[condition] = 0.0
                continue
            score_map[condition] = interaction_score(
                float(ligand_row["mean_expr"]),
                float(ligand_row["pct_expr"]),
                float(receptor_row["mean_expr"]),
                float(receptor_row["pct_expr"]),
            )
        out_rows.append(
            {
                "module": pair["module"],
                "ligand": pair["ligand"],
                "receptor": pair["receptor"],
                "source_group": pair["source_group"],
                "target_group": pair["target_group"],
                "priority_tier": pair["priority_tier"],
                "lean_score": f"{score_map['lean']:.6f}",
                "obese_score": f"{score_map['obese']:.6f}",
                "obese_enalapril_score": f"{score_map['obese_enalapril']:.6f}",
                "pattern": classify_interaction_pattern(
                    score_map["lean"], score_map["obese"], score_map["obese_enalapril"]
                ),
            }
        )
    return out_rows


def build_r_script(project_root: Path, out_dir: Path, genes: list[str]) -> str:
    gene_vector = ", ".join(f'"{gene}"' for gene in genes)
    return f"""
library(Matrix)
project_root <- "{project_root.as_posix()}"
out_dir <- "{out_dir.as_posix()}"
meta <- read.csv(file.path(project_root, "results/public_data_gse233078/cell_metadata.csv"), check.names = FALSE)
focus_clusters <- c("IC-A", "Endo", "GEC", "Mono", "Macro", "Mono/Macro", "DTL")
meta <- meta[meta$cluster %in% focus_clusters, ]
con <- gzcon(file(file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_counts_postfilter.rds.gz"), "rb"))
raw1 <- readBin(con, what="raw", n=300000000)
close(con)
raw2 <- memDecompress(raw1, type="gzip")
mat <- unserialize(raw2)
meta <- meta[match(colnames(mat), meta$cell_id), ]
idx <- which(!is.na(meta$cluster) & meta$cluster %in% focus_clusters)
meta <- meta[idx, ]
mat <- mat[, idx, drop = FALSE]
genes <- c({gene_vector})
genes <- genes[genes %in% rownames(mat)]
mat <- mat[genes, , drop = FALSE]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeMM(mat, file.path(out_dir, "cellchat_mvp_focus_whitelist_counts.mtx"))
write.table(rownames(mat), file.path(out_dir, "cellchat_mvp_focus_whitelist_genes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(colnames(mat), file.path(out_dir, "cellchat_mvp_focus_whitelist_cells.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
"""


def run_r_export(project_root: Path, out_dir: Path, genes: list[str]) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_script(project_root, out_dir, genes))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=900)
    finally:
        script_path.unlink(missing_ok=True)


def esc(text: str) -> str:
    return html.escape(str(text), quote=True)


def write_svg_text(
    svg: list[str],
    x: int,
    y: int,
    text: str,
    size: int = 22,
    weight: str = "400",
    fill: str = "#1f2937",
    anchor: str = "start",
) -> None:
    svg.append(
        f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}">{esc(text)}</text>'
    )


def module_color(module: str) -> str:
    return {
        "vascular_microcirculation": "#0f766e",
        "fibrosis_ecm": "#b91c1c",
    }.get(module, "#475569")


def build_prep_svg(
    group_counts: list[dict[str, str]],
    edge_rows: list[dict[str, str]],
    out_path: Path,
) -> None:
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    write_svg_text(svg, 60, 82, "Figure V. GSE233078 CellChat MVP preparation layer", 40, "700")
    write_svg_text(svg, 60, 126, "Prepared metadata, whitelist LR pairs, and lightweight fallback scores before formal CellChat installation.", 22, "400", "#4b5563")

    svg.append('<rect x="40" y="170" width="520" height="720" rx="24" fill="white" stroke="#e2e8f0" stroke-width="2"/>')
    write_svg_text(svg, 60, 215, "A. Focus cell groups", 28, "700")
    headers = [("Group", 60), ("Condition", 260), ("Cells", 470)]
    for text, x in headers:
        write_svg_text(svg, x, 270, text, 20, "700", "#0f172a")
    start_y = 325
    for i, row in enumerate(group_counts):
        y = start_y + i * 56
        svg.append(f'<rect x="55" y="{y-34}" width="490" height="42" rx="12" fill="#f8fafc" stroke="#e2e8f0" stroke-width="1.2"/>')
        write_svg_text(svg, 60, y, row["cellchat_group"], 18, "700")
        write_svg_text(svg, 260, y, row["condition"], 18, "700")
        write_svg_text(svg, 470, y, row["cell_count"], 18, "700")

    svg.append('<rect x="600" y="170" width="1140" height="720" rx="24" fill="white" stroke="#e2e8f0" stroke-width="2"/>')
    write_svg_text(svg, 620, 215, "B. Priority LR edges for MVP", 28, "700")
    headers = [("Axis", 620), ("Route", 870), ("Lean", 1170), ("Obese", 1320), ("Enalapril", 1490), ("Pattern", 1660)]
    for text, x in headers:
        write_svg_text(svg, x, 270, text, 20, "700", "#0f172a")

    start_y = 335
    for i, row in enumerate(edge_rows):
        y = start_y + i * 96
        fill = module_color(row["module"])
        svg.append(f'<rect x="615" y="{y-40}" width="1100" height="62" rx="16" fill="#f8fafc" stroke="#e2e8f0" stroke-width="1.2"/>')
        write_svg_text(svg, 620, y, row["module"], 18, "700", fill)
        write_svg_text(svg, 870, y, f"{row['ligand']}->{row['receptor']} | {row['source_group']}->{row['target_group']}", 17, "700")
        write_svg_text(svg, 1185, y, row["lean_score"], 17, "700")
        write_svg_text(svg, 1335, y, row["obese_score"], 17, "700")
        write_svg_text(svg, 1510, y, row["obese_enalapril_score"], 17, "700")
        write_svg_text(svg, 1660, y, row["pattern"], 15, "700", "#0f766e")

    write_svg_text(svg, 620, 820, "Prepared files:", 22, "700")
    write_svg_text(svg, 770, 820, "focus metadata + whitelist LR table + whitelist count matrix + lightweight edge-score fallback", 18, "400", "#475569")
    svg.append("</svg>")
    out_path.write_text("\n".join(svg), encoding="utf-8")


def build_summary(
    group_counts: list[dict[str, str]],
    whitelist_rows: list[dict[str, str]],
    edge_rows: list[dict[str, str]],
) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    metadata_rows = read_csv(project_root / "results" / "public_data_gse233078" / "cell_metadata.csv")
    expression_rows = read_csv(project_root / "results" / "public_data_gse233078_interaction_map" / "lr_gene_cluster_condition_means.csv")
    whitelist_rows = build_priority_lr_whitelist()
    focus_metadata_rows = build_focus_metadata_rows(metadata_rows)
    group_count_rows = build_group_condition_counts(focus_metadata_rows)
    aggregated_rows = aggregate_expression_to_cellchat_groups(expression_rows, metadata_rows)
    edge_rows = build_priority_edge_scores(aggregated_rows, whitelist_rows)

    out_dir = project_root / "results" / "public_data_gse233078_cellchat_mvp_prep"
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = sorted({row["ligand"] for row in whitelist_rows} | {row["receptor"] for row in whitelist_rows})
    run_r_export(project_root, out_dir, genes)

    write_csv(out_dir / "cellchat_mvp_focus_metadata.csv", focus_metadata_rows, list(focus_metadata_rows[0].keys()))
    write_csv(out_dir / "cellchat_mvp_group_condition_counts.csv", group_count_rows, list(group_count_rows[0].keys()))
    write_csv(out_dir / "cellchat_mvp_priority_lr_whitelist.csv", whitelist_rows, list(whitelist_rows[0].keys()))
    write_csv(out_dir / "cellchat_mvp_whitelist_group_expression.csv", aggregated_rows, list(aggregated_rows[0].keys()))
    write_csv(out_dir / "cellchat_mvp_priority_edge_scores.csv", edge_rows, list(edge_rows[0].keys()))
    write_text(out_dir / "gse233078_cellchat_mvp_prep_summary.txt", build_summary(group_count_rows, whitelist_rows, edge_rows))
    svg_path = out_dir / "FigureV_gse233078_cellchat_mvp_prep.svg"
    build_prep_svg(group_count_rows, edge_rows, svg_path)
    write_text(out_dir / "status.txt", "GSE233078 CellChat MVP preparation completed.\n")

    print(f"GSE233078 CellChat MVP preparation written to: {out_dir}")


if __name__ == "__main__":
    main()


