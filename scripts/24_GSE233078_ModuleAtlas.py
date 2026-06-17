#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import html
import shutil
import subprocess


W = 1700
H = 980


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


def classify_expression_pattern(lean: float, obese: float, obese_enalapril: float) -> str:
    if obese - lean >= 0.10 and obese - obese_enalapril >= 0.05:
        return "obese_up_with_reversal"
    if lean - obese >= 0.10 and obese_enalapril - obese >= 0.05:
        return "obese_down_with_recovery"
    if obese - lean >= 0.10:
        return "obese_up_without_clear_reversal"
    if lean - obese >= 0.10:
        return "obese_down_without_clear_recovery"
    return "stable_or_mixed"


def normalize_gene_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    maxima: dict[str, float] = {}
    for row in rows:
        gene = row["gene"]
        value = float(row["mean_expr"])
        maxima[gene] = max(maxima.get(gene, 0.0), value)

    out_rows: list[dict[str, str]] = []
    for row in rows:
        gene = row["gene"]
        max_value = maxima.get(gene, 0.0)
        normalized = float(row["mean_expr"]) / max_value if max_value > 0 else 0.0
        out_rows.append({**row, "normalized_expr": f"{normalized:.6f}"})
    return out_rows


def build_module_cluster_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[float]] = {}
    for row in rows:
        key = (row["module"], row["cluster"], row["condition"])
        grouped.setdefault(key, []).append(float(row["normalized_expr"]))
    out_rows: list[dict[str, str]] = []
    for (module, cluster, condition), values in sorted(grouped.items()):
        score = sum(values) / len(values) if values else 0.0
        out_rows.append(
            {
                "module": module,
                "cluster": cluster,
                "condition": condition,
                "module_score": f"{score:.6f}",
            }
        )
    return out_rows


def build_priority_module_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    modules = sorted({row["module"] for row in rows})
    value_map = {
        (row["module"], row["cluster"], row["condition"]): float(row["module_score"])
        for row in rows
    }
    clusters = sorted({row["cluster"] for row in rows})
    out_rows: list[dict[str, str]] = []
    for module in modules:
        priority_cluster = max(
            clusters,
            key=lambda cluster: (
                value_map.get((module, cluster, "obese"), 0.0),
                value_map.get((module, cluster, "obese_enalapril"), 0.0),
                value_map.get((module, cluster, "lean"), 0.0),
            ),
        )
        lean = value_map.get((module, priority_cluster, "lean"), 0.0)
        obese = value_map.get((module, priority_cluster, "obese"), 0.0)
        ena = value_map.get((module, priority_cluster, "obese_enalapril"), 0.0)
        out_rows.append(
            {
                "module": module,
                "priority_cluster": priority_cluster,
                "lean_module_score": f"{lean:.4f}",
                "obese_module_score": f"{obese:.4f}",
                "obese_enalapril_module_score": f"{ena:.4f}",
                "pattern": classify_expression_pattern(lean, obese, ena),
            }
        )
    return out_rows


def esc(text: str) -> str:
    return html.escape(str(text), quote=True)


def write_svg_text(svg: list[str], x: int, y: int, text: str, size: int = 22, weight: str = "400", fill: str = "#1f2937") -> None:
    svg.append(
        f'<text x="{x}" y="{y}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}">{esc(text)}</text>'
    )


def pattern_color(pattern: str) -> str:
    return {
        "obese_up_with_reversal": "#b91c1c",
        "obese_down_with_recovery": "#2563eb",
        "obese_up_without_clear_reversal": "#ea580c",
        "obese_down_without_clear_recovery": "#0891b2",
        "stable_or_mixed": "#94a3b8",
    }.get(pattern, "#64748b")


def build_module_svg(summary_rows: list[dict[str, str]], out_path: Path) -> None:
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    write_svg_text(svg, 60, 82, "Figure O. GSE233078 module-level cell atlas", 40, "700")
    write_svg_text(svg, 60, 126, "Module-normalized priority clusters in metabolic CKD single-cell reference", 22, "400", "#4b5563")

    headers = [("Module", 60), ("Priority cluster", 430), ("Lean", 850), ("Obese", 1010), ("Enalapril", 1180), ("Pattern", 1360)]
    for text, x in headers:
        write_svg_text(svg, x, 205, text, 24, "700", "#0f172a")

    start_y = 260
    row_h = 115
    for i, row in enumerate(summary_rows):
        y = start_y + i * row_h
        svg.append(f'<rect x="40" y="{y-46}" width="1620" height="86" rx="18" fill="white" stroke="#e2e8f0" stroke-width="2"/>')
        write_svg_text(svg, 60, y, row["module"], 24, "700")
        write_svg_text(svg, 430, y, row["priority_cluster"], 24, "700")
        write_svg_text(svg, 860, y, row["lean_module_score"], 22, "700")
        write_svg_text(svg, 1020, y, row["obese_module_score"], 22, "700")
        write_svg_text(svg, 1190, y, row["obese_enalapril_module_score"], 22, "700")
        fill = pattern_color(row["pattern"])
        svg.append(f'<rect x="1360" y="{y-28}" width="260" height="40" rx="12" fill="{fill}" opacity="0.92"/>')
        write_svg_text(svg, 1378, y, row["pattern"], 16, "700", "white")

    write_svg_text(svg, 60, 925, "Reading tip:", 24, "700")
    write_svg_text(svg, 190, 925, "This figure compresses gene-level information into module-level priority cell sources for the metabolic CKD branch.", 20)
    svg.append("</svg>")
    out_path.write_text("\n".join(svg), encoding="utf-8")


def export_png(svg_path: Path) -> None:
    edge = Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
    if not edge.exists():
        return
    png_path = svg_path.with_suffix(".png")
    subprocess.run(
        [
            str(edge),
            "--headless",
            "--disable-gpu",
            f"--screenshot={png_path}",
            f"--window-size={W},{H}",
            str(svg_path),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def build_summary(summary_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results" / "public_data_gse233078_module_atlas"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_csv(project_root / "results" / "public_data_gse233078_target_gene_atlas" / "target_gene_cluster_condition_means.csv")
    normalized_rows = normalize_gene_rows(rows)
    module_rows = build_module_cluster_scores(normalized_rows)
    summary_rows = build_priority_module_summary(module_rows)

    write_csv(out_dir / "target_gene_cluster_condition_means_normalized.csv", normalized_rows, list(normalized_rows[0].keys()))
    write_csv(out_dir / "module_cluster_condition_scores.csv", module_rows, list(module_rows[0].keys()))
    write_csv(out_dir / "module_priority_clusters.csv", summary_rows, list(summary_rows[0].keys()))
    write_text(out_dir / "gse233078_module_atlas_summary.txt", build_summary(summary_rows))
    svg_path = out_dir / "FigureO_gse233078_module_atlas.svg"
    build_module_svg(summary_rows, svg_path)
    export_png(svg_path)
    write_text(out_dir / "status.txt", "GSE233078 module atlas completed.\n")

    print(f"GSE233078 module atlas written to: {out_dir}")


if __name__ == "__main__":
    main()


