#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import gzip
import html
import math
import shutil
import subprocess


W = 1500
H = 950


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


def detect_groups(header: list[str]) -> dict[str, list[str]]:
    sample_cols = header[1:]
    control = [col for col in sample_cols if "control" in col.lower() or "sham" in col.lower()]
    disease = [
        col
        for col in sample_cols
        if "adenine" in col.lower() or "uuo" in col.lower() or "uuo2w" in col.lower()
    ]
    return {"control": control, "disease": disease}


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def log2fc(control_values: list[float], disease_values: list[float]) -> float:
    return math.log2((mean(disease_values) + 0.1) / (mean(control_values) + 0.1))


def parse_target_fpkm(path: Path, target_genes: set[str]) -> tuple[dict[str, dict[str, float]], dict[str, list[str]]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        groups = detect_groups(header)
        col_index = {col: idx for idx, col in enumerate(header)}
        found: dict[str, dict[str, float]] = {}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            gene = parts[0].strip().upper()
            if gene not in target_genes:
                continue
            values: dict[str, float] = {}
            for col in groups["control"] + groups["disease"]:
                idx = col_index[col]
                try:
                    values[col] = float(parts[idx])
                except (IndexError, ValueError):
                    values[col] = 0.0
            found[gene] = values
    return found, groups


def target_gene_rows() -> list[dict[str, str]]:
    return [
        {"module": "vascular_microcirculation", "gene": "NOS3"},
        {"module": "vascular_microcirculation", "gene": "VEGFA"},
        {"module": "vascular_microcirculation", "gene": "HIF1A"},
        {"module": "vascular_microcirculation", "gene": "AKT1"},
        {"module": "fibrosis_ecm", "gene": "TGFB1"},
        {"module": "fibrosis_ecm", "gene": "COL1A1"},
        {"module": "fibrosis_ecm", "gene": "FN1"},
        {"module": "fibrosis_ecm", "gene": "ACTA2"},
        {"module": "inflammation", "gene": "TLR4"},
        {"module": "inflammation", "gene": "NLRP3"},
        {"module": "metabolic_stress", "gene": "HAO2"},
    ]


def model_file_rows(data_dir: Path) -> list[dict[str, str]]:
    return [
        {
            "model_arm": "adenine",
            "file": str(data_dir / "processed" / "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz"),
        },
        {
            "model_arm": "uuo",
            "file": str(data_dir / "processed" / "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz"),
        },
    ]


def build_gene_log2fc_rows(data_dir: Path) -> list[dict[str, str]]:
    targets = target_gene_rows()
    target_set = {row["gene"].upper() for row in targets}
    module_map = {row["gene"].upper(): row["module"] for row in targets}
    output_rows: list[dict[str, str]] = []
    for model in model_file_rows(data_dir):
        found, groups = parse_target_fpkm(Path(model["file"]), target_set)
        for gene in sorted(target_set):
            values = found.get(gene, {})
            control_values = [values[col] for col in groups["control"] if col in values]
            disease_values = [values[col] for col in groups["disease"] if col in values]
            fc = log2fc(control_values, disease_values) if control_values and disease_values else 0.0
            output_rows.append(
                {
                    "model_arm": model["model_arm"],
                    "module": module_map[gene],
                    "gene": gene,
                    "control_mean_fpkm": f"{mean(control_values):.4f}",
                    "disease_mean_fpkm": f"{mean(disease_values):.4f}",
                    "log2fc_disease_vs_control": f"{fc:.4f}",
                    "direction": "up" if fc > 0.25 else "down" if fc < -0.25 else "flat",
                }
            )
    return output_rows


def build_module_summary_rows(gene_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[float]] = {}
    for row in gene_rows:
        key = (row["model_arm"], row["module"])
        grouped.setdefault(key, []).append(float(row["log2fc_disease_vs_control"]))
    rows: list[dict[str, str]] = []
    for (model_arm, module), values in sorted(grouped.items()):
        avg = mean(values)
        rows.append(
            {
                "model_arm": model_arm,
                "module": module,
                "mean_log2fc": f"{avg:.4f}",
                "direction": "up" if avg > 0.25 else "down" if avg < -0.25 else "flat",
                "gene_count": str(len(values)),
            }
        )
    return rows


def esc(text: str) -> str:
    return html.escape(str(text), quote=True)


def write_svg_text(svg: list[str], x: int, y: int, text: str, size: int = 26, weight: str = "400", fill: str = "#1f2937") -> None:
    svg.append(
        f'<text x="{x}" y="{y}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}">{esc(text)}</text>'
    )


def heat_color(value: float) -> str:
    if value >= 1.0:
        return "#991b1b"
    if value >= 0.5:
        return "#ef4444"
    if value >= 0.25:
        return "#fca5a5"
    if value <= -1.0:
        return "#1d4ed8"
    if value <= -0.5:
        return "#3b82f6"
    if value <= -0.25:
        return "#bfdbfe"
    return "#e5e7eb"


def build_module_heatmap(module_rows: list[dict[str, str]], out_path: Path) -> None:
    modules = ["vascular_microcirculation", "fibrosis_ecm", "inflammation", "metabolic_stress"]
    models = ["adenine", "uuo"]
    value_map = {(row["model_arm"], row["module"]): float(row["mean_log2fc"]) for row in module_rows}
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    write_svg_text(svg, 70, 85, "GSE216376 cross-model module direction", 40, "700")
    write_svg_text(svg, 70, 130, "Adenine and UUO public FPKM matrices mapped to kidney injury candidate modules", 22, "400", "#4b5563")

    start_x = 390
    start_y = 220
    cell_w = 280
    cell_h = 110
    for j, model in enumerate(models):
        write_svg_text(svg, start_x + j * cell_w + 80, start_y - 35, model, 28, "700")
    for i, module in enumerate(modules):
        y = start_y + i * cell_h
        write_svg_text(svg, 70, y + 68, module, 24, "700", "#334155")
        for j, model in enumerate(models):
            value = value_map.get((model, module), 0.0)
            x = start_x + j * cell_w
            svg.append(f'<rect x="{x}" y="{y}" width="{cell_w-18}" height="{cell_h-18}" rx="16" fill="{heat_color(value)}" stroke="white" stroke-width="4"/>')
            write_svg_text(svg, x + 78, y + 58, f"{value:.2f}", 30, "700", "white" if abs(value) >= 0.5 else "#111827")
    write_svg_text(svg, 70, 790, "Interpretation:", 24, "700")
    write_svg_text(svg, 220, 790, "Use this as a public-data prior. Stronger local conclusions require your own local bulk RNA or wet-lab validation.", 21)
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


def build_summary(module_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data" / "public_datasets" / "GSE216376"
    out_dir = project_root / "results" / "public_data_gse216376_module_direction"
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_rows = build_gene_log2fc_rows(data_dir)
    module_rows = build_module_summary_rows(gene_rows)
    write_csv(out_dir / "gse216376_candidate_gene_log2fc.csv", gene_rows, list(gene_rows[0].keys()))
    write_csv(out_dir / "gse216376_module_direction_summary.csv", module_rows, list(module_rows[0].keys()))
    write_text(out_dir / "gse216376_module_direction_summary.txt", build_summary(module_rows))
    svg_path = out_dir / "FigureI_gse216376_module_direction_heatmap.svg"
    build_module_heatmap(module_rows, svg_path)
    export_png(svg_path)
    write_text(out_dir / "status.txt", "GSE216376 module direction analysis completed.\n")

    print(f"GSE216376 module direction results written to: {out_dir}")


if __name__ == "__main__":
    main()


