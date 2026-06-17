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


W = 1650
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


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def log2fc(reference_values: list[float], contrast_values: list[float]) -> float:
    return math.log2((mean(contrast_values) + 0.1) / (mean(reference_values) + 0.1))


def contrast_rows() -> list[dict[str, str]]:
    return [
        {
            "contrast": "doca_vs_control",
            "reference_group": "control",
            "contrast_group": "doca",
            "contrast_role": "disease_induction",
        },
        {
            "contrast": "finerenone_vs_doca",
            "reference_group": "doca",
            "contrast_group": "finerenone",
            "contrast_role": "treatment_reversal",
        },
        {
            "contrast": "spironolactone_vs_doca",
            "reference_group": "doca",
            "contrast_group": "spironolactone",
            "contrast_role": "treatment_reversal",
        },
    ]


def load_group_to_samples(metadata_rows: list[dict[str, str]]) -> dict[str, list[str]]:
    group_map: dict[str, list[str]] = {}
    for row in metadata_rows:
        group_map.setdefault(row["normalized_group"], []).append(row["sample_id"])
    return group_map


def parse_target_tpm(
    tpm_path: Path,
    target_map: list[dict[str, str]],
) -> dict[str, dict[str, float]]:
    id_to_target = {row["rat_ensembl_id"].upper(): row["gene"] for row in target_map}
    found: dict[str, dict[str, float]] = {}
    with gzip.open(tpm_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col_index = {col: idx for idx, col in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            gene_id = parts[0].strip().upper()
            gene_symbol = id_to_target.get(gene_id)
            if not gene_symbol:
                continue
            values: dict[str, float] = {}
            for sample_id, idx in col_index.items():
                if sample_id == "gene":
                    continue
                try:
                    values[sample_id] = float(parts[idx])
                except (IndexError, ValueError):
                    values[sample_id] = 0.0
            found[gene_symbol] = values
    return found


def build_gene_contrast_rows(project_root: Path) -> list[dict[str, str]]:
    metadata_rows = read_csv(project_root / "results" / "public_data_gse183841" / "sample_metadata.csv")
    target_map = read_csv(project_root / "results" / "public_data_gse183841" / "target_gene_ensembl_map.csv")
    group_map = load_group_to_samples(metadata_rows)
    tpm_path = project_root / "data" / "public_datasets" / "GSE183841" / "processed" / "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz"
    gene_values = parse_target_tpm(tpm_path, target_map)
    gene_to_module = {row["gene"]: row["module"] for row in target_map}

    output_rows: list[dict[str, str]] = []
    for contrast in contrast_rows():
        reference_samples = group_map.get(contrast["reference_group"], [])
        contrast_samples = group_map.get(contrast["contrast_group"], [])
        for gene, sample_values in gene_values.items():
            reference_values = [sample_values[sid] for sid in reference_samples if sid in sample_values]
            contrast_values = [sample_values[sid] for sid in contrast_samples if sid in sample_values]
            fc = log2fc(reference_values, contrast_values)
            output_rows.append(
                {
                    "contrast": contrast["contrast"],
                    "contrast_role": contrast["contrast_role"],
                    "reference_group": contrast["reference_group"],
                    "contrast_group": contrast["contrast_group"],
                    "module": gene_to_module[gene],
                    "gene": gene,
                    "reference_mean_tpm": f"{mean(reference_values):.4f}",
                    "contrast_mean_tpm": f"{mean(contrast_values):.4f}",
                    "log2fc_contrast_vs_reference": f"{fc:.4f}",
                    "direction": "up" if fc > 0.25 else "down" if fc < -0.25 else "flat",
                }
            )
    return output_rows


def build_module_summary_rows(gene_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[float]] = {}
    contrast_role_map: dict[str, str] = {}
    for row in gene_rows:
        key = (row["contrast"], row["module"])
        grouped.setdefault(key, []).append(float(row["log2fc_contrast_vs_reference"]))
        contrast_role_map[row["contrast"]] = row.get("contrast_role", "")

    rows: list[dict[str, str]] = []
    for (contrast, module), values in sorted(grouped.items()):
        avg = mean(values)
        rows.append(
            {
                "contrast": contrast,
                "contrast_role": contrast_role_map.get(contrast, ""),
                "module": module,
                "mean_log2fc": f"{avg:.4f}",
                "direction": "up" if avg > 0.25 else "down" if avg < -0.25 else "flat",
                "gene_count": str(len(values)),
            }
        )
    return rows


def build_rescue_interpretation_rows(module_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    disease_direction_map = {
        row["module"]: float(row["mean_log2fc"])
        for row in module_rows
        if row["contrast"] == "doca_vs_control"
    }
    rows: list[dict[str, str]] = []
    for row in module_rows:
        value = float(row["mean_log2fc"])
        module = row["module"]
        contrast = row["contrast"]
        if contrast == "doca_vs_control":
            if value > 0.25:
                call = "disease_shift_up"
            elif value < -0.25:
                call = "disease_shift_down"
            else:
                call = "not_strong"
        else:
            disease_direction = disease_direction_map.get(module, 0.0)
            if abs(value) <= 0.25:
                call = "not_strong"
            elif disease_direction > 0.25 and value < -0.25:
                call = "normalization_consistent"
            elif disease_direction < -0.25 and value > 0.25:
                call = "normalization_consistent"
            elif disease_direction > 0.25 and value > 0.25:
                call = "same_direction_as_disease"
            elif disease_direction < -0.25 and value < -0.25:
                call = "same_direction_as_disease"
            else:
                call = "mixed_or_context_dependent"
        rows.append(
            {
                "contrast": contrast,
                "module": module,
                "mean_log2fc": row["mean_log2fc"],
                "project_call": call,
            }
        )
    return rows


def esc(text: str) -> str:
    return html.escape(str(text), quote=True)


def write_svg_text(
    svg: list[str],
    x: int,
    y: int,
    text: str,
    size: int = 26,
    weight: str = "400",
    fill: str = "#1f2937",
) -> None:
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
    contrasts = ["doca_vs_control", "finerenone_vs_doca", "spironolactone_vs_doca"]
    value_map = {(row["contrast"], row["module"]): float(row["mean_log2fc"]) for row in module_rows}

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    write_svg_text(svg, 70, 85, "Figure K. GSE183841 module-direction heatmap", 40, "700")
    write_svg_text(svg, 70, 130, "DOCA disease induction plus finerenone/spironolactone reversal projected onto kidney injury modules", 22, "400", "#4b5563")

    start_x = 420
    start_y = 220
    cell_w = 300
    cell_h = 110
    for j, contrast in enumerate(contrasts):
        write_svg_text(svg, start_x + j * cell_w + 25, start_y - 35, contrast, 24, "700")
    for i, module in enumerate(modules):
        y = start_y + i * cell_h
        write_svg_text(svg, 70, y + 68, module, 24, "700", "#334155")
        for j, contrast in enumerate(contrasts):
            value = value_map.get((contrast, module), 0.0)
            x = start_x + j * cell_w
            svg.append(
                f'<rect x="{x}" y="{y}" width="{cell_w-18}" height="{cell_h-18}" rx="16" '
                f'fill="{heat_color(value)}" stroke="white" stroke-width="4"/>'
            )
            write_svg_text(
                svg,
                x + 78,
                y + 58,
                f"{value:.2f}",
                30,
                "700",
                "white" if abs(value) >= 0.5 else "#111827",
            )

    write_svg_text(svg, 70, 815, "Interpretation:", 24, "700")
    write_svg_text(svg, 220, 815, "Use DOCA vs Control as disease prior, and Finerenone/Spironolactone vs DOCA as public reversal references for pathway plausibility.", 21)
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


def build_summary(module_rows: list[dict[str, str]], rescue_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results" / "public_data_gse183841_module_direction"
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_rows = build_gene_contrast_rows(project_root)
    module_rows = build_module_summary_rows(gene_rows)
    rescue_rows = build_rescue_interpretation_rows(module_rows)

    write_csv(out_dir / "gse183841_candidate_gene_log2fc.csv", gene_rows, list(gene_rows[0].keys()))
    write_csv(out_dir / "gse183841_module_direction_summary.csv", module_rows, list(module_rows[0].keys()))
    write_csv(out_dir / "gse183841_rescue_interpretation.csv", rescue_rows, list(rescue_rows[0].keys()))
    write_text(out_dir / "gse183841_module_direction_summary.txt", build_summary(module_rows, rescue_rows))

    svg_path = out_dir / "FigureK_gse183841_module_direction_heatmap.svg"
    build_module_heatmap(module_rows, svg_path)
    export_png(svg_path)
    write_text(out_dir / "status.txt", "GSE183841 module direction analysis completed.\n")

    print(f"GSE183841 module direction results written to: {out_dir}")


if __name__ == "__main__":
    main()


