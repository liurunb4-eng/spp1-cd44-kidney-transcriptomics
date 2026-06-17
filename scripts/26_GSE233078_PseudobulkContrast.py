#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import html
import shutil
import subprocess
import tempfile


W = 1720
H = 1060
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


def build_target_gene_rows() -> list[dict[str, str]]:
    return [
        {"module": "vascular_microcirculation", "gene": "Nos3"},
        {"module": "vascular_microcirculation", "gene": "Vegfa"},
        {"module": "vascular_microcirculation", "gene": "Hif1a"},
        {"module": "vascular_microcirculation", "gene": "Akt1"},
        {"module": "fibrosis_ecm", "gene": "Tgfb1"},
        {"module": "fibrosis_ecm", "gene": "Col1a1"},
        {"module": "fibrosis_ecm", "gene": "Fn1"},
        {"module": "fibrosis_ecm", "gene": "Acta2"},
        {"module": "inflammation", "gene": "Tlr4"},
        {"module": "inflammation", "gene": "Nlrp3"},
        {"module": "metabolic_stress", "gene": "Hao2"},
    ]


def compute_group_stats(values: list[float]) -> dict[str, float | int]:
    if not values:
        return {"n": 0, "mean": 0.0}
    return {"n": len(values), "mean": round(sum(values) / len(values), 10)}


def direction(value: float, threshold: float = 0.05) -> str:
    if value > threshold:
        return "up"
    if value < -threshold:
        return "down"
    return "flat"


def build_contrast_rows(sample_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in sample_rows:
        grouped.setdefault((row["module"], row["cluster"]), []).append(row)

    out_rows: list[dict[str, str]] = []
    for (module, cluster), rows in sorted(grouped.items()):
        lean_values = [float(row["module_score"]) for row in rows if row["condition"] == "lean"]
        obese_values = [float(row["module_score"]) for row in rows if row["condition"] == "obese"]
        ena_values = [float(row["module_score"]) for row in rows if row["condition"] == "obese_enalapril"]
        lean_stats = compute_group_stats(lean_values)
        obese_stats = compute_group_stats(obese_values)
        ena_stats = compute_group_stats(ena_values)
        obese_vs_lean = obese_stats["mean"] - lean_stats["mean"]
        ena_vs_obese = ena_stats["mean"] - obese_stats["mean"]
        if direction(obese_vs_lean) == "up" and direction(ena_vs_obese) == "down":
            call = "reversal_consistent"
        elif direction(obese_vs_lean) == "down" and direction(ena_vs_obese) == "up":
            call = "recovery_consistent"
        else:
            call = "mixed_or_flat"
        out_rows.append(
            {
                "module": module,
                "cluster": cluster,
                "lean_n": str(lean_stats["n"]),
                "obese_n": str(obese_stats["n"]),
                "obese_enalapril_n": str(ena_stats["n"]),
                "lean_mean": f"{lean_stats['mean']:.4f}",
                "obese_mean": f"{obese_stats['mean']:.4f}",
                "obese_enalapril_mean": f"{ena_stats['mean']:.4f}",
                "obese_vs_lean_delta": f"{obese_vs_lean:.4f}",
                "enalapril_vs_obese_delta": f"{ena_vs_obese:.4f}",
                "obese_vs_lean_direction": direction(obese_vs_lean),
                "enalapril_vs_obese_direction": direction(ena_vs_obese),
                "contrast_call": call,
            }
        )
    return out_rows


def select_priority_contrasts(contrast_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    modules = sorted({row["module"] for row in contrast_rows})
    out_rows: list[dict[str, str]] = []
    for module in modules:
        candidates = [row for row in contrast_rows if row["module"] == module]
        candidates.sort(
            key=lambda row: (
                1 if row["contrast_call"] == "reversal_consistent" else 0,
                abs(float(row["obese_vs_lean_delta"])) + abs(float(row["enalapril_vs_obese_delta"])),
                float(row["obese_mean"]),
            ),
            reverse=True,
        )
        out_rows.append(candidates[0])
    return out_rows


def build_r_script(project_root: Path, out_csv: Path) -> str:
    genes = [row["gene"] for row in build_target_gene_rows()]
    gene_vector = ", ".join(f'"{gene}"' for gene in genes)
    return f"""
library(Matrix)
project_root <- "{project_root.as_posix()}"
meta <- read.csv(file.path(project_root, "results/public_data_gse233078/cell_metadata.csv"), check.names = FALSE)
con <- gzcon(file(file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_counts_postfilter.rds.gz"), "rb"))
raw1 <- readBin(con, what="raw", n=300000000)
close(con)
raw2 <- memDecompress(raw1, type="gzip")
mat <- unserialize(raw2)
meta <- meta[match(colnames(mat), meta$cell_id), ]
genes <- c({gene_vector})
genes <- genes[genes %in% rownames(mat)]
all_rows <- data.frame()
for (sample_id in unique(meta$sample_id)) {{
  sample_meta <- meta[meta$sample_id == sample_id, ]
  for (cluster in unique(sample_meta$cluster)) {{
    idx <- which(meta$sample_id == sample_id & meta$cluster == cluster)
    if (length(idx) < 20) next
    cond <- unique(meta$normalized_condition[idx])[1]
    for (gene in genes) {{
      vals <- as.numeric(mat[gene, idx])
      all_rows <- rbind(all_rows, data.frame(
        sample_id = sample_id,
        condition = cond,
        cluster = cluster,
        gene = gene,
        mean_expr = mean(vals),
        cell_count = length(idx)
      ))
    }}
  }}
}}
write.csv(all_rows, "{out_csv.as_posix()}", row.names = FALSE)
"""


def run_r_export(project_root: Path, out_csv: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_script(project_root, out_csv))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=900)
    finally:
        script_path.unlink(missing_ok=True)


def normalize_gene_sample_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    maxima: dict[str, float] = {}
    for row in rows:
        gene = row["gene"]
        value = float(row["mean_expr"])
        maxima[gene] = max(maxima.get(gene, 0.0), value)
    out: list[dict[str, str]] = []
    for row in rows:
        gene = row["gene"]
        normalized = float(row["mean_expr"]) / maxima[gene] if maxima.get(gene, 0.0) > 0 else 0.0
        out.append({**row, "normalized_expr": f"{normalized:.6f}"})
    return out


def build_sample_module_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str, str], list[float]] = {}
    module_map = {row["gene"]: row["module"] for row in build_target_gene_rows()}
    for row in rows:
        key = (module_map[row["gene"]], row["sample_id"], row["condition"], row["cluster"])
        grouped.setdefault(key, []).append(float(row["normalized_expr"]))
    out_rows: list[dict[str, str]] = []
    for (module, sample_id, condition, cluster), values in sorted(grouped.items()):
        out_rows.append(
            {
                "module": module,
                "sample_id": sample_id,
                "condition": condition,
                "cluster": cluster,
                "module_score": f"{(sum(values) / len(values)):.6f}",
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


def contrast_color(call: str) -> str:
    return {
        "reversal_consistent": "#047857",
        "recovery_consistent": "#2563eb",
        "mixed_or_flat": "#94a3b8",
    }.get(call, "#64748b")


def build_contrast_svg(rows: list[dict[str, str]], out_path: Path) -> None:
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    write_svg_text(svg, 60, 82, "Figure Q. GSE233078 sample-aware pseudobulk contrast", 40, "700")
    write_svg_text(svg, 60, 126, "Module-cluster contrasts across lean, obese, and obese+enalapril using sample-aware aggregation", 22, "400", "#4b5563")

    headers = [("Module", 60), ("Cluster", 310), ("Obese-Lean", 620), ("Enalapril-Obese", 860), ("Call", 1130)]
    for text, x in headers:
        write_svg_text(svg, x, 205, text, 24, "700", "#0f172a")

    start_y = 260
    row_h = 125
    for i, row in enumerate(rows):
        y = start_y + i * row_h
        svg.append(f'<rect x="40" y="{y-46}" width="1580" height="86" rx="18" fill="white" stroke="#e2e8f0" stroke-width="2"/>')
        write_svg_text(svg, 60, y, row["module"], 24, "700")
        write_svg_text(svg, 310, y, row["cluster"], 24, "700")
        write_svg_text(svg, 650, y, row["obese_vs_lean_delta"], 22, "700")
        write_svg_text(svg, 910, y, row["enalapril_vs_obese_delta"], 22, "700")
        fill = contrast_color(row["contrast_call"])
        svg.append(f'<rect x="1120" y="{y-28}" width="260" height="40" rx="12" fill="{fill}" opacity="0.92"/>')
        write_svg_text(svg, 1138, y, row["contrast_call"], 16, "700", "white")

    write_svg_text(svg, 60, 940, "Reading tip:", 24, "700")
    write_svg_text(svg, 190, 940, "This is closer to a manuscript-style comparison because it aggregates signals at the sample level before comparing conditions.", 20)
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


def build_summary(rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results" / "public_data_gse233078_pseudobulk"
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_csv = out_dir / "target_gene_sample_cluster_means.csv"
    run_r_export(project_root, raw_csv)
    raw_rows = read_csv(raw_csv)
    normalized_rows = normalize_gene_sample_rows(raw_rows)
    sample_module_rows = build_sample_module_rows(normalized_rows)
    contrast_rows = build_contrast_rows(sample_module_rows)
    priority_rows = select_priority_contrasts(contrast_rows)

    write_csv(raw_csv, normalized_rows, list(normalized_rows[0].keys()))
    write_csv(out_dir / "sample_module_scores.csv", sample_module_rows, list(sample_module_rows[0].keys()))
    write_csv(out_dir / "module_cluster_contrast_summary.csv", contrast_rows, list(contrast_rows[0].keys()))
    write_csv(out_dir / "module_cluster_contrast_priority.csv", priority_rows, list(priority_rows[0].keys()))
    write_text(out_dir / "gse233078_pseudobulk_summary.txt", build_summary(priority_rows))
    svg_path = out_dir / "FigureQ_gse233078_sample_aware_pseudobulk.svg"
    build_contrast_svg(priority_rows, svg_path)
    export_png(svg_path)
    write_text(out_dir / "status.txt", "GSE233078 sample-aware pseudobulk completed.\n")

    print(f"GSE233078 sample-aware pseudobulk written to: {out_dir}")


if __name__ == "__main__":
    main()


