#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import html
import shutil
import subprocess


W = 1720
H = 1040


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


def build_condition_cluster_proportions(metadata_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    totals: dict[str, int] = {}
    counts: dict[tuple[str, str], int] = {}
    for row in metadata_rows:
        condition = row["normalized_condition"]
        cluster = row["cluster"]
        totals[condition] = totals.get(condition, 0) + 1
        counts[(condition, cluster)] = counts.get((condition, cluster), 0) + 1

    rows: list[dict[str, str]] = []
    for (condition, cluster), cell_count in sorted(counts.items()):
        proportion = cell_count / totals[condition] if totals[condition] else 0.0
        rows.append(
            {
                "condition": condition,
                "cluster": cluster,
                "cell_count": str(cell_count),
                "proportion": f"{proportion:.6f}",
            }
        )
    return rows


def classify_cluster_shift(cluster: str, lean: float, obese: float, obese_enalapril: float) -> dict[str, str]:
    if obese - lean >= 0.01 and obese - obese_enalapril >= 0.01:
        shift = "obese_expansion_with_reversal"
    elif lean - obese >= 0.01 and obese_enalapril - obese >= 0.005:
        shift = "obese_loss_with_recovery"
    elif obese - lean >= 0.01:
        shift = "obese_expansion_without_clear_reversal"
    elif lean - obese >= 0.01:
        shift = "obese_loss_without_clear_recovery"
    else:
        shift = "stable_or_minor_shift"
    return {
        "cluster": cluster,
        "lean_proportion": f"{lean:.6f}",
        "obese_proportion": f"{obese:.6f}",
        "obese_enalapril_proportion": f"{obese_enalapril:.6f}",
        "shift_class": shift,
    }


def build_cluster_shift_rows(proportion_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    clusters = sorted({row["cluster"] for row in proportion_rows})
    value_map = {
        (row["condition"], row["cluster"]): float(row["proportion"])
        for row in proportion_rows
    }
    rows: list[dict[str, str]] = []
    for cluster in clusters:
        rows.append(
            classify_cluster_shift(
                cluster,
                value_map.get(("lean", cluster), 0.0),
                value_map.get(("obese", cluster), 0.0),
                value_map.get(("obese_enalapril", cluster), 0.0),
            )
        )
    return rows


def top_priority_rows(shift_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    kept = [
        row for row in shift_rows
        if row["shift_class"] in {"obese_expansion_with_reversal", "obese_loss_with_recovery"}
    ]
    def score(row: dict[str, str]) -> float:
        lean = float(row["lean_proportion"])
        obese = float(row["obese_proportion"])
        ena = float(row["obese_enalapril_proportion"])
        return abs(obese - lean) + abs(obese - ena)
    kept.sort(key=score, reverse=True)
    return kept[:12]


def esc(text: str) -> str:
    return html.escape(str(text), quote=True)


def write_svg_text(svg: list[str], x: int, y: int, text: str, size: int = 24, weight: str = "400", fill: str = "#1f2937") -> None:
    svg.append(
        f'<text x="{x}" y="{y}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}">{esc(text)}</text>'
    )


def shift_color(shift_class: str) -> str:
    return {
        "obese_expansion_with_reversal": "#b91c1c",
        "obese_loss_with_recovery": "#2563eb",
        "obese_expansion_without_clear_reversal": "#ea580c",
        "obese_loss_without_clear_recovery": "#0891b2",
        "stable_or_minor_shift": "#94a3b8",
    }.get(shift_class, "#64748b")


def build_shift_svg(priority_rows: list[dict[str, str]], out_path: Path) -> None:
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
    ]
    write_svg_text(svg, 60, 82, "Figure M. GSE233078 metabolic CKD cell-composition atlas", 40, "700")
    write_svg_text(svg, 60, 126, "Top condition-sensitive clusters from Lean -> Obese -> Obese+enalapril", 22, "400", "#4b5563")

    headers = [("Cluster", 60), ("Lean", 560), ("Obese", 800), ("Obese+enalapril", 1080), ("Shift class", 1360)]
    for text, x in headers:
        write_svg_text(svg, x, 205, text, 24, "700", "#0f172a")

    start_y = 248
    row_h = 62
    for i, row in enumerate(priority_rows):
        y = start_y + i * row_h
        svg.append(f'<rect x="40" y="{y-30}" width="1640" height="46" rx="16" fill="white" stroke="#e2e8f0" stroke-width="2"/>')
        write_svg_text(svg, 60, y, row["cluster"], 22, "700")
        write_svg_text(svg, 560, y, f'{float(row["lean_proportion"]):.3f}', 22, "700")
        write_svg_text(svg, 800, y, f'{float(row["obese_proportion"]):.3f}', 22, "700")
        write_svg_text(svg, 1080, y, f'{float(row["obese_enalapril_proportion"]):.3f}', 22, "700")
        fill = shift_color(row["shift_class"])
        svg.append(f'<rect x="1360" y="{y-24}" width="260" height="34" rx="12" fill="{fill}" opacity="0.92"/>')
        write_svg_text(svg, 1376, y, row["shift_class"], 16, "700", "white")

    write_svg_text(svg, 60, 960, "Reading tip:", 24, "700")
    write_svg_text(svg, 190, 960, "Immune expansion and endothelial/tubular composition shifts are especially useful as metabolic-CKD reference language for future kidney projects.", 20)
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


def build_summary(shift_rows: list[dict[str, str]], priority_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    metadata_rows = read_csv(project_root / "results" / "public_data_gse233078" / "cell_metadata.csv")
    out_dir = project_root / "results" / "public_data_gse233078_metadata_atlas"
    out_dir.mkdir(parents=True, exist_ok=True)

    proportion_rows = build_condition_cluster_proportions(metadata_rows)
    shift_rows = build_cluster_shift_rows(proportion_rows)
    priority_rows = top_priority_rows(shift_rows)

    write_csv(out_dir / "condition_cluster_proportions.csv", proportion_rows, list(proportion_rows[0].keys()))
    write_csv(out_dir / "cluster_shift_summary.csv", shift_rows, list(shift_rows[0].keys()))
    write_csv(out_dir / "priority_cluster_panel.csv", priority_rows, list(priority_rows[0].keys()))
    write_text(out_dir / "gse233078_metadata_atlas_summary.txt", build_summary(shift_rows, priority_rows))
    svg_path = out_dir / "FigureM_gse233078_metadata_atlas.svg"
    build_shift_svg(priority_rows, svg_path)
    export_png(svg_path)
    write_text(out_dir / "status.txt", "GSE233078 metadata atlas completed.\n")

    print(f"GSE233078 metadata atlas written to: {out_dir}")


if __name__ == "__main__":
    main()


