#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
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


def build_condition_order() -> list[str]:
    return ["lean", "obese", "obese_enalapril"]


def build_priority_axis_rows() -> list[dict[str, str]]:
    return [
        {
            "module": "vascular_microcirculation",
            "ligand_symbol": "Vegfa",
            "receptor_symbol": "Kdr/Flt1",
            "display_axis": "Vegfa -> Kdr/Flt1",
            "priority_tier": "core",
        },
        {
            "module": "fibrosis_ecm",
            "ligand_symbol": "Tgfb1",
            "receptor_symbol": "Tgfbr2",
            "display_axis": "Tgfb1 -> Tgfbr2-complex",
            "priority_tier": "core",
        },
        {
            "module": "fibrosis_ecm",
            "ligand_symbol": "Spp1",
            "receptor_symbol": "Cd44",
            "display_axis": "Spp1 -> Cd44",
            "priority_tier": "formal_retained",
        },
    ]


def build_layer_specs() -> list[dict[str, str]]:
    return [
        {
            "layer": "strict",
            "display_label": "Strict formal",
            "min_cells": "10",
            "p_threshold": "0.05",
            "interpretation": "formal_retained",
        },
        {
            "layer": "relaxed",
            "display_label": "Relaxed exploratory",
            "min_cells": "5",
            "p_threshold": "0.20",
            "interpretation": "exploratory_supported",
        },
    ]


def match_axis_rows(
    comm_rows: list[dict[str, str]],
    axis: dict[str, str],
    layer: str,
    condition: str,
) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in comm_rows:
        if row.get("layer") != layer or row.get("condition") != condition:
            continue
        ligand = row.get("ligand", "")
        receptor = row.get("receptor", "")
        if axis["ligand_symbol"] == "Vegfa" and ligand == "Vegfa" and receptor in {"Kdr", "Flt1"}:
            out.append(row)
        elif axis["ligand_symbol"] == "Tgfb1" and ligand == "Tgfb1" and receptor == "Tgfbr2":
            out.append(row)
        elif axis["ligand_symbol"] == "Spp1" and ligand == "Spp1" and receptor == "Cd44":
            out.append(row)
    return out


def build_layer_plot_rows(comm_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for layer_spec in build_layer_specs():
        for axis in build_priority_axis_rows():
            for condition in build_condition_order():
                matched = match_axis_rows(comm_rows, axis, layer_spec["layer"], condition)
                if matched:
                    best = max(matched, key=lambda row: float(row.get("prob", "0") or 0))
                    out_rows.append(
                        {
                            "layer": layer_spec["layer"],
                            "display_label": layer_spec["display_label"],
                            "condition": condition,
                            "display_axis": axis["display_axis"],
                            "prob": f"{float(best.get('prob', '0')):.6f}",
                            "pval": best.get("pval", "NA"),
                            "source_target": f"{best.get('source','NA')} -> {best.get('target','NA')}",
                            "call": layer_spec["interpretation"],
                            "present": "yes",
                        }
                    )
                else:
                    out_rows.append(
                        {
                            "layer": layer_spec["layer"],
                            "display_label": layer_spec["display_label"],
                            "condition": condition,
                            "display_axis": axis["display_axis"],
                            "prob": "0.000000",
                            "pval": "NA",
                            "source_target": "not retained",
                            "call": "not_retained",
                            "present": "no",
                        }
                    )
    return out_rows


def bubble_radius(value: float, max_value: float) -> float:
    if max_value <= 0:
        return 12.0
    return 12.0 + 28.0 * (value / max_value) ** 0.5


def build_dual_layer_svg(rows: list[dict[str, str]], out_path: Path) -> None:
    max_value = max((float(row["prob"]) for row in rows), default=0.0)
    x_map = {"lean": 820, "obese": 1110, "obese_enalapril": 1400}
    layer_y0 = {"strict": 285, "relaxed": 760}
    axes = [row["display_axis"] for row in build_priority_axis_rows()]
    axis_offset = {axis: idx * 130 for idx, axis in enumerate(axes)}
    layer_fill = {"strict": "#ecfeff", "relaxed": "#fff7ed"}
    call_color = {
        "formal_retained": "#0f766e",
        "exploratory_supported": "#c2410c",
        "not_retained": "#94a3b8",
    }

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1820" height="1220" viewBox="0 0 1820 1220">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="60" y="80" font-family="Segoe UI, Arial, sans-serif" font-size="38" font-weight="700" fill="#0f172a">Figure AA. GSE233078 CellChat strict-vs-relaxed layered summary</text>',
        '<text x="60" y="124" font-family="Segoe UI, Arial, sans-serif" font-size="21" fill="#475569">Strict layer uses min.cells = 10 and p &lt;= 0.05; relaxed layer uses min.cells = 5 and p &lt;= 0.20 as exploratory support.</text>',
    ]

    for layer_spec in build_layer_specs():
        y0 = layer_y0[layer_spec["layer"]]
        lines.append(
            f'<rect x="40" y="{y0-95}" width="1740" height="380" rx="24" fill="{layer_fill[layer_spec["layer"]]}" stroke="#cbd5e1" stroke-width="2"/>'
        )
        lines.append(
            f'<text x="60" y="{y0-42}" font-family="Segoe UI, Arial, sans-serif" font-size="28" font-weight="700" fill="#0f172a">{layer_spec["display_label"]}</text>'
        )
        lines.append(
            f'<text x="60" y="{y0-10}" font-family="Segoe UI, Arial, sans-serif" font-size="18" fill="#475569">min.cells = {layer_spec["min_cells"]}; p-threshold = {layer_spec["p_threshold"]}</text>'
        )
        for condition, x in x_map.items():
            lines.append(
                f'<text x="{x}" y="{y0-10}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="22" font-weight="700" fill="#0f172a">{condition}</text>'
            )
        for axis in axes:
            y = y0 + axis_offset[axis]
            lines.append(
                f'<text x="60" y="{y+8}" font-family="Segoe UI, Arial, sans-serif" font-size="22" font-weight="700" fill="#1f2937">{axis}</text>'
            )

    for row in rows:
        x = x_map[row["condition"]]
        y = layer_y0[row["layer"]] + axis_offset[row["display_axis"]]
        prob = float(row["prob"])
        radius = bubble_radius(prob, max_value)
        fill = call_color[row["call"]] if row["present"] == "yes" else "#e2e8f0"
        stroke = call_color[row["call"]]
        dash = "none" if row["present"] == "yes" else "7,5"
        label = f"{prob:.3f}" if row["present"] == "yes" else "NA"
        label_fill = "white" if row["present"] == "yes" else "#475569"
        lines.append(
            f'<circle cx="{x}" cy="{y}" r="{radius:.2f}" fill="{fill}" fill-opacity="0.82" stroke="{stroke}" stroke-width="2.5" stroke-dasharray="{dash}"/>'
        )
        lines.append(
            f'<text x="{x}" y="{y+6}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="700" fill="{label_fill}">{label}</text>'
        )
        lines.append(
            f'<text x="{x}" y="{y+52}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">{row["source_target"]}</text>'
        )

    lines.extend(
        [
            '<text x="60" y="1145" font-family="Segoe UI, Arial, sans-serif" font-size="19" font-weight="700" fill="#0f172a">Reading rule:</text>',
            '<text x="182" y="1145" font-family="Segoe UI, Arial, sans-serif" font-size="17" fill="#475569">Strict layer defines formal retained edges; relaxed layer is exploratory support and should not be written as formal confirmation.</text>',
            '</svg>',
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_summary(
    plot_rows: list[dict[str, str]],
    sensitivity_rows: list[dict[str, str]],
    sample_rows: list[dict[str, str]],
) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def build_r_script(project_root: Path, out_dir: Path) -> str:
    return f"""
options(stringsAsFactors = FALSE)
Sys.setenv(PATH=paste('C:/rtools45/usr/bin','C:/rtools45/x86_64-w64-mingw32.static.posix/bin',Sys.getenv('PATH'),sep=';'))
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

project_root <- "{project_root.as_posix()}"
out_dir <- "{out_dir.as_posix()}"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path(project_root, "results/public_data_gse233078_cellchat_mvp_prep/cellchat_mvp_focus_metadata.csv"), check.names = FALSE)
rownames(meta) <- meta$cell_id

sample_group_counts <- as.data.frame(table(meta$sample_id, meta$condition, meta$cellchat_group), stringsAsFactors = FALSE)
colnames(sample_group_counts) <- c("sample_id", "condition", "cellchat_group", "n_cells")
sample_group_counts$n_cells <- as.integer(sample_group_counts$n_cells)
sample_total_counts <- aggregate(n_cells ~ sample_id + condition, data = sample_group_counts, sum)
sample_total_counts$cellchat_group <- "TOTAL"
sample_counts <- rbind(sample_group_counts[, c("sample_id", "condition", "cellchat_group", "n_cells")], sample_total_counts[, c("sample_id", "condition", "cellchat_group", "n_cells")])
write.csv(sample_counts, file.path(out_dir, "cellchat_focus_sample_cell_counts.csv"), row.names = FALSE)

con <- gzcon(file(file.path(project_root, "data/public_datasets/GSE233078/processed/GSE233078_EXPORT_GEO_counts_postfilter.rds.gz"), "rb"))
raw1 <- readBin(con, what="raw", n=300000000)
close(con)
raw2 <- memDecompress(raw1, type="gzip")
mat <- unserialize(raw2)
focus_cells <- intersect(meta$cell_id, colnames(mat))
meta <- meta[focus_cells, , drop = FALSE]
mat <- mat[, focus_cells, drop = FALSE]

data(CellChatDB.mouse)
custom_db <- CellChatDB.mouse
custom_db$interaction <- subset(
  CellChatDB.mouse$interaction,
  (ligand.symbol == "Vegfa" & grepl("Kdr|Flt1", interaction_name_2, ignore.case = TRUE)) |
  (ligand.symbol == "Tgfb1" & grepl("Tgfbr2", interaction_name_2, ignore.case = TRUE)) |
  (ligand.symbol == "Spp1" & grepl("Cd44", interaction_name_2, ignore.case = TRUE))
)
write.csv(custom_db$interaction, file.path(out_dir, "cellchat_dual_custom_db_interactions.csv"), row.names = FALSE)

future::plan("sequential")
condition_order <- c("lean", "obese", "obese_enalapril")
threshold_grid <- c(0.001, 0.010, 0.050, 0.100, 0.200, 1.000)
layers <- list(
  strict = list(min_cells = 10L, p_threshold = 0.05),
  relaxed = list(min_cells = 5L, p_threshold = 0.20)
)

all_sig_rows <- list()
all_sensitivity_rows <- list()

for (layer_name in names(layers)) {{
  layer_cfg <- layers[[layer_name]]
  for (cond in condition_order) {{
    meta_cond <- meta[meta$condition == cond, , drop = FALSE]
    cells_cond <- rownames(meta_cond)
    expr_cond <- mat[, cells_cond, drop = FALSE]

    seu <- CreateSeuratObject(counts = expr_cond, meta.data = meta_cond)
    seu <- NormalizeData(seu, verbose = FALSE)
    data_input <- GetAssayData(seu, assay = "RNA", layer = "data")
    meta_input <- data.frame(group = meta_cond$cellchat_group, samples = as.factor(meta_cond$sample_id), row.names = rownames(meta_cond))

    cellchat <- createCellChat(object = data_input, meta = meta_input, group.by = "group")
    cellchat@DB <- custom_db
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = layer_cfg$min_cells)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    comm_all <- tryCatch(subsetCommunication(cellchat, thresh = 1), error = function(e) data.frame())
    if (nrow(comm_all) > 0) {{
      comm_all$condition <- cond
      comm_all$layer <- layer_name
      comm_all$min_cells <- layer_cfg$min_cells
      comm_all$p_threshold <- layer_cfg$p_threshold
      write.csv(comm_all, file.path(out_dir, paste0("cellchat_", layer_name, "_", cond, "_all.csv")), row.names = FALSE)

      comm_sig <- subset(comm_all, pval <= layer_cfg$p_threshold)
      if (nrow(comm_sig) > 0) {{
        write.csv(comm_sig, file.path(out_dir, paste0("cellchat_", layer_name, "_", cond, "_significant.csv")), row.names = FALSE)
        all_sig_rows[[paste(layer_name, cond, sep = "_")]] <- comm_sig
      }}

      for (thr in threshold_grid) {{
        rows_thr <- subset(comm_all, pval <= thr)
        all_sensitivity_rows[[length(all_sensitivity_rows) + 1]] <- data.frame(
          layer = layer_name,
          condition = cond,
          min_cells = layer_cfg$min_cells,
          threshold = sprintf("%.3f", thr),
          retained_edges = nrow(rows_thr),
          retained_prob_sum = ifelse(nrow(rows_thr) > 0, sum(rows_thr$prob), 0)
        )
      }}
    }} else {{
      for (thr in threshold_grid) {{
        all_sensitivity_rows[[length(all_sensitivity_rows) + 1]] <- data.frame(
          layer = layer_name,
          condition = cond,
          min_cells = layer_cfg$min_cells,
          threshold = sprintf("%.3f", thr),
          retained_edges = 0,
          retained_prob_sum = 0
        )
      }}
    }}
  }}
}}

if (length(all_sig_rows) > 0) {{
  sig_df <- do.call(rbind, all_sig_rows)
}} else {{
  sig_df <- data.frame(
    source = character(),
    target = character(),
    ligand = character(),
    receptor = character(),
    prob = numeric(),
    pval = numeric(),
    interaction_name = character(),
    interaction_name_2 = character(),
    pathway_name = character(),
    annotation = character(),
    evidence = character(),
    condition = character(),
    layer = character(),
    min_cells = integer(),
    p_threshold = numeric()
  )
}}
write.csv(sig_df, file.path(out_dir, "cellchat_strict_relaxed_significant_communications.csv"), row.names = FALSE)

sensitivity_df <- do.call(rbind, all_sensitivity_rows)
write.csv(sensitivity_df, file.path(out_dir, "cellchat_threshold_sensitivity.csv"), row.names = FALSE)
"""


def run_r_script(project_root: Path, out_dir: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_script(project_root, out_dir))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=7200)
    finally:
        script_path.unlink(missing_ok=True)


def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results" / "public_data_gse233078_cellchat_strict_relaxed"
    out_dir.mkdir(parents=True, exist_ok=True)

    write_csv(out_dir / "cellchat_priority_axes.csv", build_priority_axis_rows(), list(build_priority_axis_rows()[0].keys()))
    write_csv(out_dir / "cellchat_layer_specs.csv", build_layer_specs(), list(build_layer_specs()[0].keys()))
    run_r_script(project_root, out_dir)

    sig_rows = read_csv(out_dir / "cellchat_strict_relaxed_significant_communications.csv")
    sample_rows = read_csv(out_dir / "cellchat_focus_sample_cell_counts.csv")
    sensitivity_rows = read_csv(out_dir / "cellchat_threshold_sensitivity.csv")
    plot_rows = build_layer_plot_rows(sig_rows)
    write_csv(out_dir / "cellchat_strict_relaxed_plot_rows.csv", plot_rows, list(plot_rows[0].keys()))
    build_dual_layer_svg(plot_rows, out_dir / "FigureAA_gse233078_cellchat_strict_relaxed.svg")
    write_text(out_dir / "gse233078_cellchat_strict_relaxed_summary.txt", build_summary(plot_rows, sensitivity_rows, sample_rows))
    write_text(out_dir / "status.txt", "GSE233078 CellChat strict-vs-relaxed completed.\n")

    print(f"GSE233078 CellChat strict-vs-relaxed written to: {out_dir}")


if __name__ == "__main__":
    main()


