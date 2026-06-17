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
            "receptor_pattern": "Kdr|Flt1",
            "display_axis": "Vegfa -> Kdr/Flt1",
            "priority_tier": "core",
        },
        {
            "module": "fibrosis_ecm",
            "ligand_symbol": "Tgfb1",
            "receptor_pattern": "Tgfbr2",
            "display_axis": "Tgfb1 -> Tgfbr2-complex",
            "priority_tier": "core",
        },
        {
            "module": "fibrosis_ecm",
            "ligand_symbol": "Spp1",
            "receptor_pattern": "Cd44",
            "display_axis": "Spp1 -> Cd44",
            "priority_tier": "hypothesis",
        },
    ]


def match_axis_rows(comm_rows: list[dict[str, str]], axis: dict[str, str], condition: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in comm_rows:
        if row.get("condition") != condition:
            continue
        ligand = row.get("ligand", "")
        receptor = row.get("receptor", "")
        if axis["ligand_symbol"] == "Vegfa" and ligand == "Vegfa" and receptor in {"Kdr", "Flt1"}:
            out.append(row)
        elif axis["ligand_symbol"] == "Tgfb1" and ligand == "Tgfb1":
            out.append(row)
        elif axis["ligand_symbol"] == "Spp1" and ligand == "Spp1" and receptor == "Cd44":
            out.append(row)
    return out


def build_formal_plot_rows(comm_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for axis in build_priority_axis_rows():
        for condition in build_condition_order():
            matched = match_axis_rows(comm_rows, axis, condition)
            if matched:
                best = max(matched, key=lambda row: float(row.get("prob", "0") or 0))
                out_rows.append(
                    {
                        "display_axis": axis["display_axis"],
                        "condition": condition,
                        "prob": f"{float(best.get('prob', '0')):.6f}",
                        "pval": best.get("pval", "NA"),
                        "threshold_mode": best.get("threshold_mode", "significant_only"),
                        "source_target": f"{best.get('source','NA')} -> {best.get('target','NA')}",
                        "present": "yes",
                    }
                )
            else:
                out_rows.append(
                    {
                        "display_axis": axis["display_axis"],
                        "condition": condition,
                        "prob": "0.000000",
                        "pval": "NA",
                        "threshold_mode": "not_retained",
                        "source_target": "not retained",
                        "present": "no",
                    }
                )
    return out_rows


def bubble_radius(value: float, max_value: float) -> float:
    if max_value <= 0:
        return 10.0
    return 10.0 + 30.0 * (value / max_value) ** 0.5


def build_formal_bubble_svg(rows: list[dict[str, str]], out_path: Path) -> None:
    conditions = build_condition_order()
    axes = [row["display_axis"] for row in build_priority_axis_rows()]
    max_value = max(float(row["prob"]) for row in rows) if rows else 0.0
    x_map = {"lean": 760, "obese": 1040, "obese_enalapril": 1320}
    y_map = {
        axis: 320 + idx * 190
        for idx, axis in enumerate(axes)
    }

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1700" height="980" viewBox="0 0 1700 980">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="60" y="82" font-family="Segoe UI, Arial, sans-serif" font-size="38" font-weight="700" fill="#1f2937">Figure X. GSE233078 formal CellChat MVP bubble summary</text>',
        '<text x="60" y="126" font-family="Segoe UI, Arial, sans-serif" font-size="22" font-weight="400" fill="#4b5563">Bubble sizes come from formal CellChat communication probability; blank cells mean the axis was not retained at the current threshold.</text>',
    ]
    for condition in conditions:
        lines.append(f'<text x="{x_map[condition]}" y="220" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">{condition}</text>')
    for axis in axes:
        lines.append(f'<text x="60" y="{y_map[axis]+8}" font-family="Segoe UI, Arial, sans-serif" font-size="22" font-weight="700" fill="#1f2937">{axis}</text>')
    for row in rows:
        x = x_map[row["condition"]]
        y = y_map[row["display_axis"]]
        prob = float(row["prob"])
        radius = bubble_radius(prob, max_value)
        present = row["present"] == "yes"
        fill = "#0f766e" if present else "#e2e8f0"
        stroke = "#0f766e" if row["threshold_mode"] == "significant_only" else "#94a3b8"
        dash = "none" if row["threshold_mode"] == "significant_only" else "6,4"
        lines.append(
            f'<circle cx="{x}" cy="{y}" r="{radius:.2f}" fill="{fill}" fill-opacity="0.78" stroke="{stroke}" stroke-width="2" stroke-dasharray="{dash}"/>'
        )
        label = f'{prob:.3f}' if present else 'NA'
        label_fill = "white" if present else "#475569"
        lines.append(f'<text x="{x}" y="{y+6}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="700" fill="{label_fill}">{label}</text>')
        lines.append(f'<text x="{x}" y="{y+52}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" font-weight="400" fill="#475569">{row["source_target"]}</text>')
    lines.append('<text x="60" y="930" font-family="Segoe UI, Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2937">Note:</text>')
    lines.append('<text x="130" y="930" font-family="Segoe UI, Arial, sans-serif" font-size="17" font-weight="400" fill="#475569">Solid outline = retained as significant by formal CellChat. Dashed outline = not retained or only exploratory fallback.</text>')
    lines.append('</svg>')
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_r_script(project_root: Path, out_dir: Path) -> str:
    return f"""
options(stringsAsFactors = FALSE)
Sys.setenv(PATH=paste('C:/rtools45/usr/bin','C:/rtools45/x86_64-w64-mingw32.static.posix/bin',Sys.getenv('PATH'),sep=';'))
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(svglite))

project_root <- "{project_root.as_posix()}"
out_dir <- "{out_dir.as_posix()}"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path(project_root, "results/public_data_gse233078_cellchat_mvp_prep/cellchat_mvp_focus_metadata.csv"), check.names = FALSE)
rownames(meta) <- meta$cell_id

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
write.csv(custom_db$interaction, file.path(out_dir, "cellchat_formal_custom_db_interactions.csv"), row.names = FALSE)

future::plan("sequential")
condition_order <- c("lean", "obese", "obese_enalapril")
cellchat_list <- list()
comm_rows_all <- list()

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
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  comm <- tryCatch(
    {{
      tmp <- subsetCommunication(cellchat, thresh = 0.05)
      tmp$threshold_mode <- "significant_only"
      tmp
    }},
    error = function(e) {{
      tmp <- tryCatch(subsetCommunication(cellchat, thresh = 1), error = function(e2) data.frame())
      if (nrow(tmp) > 0) {{
        tmp$threshold_mode <- "all_inferred"
      }}
      tmp
    }}
  )
  if (nrow(comm) > 0) {{
    comm$condition <- cond
    write.csv(comm, file.path(out_dir, paste0("cellchat_", cond, "_communications.csv")), row.names = FALSE)
  }}

  saveRDS(cellchat, file.path(out_dir, paste0("cellchat_", cond, ".rds")))
  cellchat_list[[cond]] <- cellchat
  comm_rows_all[[cond]] <- comm
}}

merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
saveRDS(merged, file.path(out_dir, "cellchat_3cond_object.rds"))

non_empty <- Filter(function(x) !is.null(x) && nrow(x) > 0, comm_rows_all)
if (length(non_empty) > 0) {{
  all_comm <- do.call(rbind, non_empty)
}} else {{
  all_comm <- data.frame(
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
    threshold_mode = character(),
    condition = character()
  )
}}
write.csv(all_comm, file.path(out_dir, "cellchat_formal_combined_communications.csv"), row.names = FALSE)
"""


def run_r_script(project_root: Path, out_dir: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(build_r_script(project_root, out_dir))
        script_path = Path(handle.name)
    try:
        subprocess.run([str(RSCRIPT), str(script_path)], check=True, cwd=project_root, timeout=5400)
    finally:
        script_path.unlink(missing_ok=True)


def build_summary(comm_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "results" / "public_data_gse233078_cellchat_formal_mvp"
    out_dir.mkdir(parents=True, exist_ok=True)

    write_csv(out_dir / "cellchat_formal_priority_axes.csv", build_priority_axis_rows(), list(build_priority_axis_rows()[0].keys()))
    run_r_script(project_root, out_dir)

    combined_csv = out_dir / "cellchat_formal_combined_communications.csv"
    comm_rows = read_csv(combined_csv) if combined_csv.exists() else []
    formal_plot_rows = build_formal_plot_rows(comm_rows)
    write_csv(out_dir / "cellchat_formal_plot_rows.csv", formal_plot_rows, list(formal_plot_rows[0].keys()))
    build_formal_bubble_svg(formal_plot_rows, out_dir / "FigureX_gse233078_cellchat_formal_mvp_bubble.svg")
    write_text(out_dir / "gse233078_cellchat_formal_mvp_summary.txt", build_summary(comm_rows))
    write_text(out_dir / "status.txt", "GSE233078 formal CellChat MVP completed.\n")

    print(f"GSE233078 formal CellChat MVP written to: {out_dir}")


if __name__ == "__main__":
    main()


