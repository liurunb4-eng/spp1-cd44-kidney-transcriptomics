# -*- coding: utf-8 -*-
"""
Step 71: rebuild Figure 5 cellular context with fully marker-excluded PT injury-state coordinates.

Figure contract:
Core conclusion: selected injury-related components align with PT/PT-like injury-state context
when the tested marker is removed both from the direct PCA feature set and from the precomputed
injury score; a target-free PT-state coordinate is also constructed from the residual injury
genes and tubular identity features.
Archetype: schematic-led quantitative grid.
Backend: Python/matplotlib only.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

REPO_ROOT = Path(__file__).resolve().parents[3]
HERE = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
HERE.mkdir(parents=True, exist_ok=True)
ROOT = Path(os.environ.get("AJP_DATA_ROOT", REPO_ROOT))
DEPOSITION_SOURCE = Path(__file__).resolve().parents[2] / "source_data" / "figure5"
GSE131_CELLS = ROOT / "results" / "human_ckd_gse131882_atlas_lite" / "human_ckd_gse131882_atlas_lite_cells.csv"
GSE195_CELLS = ROOT / "results" / "human_ckd_gse195460_first_pass" / "human_ckd_gse195460_celllevel_markers.csv"
PREFIX = HERE / "71_main_figure5_marker_excluded_context"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

COL_TUB = "#4C78A8"
COL_TUB_LIGHT = "#DCE8F5"
COL_IMM = "#B85C5C"
COL_IMM_LIGHT = "#F2D8D6"
COL_ACCENT = "#D55E00"
COL_NEUTRAL = "#6C757D"
COL_LIGHT = "#E9ECEF"
COL_DARK = "#2F3A45"
COL_NEG = "#5B8CC0"
COL_POS = "#C95A5A"

TEST_MARKERS = ["SPP1_logcp10k", "HAVCR1_logcp10k", "LCN2_logcp10k", "VCAM1_logcp10k"]
MARKER_LABEL = {
    "SPP1_logcp10k": "SPP1",
    "HAVCR1_logcp10k": "HAVCR1",
    "LCN2_logcp10k": "LCN2",
    "VCAM1_logcp10k": "VCAM1",
}
MARKER_COLOR = {
    "SPP1": "#3B6EA8",
    "HAVCR1": "#B5533D",
    "LCN2": "#7B6BAE",
    "VCAM1": "#4B8B7A",
}
BASE_FEATURES = ["marker_excluded_injury_score", "SPP1_logcp10k", "HAVCR1_logcp10k", "LCN2_logcp10k", "VCAM1_logcp10k", "LRP2_logcp10k", "SLC34A1_logcp10k"]
TARGET_FREE_FEATURES = ["target_free_injury_score", "PT_score", "LRP2_logcp10k", "SLC34A1_logcp10k"]

# The precomputed injury score is the arithmetic mean of these logCP10K values.
INJURY_SCORE_DEFINITION = {
    "GSE131882": ["SPP1_logcp10k", "HAVCR1_logcp10k", "LCN2_logcp10k", "KRT8_logcp10k", "KRT18_logcp10k"],
    "GSE195460": ["SPP1_logcp10k", "HAVCR1_logcp10k", "LCN2_logcp10k", "VCAM1_logcp10k", "VIM_logcp10k", "KRT8_logcp10k", "KRT18_logcp10k"],
}


def resolve_context_file(*names: str) -> Path | None:
    for base in (HERE, DEPOSITION_SOURCE):
        for name in names:
            candidate = base / name
            if candidate.is_file():
                return candidate
    return None



def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "(no rows)"
    cols = list(df.columns)
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                vals.append(f"{v:.4g}" if np.isfinite(v) else "NA")
            elif pd.isna(v):
                vals.append("NA")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)
def load_pt_cells() -> pd.DataFrame:
    g131_cols = [
        "sample_id", "condition", "PT_score", "Injury_score", "SPP1_logcp10k", "HAVCR1_logcp10k",
        "LCN2_logcp10k", "VCAM1_logcp10k", "LRP2_logcp10k", "SLC34A1_logcp10k", "atlas_lite_celltype", "injury_status",
    ]
    g131 = pd.read_csv(GSE131_CELLS, usecols=g131_cols)
    g131 = g131[g131["atlas_lite_celltype"].eq("PT")].copy()
    g131["dataset"] = "GSE131882"
    g131["compartment"] = "PT"

    g195_cols = [
        "sample_id", "condition", "PT_score", "Injury_score", "SPP1_logcp10k", "HAVCR1_logcp10k",
        "LCN2_logcp10k", "VCAM1_logcp10k", "LRP2_logcp10k", "SLC34A1_logcp10k", "compartment",
    ]
    g195 = pd.read_csv(GSE195_CELLS, usecols=g195_cols)
    g195 = g195[g195["compartment"].eq("PT_like")].copy()
    g195["dataset"] = "GSE195460"
    g195["compartment"] = "PT_like"
    g195["injury_status"] = np.where(g195["Injury_score"] > 0.5, "tubular_injury_like", "tubular_non_injured")

    common = ["dataset", "sample_id", "condition", "compartment", "PT_score", "Injury_score"] + TEST_MARKERS + ["LRP2_logcp10k", "SLC34A1_logcp10k", "injury_status"]
    return pd.concat([g131[common], g195[common]], ignore_index=True)


def rank01(coord: np.ndarray) -> np.ndarray:
    rank = pd.Series(coord).rank(method="average").to_numpy()
    denom = np.nanmax(rank) - np.nanmin(rank)
    if denom <= 0 or not np.isfinite(denom):
        return np.zeros_like(rank, dtype=float)
    return (rank - np.nanmin(rank)) / denom


def pca_coordinate(sub: pd.DataFrame, features: list[str], orient_col: str) -> tuple[np.ndarray, dict]:
    x = sub[features].fillna(0.0).to_numpy(dtype=float)
    x_scaled = StandardScaler().fit_transform(x)
    pca = PCA(n_components=min(2, x_scaled.shape[1]), random_state=20261007)
    pcs = pca.fit_transform(x_scaled)
    coord = pcs[:, 0]
    orient = np.corrcoef(coord, sub[orient_col].fillna(0).to_numpy(dtype=float))[0, 1]
    loading = pca.components_[0].copy()
    if np.isfinite(orient) and orient < 0:
        coord = -coord
        loading = -loading
    meta = {
        "features": ";".join(features),
        "pc1_explained_variance": float(pca.explained_variance_ratio_[0]),
        "oriented_to": orient_col,
    }
    for f, l in zip(features, loading):
        meta[f"loading_{f}"] = float(l)
    return rank01(coord), meta


def spearman(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    d = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 3 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return int(len(d)), np.nan, np.nan
    rho, p = stats.spearmanr(d["x"], d["y"])
    return int(len(d)), float(rho), float(p)


def add_independent_coordinate(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = []
    meta_rows = []
    for dataset, sub in df.groupby("dataset", sort=False):
        sub = sub.copy()
        score_markers = INJURY_SCORE_DEFINITION[dataset]
        excluded_markers = [m for m in TEST_MARKERS if m in score_markers]
        residual_n = len(score_markers) - len(excluded_markers)
        if residual_n <= 0:
            raise ValueError(f"No residual injury genes remain for {dataset}")
        sub["target_free_injury_score"] = (
            len(score_markers) * sub["Injury_score"]
            - sub[excluded_markers].sum(axis=1)
        ) / residual_n
        coord, meta = pca_coordinate(sub, TARGET_FREE_FEATURES, orient_col="target_free_injury_score")
        sub["target_free_pt_state_coordinate"] = coord
        sub["target_free_bin"] = pd.qcut(sub["target_free_pt_state_coordinate"], q=5, labels=["Q1 low", "Q2", "Q3", "Q4", "Q5 high"], duplicates="drop")
        sub["target_free_tertile"] = pd.qcut(sub["target_free_pt_state_coordinate"], q=3, labels=["T1 low", "T2 mid", "T3 high"], duplicates="drop")
        meta["dataset"] = dataset
        meta["coordinate"] = "target_free_pt_state_coordinate"
        meta["original_injury_gene_n"] = len(score_markers)
        meta["excluded_target_markers"] = ";".join(excluded_markers)
        meta["residual_injury_gene_n"] = residual_n
        meta_rows.append(meta)
        out.append(sub)
    return pd.concat(out, ignore_index=True), pd.DataFrame(meta_rows)


def leave_one_marker_out_stats(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    meta_rows = []
    highlow_rows = []
    for dataset, sub0 in df.groupby("dataset", sort=False):
        for marker in TEST_MARKERS:
            sub = sub0.copy()
            score_markers = INJURY_SCORE_DEFINITION[dataset]
            if marker in score_markers:
                sub["marker_excluded_injury_score"] = (
                    len(score_markers) * sub["Injury_score"] - sub[marker]
                ) / (len(score_markers) - 1)
            else:
                sub["marker_excluded_injury_score"] = sub["Injury_score"]
            features = [f for f in BASE_FEATURES if f != marker]
            coord, meta = pca_coordinate(sub, features, orient_col="marker_excluded_injury_score")
            coord_name = f"marker_excluded_coordinate_for_{MARKER_LABEL[marker]}"
            sub[coord_name] = coord
            sub["marker_excluded_tertile"] = pd.qcut(sub[coord_name], q=3, labels=["T1 low", "T2 mid", "T3 high"], duplicates="drop")
            n, rho, p = spearman(sub[coord_name], sub[marker])
            # sample-level high-low differences using the marker-specific leave-one-out coordinate.
            diffs = []
            for (sample, condition), ss in sub.groupby(["sample_id", "condition"], sort=False):
                low = ss[ss["marker_excluded_tertile"].astype(str).eq("T1 low")]
                high = ss[ss["marker_excluded_tertile"].astype(str).eq("T3 high")]
                if low.empty or high.empty:
                    continue
                diff = float(high[marker].mean() - low[marker].mean())
                diffs.append(diff)
                highlow_rows.append({
                    "dataset": dataset,
                    "sample_id": sample,
                    "condition": condition,
                    "marker": marker,
                    "marker_label": MARKER_LABEL[marker],
                    "coordinate_type": "fully_marker_excluded",
                    "coordinate_features": ";".join(features),
                    "low_mean": float(low[marker].mean()),
                    "high_mean": float(high[marker].mean()),
                    "high_minus_low": diff,
                    "low_n_cells": int(len(low)),
                    "high_n_cells": int(len(high)),
                })
            try:
                p_sample = stats.wilcoxon(diffs).pvalue if len(diffs) >= 3 and any(abs(d) > 0 for d in diffs) else np.nan
            except Exception:
                p_sample = np.nan
            rows.append({
                "dataset": dataset,
                "marker": marker,
                "marker_label": MARKER_LABEL[marker],
                "coordinate_type": "fully_marker_excluded",
                "coordinate_features": ";".join(features),
                "marker_removed_from_injury_score": marker in score_markers,
                "cell_n": n,
                "cell_level_spearman_rho": rho,
                "cell_level_spearman_p": p,
                "sample_n": len(diffs),
                "sample_mean_high_minus_low": float(np.nanmean(diffs)) if diffs else np.nan,
                "sample_wilcoxon_p": float(p_sample) if np.isfinite(p_sample) else np.nan,
            })
            meta["dataset"] = dataset
            meta["marker_tested"] = marker
            meta["marker_label"] = MARKER_LABEL[marker]
            meta["coordinate_type"] = "fully_marker_excluded"
            meta["marker_removed_from_injury_score"] = marker in score_markers
            meta_rows.append(meta)
    return pd.DataFrame(rows), pd.DataFrame(meta_rows), pd.DataFrame(highlow_rows)


def independent_bin_trends(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, bin_name), sub in df.groupby(["dataset", "target_free_bin"], observed=True, sort=False):
        row = {
            "dataset": dataset,
            "coordinate_type": "target_free_pt_state_coordinate",
            "coordinate_features": ";".join(TARGET_FREE_FEATURES),
            "trajectory_bin": str(bin_name),
            "n_cells": int(len(sub)),
            "n_samples": int(sub["sample_id"].nunique()),
            "mean_coordinate": float(sub["target_free_pt_state_coordinate"].mean()),
            "mean_target_free_injury_score": float(sub["target_free_injury_score"].mean()),
            "mean_Injury_score": float(sub["Injury_score"].mean()),
            "mean_PT_score": float(sub["PT_score"].mean()),
            "mean_LRP2": float(sub["LRP2_logcp10k"].mean()),
            "mean_SLC34A1": float(sub["SLC34A1_logcp10k"].mean()),
        }
        for marker in TEST_MARKERS:
            label = MARKER_LABEL[marker]
            row[f"mean_{label}"] = float(sub[marker].mean())
            row[f"pct_{label}_positive"] = float((sub[marker] > 0).mean())
        rows.append(row)
    trend = pd.DataFrame(rows)
    order = {"Q1 low": 1, "Q2": 2, "Q3": 3, "Q4": 4, "Q5 high": 5}
    trend["bin_order"] = trend["trajectory_bin"].map(order).astype(int)
    return trend.sort_values(["dataset", "bin_order"])


def make_panelC_long(trend: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in trend.iterrows():
        for marker in TEST_MARKERS:
            gene = MARKER_LABEL[marker]
            rows.append({
                "dataset": row["dataset"],
                "coordinate_type": row["coordinate_type"],
                "coordinate_features": row["coordinate_features"],
                "trajectory_bin": row["trajectory_bin"],
                "bin_order": row["bin_order"],
                "mean_coordinate": row["mean_coordinate"],
                "gene": gene,
                "mean_logcp10k": row[f"mean_{gene}"],
                "pct_positive": row[f"pct_{gene}_positive"],
                "n_cells": row["n_cells"],
                "n_samples": row["n_samples"],
            })
    long = pd.DataFrame(rows)
    long["mean_logcp10k_scaled"] = np.nan
    for (dataset, gene), idx in long.groupby(["dataset", "gene"]).groups.items():
        vals = long.loc[idx, "mean_logcp10k"].astype(float)
        if vals.max() > vals.min():
            long.loc[idx, "mean_logcp10k_scaled"] = (vals - vals.min()) / (vals.max() - vals.min())
        else:
            long.loc[idx, "mean_logcp10k_scaled"] = 0.0
    return long


def panel_label(ax: plt.Axes, label: str, title: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=9, fontweight="bold")
    ax.set_title(title, loc="left", fontsize=8.2, fontweight="bold", pad=10)


def draw_schematic(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    panel_label(ax, "A", "Marker-excluded cellular-context analysis")
    boxes = [
        (0.05, 0.58, 0.28, 0.22, "PT / PT-like cells\nGSE131882 | GSE195460", COL_TUB_LIGHT, COL_TUB),
        (0.39, 0.58, 0.25, 0.22, "Target excluded\nfrom inputs and score", "#F7F1DA", COL_ACCENT),
        (0.71, 0.58, 0.24, 0.22, "Marker association\nand sample contrast", "#E8F0EA", "#4B8B7A"),
        (0.39, 0.18, 0.25, 0.22, "CD44 immune context\nkept as separate evidence", COL_IMM_LIGHT, COL_IMM),
    ]
    for x, y, w, h, txt, fc, ec in boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.025", facecolor=fc, edgecolor=ec, linewidth=1.1)
        ax.add_patch(patch)
        ax.text(x+w/2, y+h/2, txt, ha="center", va="center", fontsize=6.2, color=COL_DARK)
    for x0, y0, x1, y1 in [(0.33,0.69,0.39,0.69),(0.64,0.69,0.71,0.69),(0.515,0.58,0.515,0.40)]:
        ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1), arrowstyle="-|>", mutation_scale=8, lw=0.8, color=COL_NEUTRAL))
    ax.text(0.05, 0.08, "Each displayed association uses a coordinate that excludes the tested marker directly and from the injury score.", fontsize=6.2, color="#555555")


def draw_figure(panelB: pd.DataFrame, panelC: pd.DataFrame, panelD: pd.DataFrame | None) -> None:
    fig = plt.figure(figsize=(8.4, 5.15), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.98, top=0.89, bottom=0.16, wspace=0.42, hspace=0.58)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    draw_schematic(axA)

    # Panel B heatmap: fully marker-excluded rho.
    panel_label(axB, "B", "Cell-level correlations with marker-excluded coordinates")
    heat = panelB.pivot(index="marker_label", columns="dataset", values="cell_level_spearman_rho").reindex(["SPP1", "HAVCR1", "LCN2", "VCAM1"])
    vals = heat.to_numpy(dtype=float)
    vmax = max(0.65, float(np.nanmax(np.abs(vals))))
    im = axB.imshow(vals, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    axB.set_xticks(range(len(heat.columns)), heat.columns, rotation=0)
    axB.set_yticks(range(len(heat.index)), heat.index)
    for i, gene in enumerate(heat.index):
        for j, ds in enumerate(heat.columns):
            v = heat.loc[gene, ds]
            txt = "NA" if pd.isna(v) else f"{v:.2f}"
            axB.text(j, i, txt, ha="center", va="center", fontsize=7, color="white" if pd.notna(v) and abs(v) > vmax*0.45 else "#222222")
    cbar = fig.colorbar(im, ax=axB, fraction=0.045, pad=0.04)
    cbar.set_label("Spearman rho", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)
    axB.tick_params(length=0)
    for s in axB.spines.values():
        s.set_visible(False)

    # Panel C: sample-level high-low contrast from marker-excluded coordinate.
    panel_label(axC, "C", "Sample-level high-minus-low marker contrasts")
    plot = panelB.copy()
    genes = ["SPP1", "HAVCR1", "LCN2", "VCAM1"]
    datasets = ["GSE131882", "GSE195460"]
    x = np.arange(len(genes))
    width = 0.36
    for k, ds in enumerate(datasets):
        sub = plot[plot["dataset"].eq(ds)].set_index("marker_label").reindex(genes)
        vals = sub["sample_mean_high_minus_low"].astype(float).to_numpy()
        offset = (-0.5 + k) * width
        bars = axC.bar(x + offset, vals, width=width, color=COL_TUB if k == 0 else "#7AA374", edgecolor="white", label=ds)
        for xi, v, p in zip(x + offset, vals, sub["sample_wilcoxon_p"].to_numpy()):
            if np.isfinite(p) and p < 0.05:
                axC.text(xi, v + max(vals.max(), 0.1)*0.04, "*", ha="center", va="bottom", fontsize=8, color=COL_ACCENT)
    axC.axhline(0, color="#333333", lw=0.8)
    axC.set_xticks(x, genes)
    axC.set_ylabel("Mean high-minus-low logCP10K")
    axC.set_xlabel("Tested marker")
    axC.grid(axis="y", color=COL_LIGHT, lw=0.6)
    axC.legend(loc="upper right", fontsize=6.2)

    # Panel D: keep CD44 as separate context, lighter emphasis.
    panel_label(axD, "D", "Separate CD44 immune-context summary")
    if panelD is not None and not panelD.empty:
        d = panelD.copy().dropna(subset=["effect_delta"])
        d = d.head(8)
        labels = d["display_short"].tolist()[::-1]
        vals = d["effect_delta"].astype(float).tolist()[::-1]
        y = np.arange(len(vals))
        colors = [COL_POS if v >= 0 else COL_NEG for v in vals]
        axD.barh(y, vals, color=colors, height=0.65, edgecolor="white")
        axD.axvline(0, color="#333333", lw=0.8)
        axD.set_yticks(y, labels)
        axD.set_xlabel("Disease-control delta")
        axD.grid(axis="x", color=COL_LIGHT, lw=0.6)
    else:
        axD.axis("off")
        axD.text(0.02, 0.5, "CD44 immune-context data not available", fontsize=7)

    fig.suptitle("Figure 5. Marker-excluded cellular context for injury-related remodeling components", x=0.01, y=0.985, ha="left", fontsize=10, fontweight="bold")
    fig.text(0.01, 0.020, "Panels B-C exclude each tested marker from both the direct feature set and its injury-score component. CD44 immune context is summarized separately.", ha="left", fontsize=6.1, color="#555555")
    for ext in ["svg", "pdf", "png", "tiff"]:
        if ext in ("png", "tiff"):
            fig.savefig(f"{PREFIX}.{ext}", dpi=600, bbox_inches="tight")
        else:
            fig.savefig(f"{PREFIX}.{ext}", bbox_inches="tight")
    plt.close(fig)

def load_cd44_context() -> pd.DataFrame:
    p = resolve_context_file("38_main_figure5_panelD_cd44_immune_context.csv")
    if p is None:
        return pd.DataFrame()
    d = pd.read_csv(p)
    # Harmonize label and effect column names from existing panelD source.
    if "effect_delta" not in d.columns:
        for cand in ["value", "delta", "mean_delta", "disease_control_delta", "effect"]:
            if cand in d.columns:
                d["effect_delta"] = d[cand]
                break
    if "effect_delta" not in d.columns:
        # Try common source-data names; if unavailable return empty.
        numeric = [c for c in d.columns if c.lower() in ("delta", "effect_delta")]
        if numeric:
            d["effect_delta"] = d[numeric[0]]
        else:
            return pd.DataFrame()
    def short(row):
        raw = str(row.get("display_label", ""))
        if raw in ("", "nan"):
            raw = f"{row.get('source', '')} {row.get('category', '')} {row.get('metric', '')}"
        raw = raw.replace("GSE", "G").replace("mean_", "").replace("nan ", "").strip()
        raw_lower = raw.lower()
        if "g211785" in raw_lower and "mac" in raw_lower and "cd44" in raw_lower:
            if "sc_rna_only" in raw_lower:
                return "G211785 Mac CD44 sc"
            if "all_tech" in raw_lower:
                return "G211785 Mac CD44 all"
        raw = raw.replace("all_tech ", "").replace("sc_rna_only ", "")
        return raw.replace("  ", " ").strip()
    d["display_short"] = d.apply(short, axis=1)
    return d


def main() -> None:
    df = load_pt_cells()
    df_ind, ind_meta = add_independent_coordinate(df)
    loo_stats, loo_meta, loo_highlow = leave_one_marker_out_stats(df_ind)
    trend = independent_bin_trends(df_ind)
    panelC = make_panelC_long(trend)
    panelD = load_cd44_context()

    df_ind.to_csv(HERE / "71_marker_excluded_pt_state_cells.csv.gz", index=False, compression="gzip")
    ind_meta.to_csv(HERE / "71_target_free_coordinate_loadings.csv", index=False, encoding="utf-8-sig")
    loo_stats.to_csv(HERE / "71_figure5_panelB_marker_excluded_correlations.csv", index=False, encoding="utf-8-sig")
    loo_meta.to_csv(HERE / "71_marker_excluded_coordinate_loadings.csv", index=False, encoding="utf-8-sig")
    loo_highlow.to_csv(HERE / "71_marker_excluded_sample_high_low.csv", index=False, encoding="utf-8-sig")
    trend.to_csv(HERE / "71_target_free_pt_state_bin_summary.csv", index=False, encoding="utf-8-sig")
    panelC.to_csv(HERE / "71_figure5_panelC_target_free_coordinate_trends.csv", index=False, encoding="utf-8-sig")

    draw_figure(loo_stats, panelC, panelD)

    old_path = resolve_context_file(
        "38_main_figure5_panelB_injury_coordinate_correlations.csv",
        "71_old_vs_marker_excluded_correlations.csv",
    )
    if old_path is None:
        old = pd.DataFrame(columns=["dataset", "marker", "old_circular_rho"])
    else:
        old_stats = pd.read_csv(old_path)
        if "old_circular_rho" in old_stats.columns:
            old = old_stats[["dataset", "marker", "old_circular_rho"]].drop_duplicates()
        else:
            old = old_stats[old_stats["marker"].isin(TEST_MARKERS)][["dataset", "marker", "cell_level_spearman_rho"]].rename(columns={"cell_level_spearman_rho": "old_circular_rho"})
    comp = loo_stats.merge(old, on=["dataset", "marker"], how="left")
    comp.to_csv(HERE / "71_old_vs_marker_excluded_correlations.csv", index=False, encoding="utf-8-sig")

    md = []
    md.append("# Figure 5 marker-excluded rebuild, 2026-10-09\n\n")
    md.append("## Core conclusion\n\n")
    md.append("Selected injury-related components align with PT/PT-like cellular injury-state context after the tested marker is removed both directly and from the precomputed injury score.\n\n")
    md.append("## What changed\n\n")
    md.append("- Panels B-C use marker-excluded coordinates: the tested marker is removed from both the direct PCA features and its contribution to the injury score.\n")
    md.append("- A target-free PT-state coordinate is built from the residual injury genes, PT score, LRP2 and SLC34A1, excluding SPP1, HAVCR1, LCN2 and VCAM1.\n")
    md.append("- CD44 is summarized as separate immune-context evidence and is not used to define the tubular coordinate.\n\n")
    md.append("## Fully marker-excluded correlations\n\n")
    md.append(md_table(loo_stats[["dataset","marker_label","cell_level_spearman_rho","cell_level_spearman_p","sample_n","sample_mean_high_minus_low","sample_wilcoxon_p"]]))
    md.append("\n\n## Old versus fully marker-excluded comparison\n\n")
    md.append(md_table(comp[["dataset","marker_label","old_circular_rho","cell_level_spearman_rho"]]))
    md.append("\n\n## Output files\n\n")
    for fn in [
        "71_main_figure5_marker_excluded_context.svg",
        "71_main_figure5_marker_excluded_context.pdf",
        "71_main_figure5_marker_excluded_context.png",
        "71_main_figure5_marker_excluded_context.tiff",
        "71_figure5_panelB_marker_excluded_correlations.csv",
        "71_figure5_panelC_target_free_coordinate_trends.csv",
        "71_marker_excluded_sample_high_low.csv",
        "71_target_free_pt_state_bin_summary.csv",
        "71_target_free_coordinate_loadings.csv",
        "71_marker_excluded_coordinate_loadings.csv",
        "71_old_vs_marker_excluded_correlations.csv",
    ]:
        md.append(f"- `{fn}`\n")
    (HERE / "71_figure5_marker_excluded_rebuild_report.md").write_text("".join(md), encoding="utf-8")

    print("DONE")
    print(PREFIX.with_suffix(".png"))
    print(HERE / "71_figure5_marker_excluded_rebuild_report.md")
    print(loo_stats[["dataset","marker_label","cell_level_spearman_rho","sample_mean_high_minus_low"]].to_string(index=False))


if __name__ == "__main__":
    main()



# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
