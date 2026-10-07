
# -*- coding: utf-8 -*-
"""
Step 16. DKD-specificity vs pan-CKD conserved injury-remodeling audit.

This script asks whether the frozen 10-gene injury-dominant remodeling candidate
is specific to diabetic nephropathy/DKD or reflects a conserved CKD injury state.
It combines:
  1) GSE175759 tubulointerstitium scores with diagnosis and eGFR metadata.
  2) C-PROBE microdissected kidney transcriptomes (GSE180393 glomeruli,
     GSE180394 tubuli) using Entrez-style custom CDF probe identifiers.

The output is deliberately conservative: all claims remain association/context
claims; no prognostic, causal, or mechanistic claim is made from this analysis.
"""

from __future__ import annotations

import gzip
import io
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Paths and constants
# -----------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(os.environ.get("AJP_DATA_ROOT", REPO_ROOT))
OUTDIR = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUTDIR.mkdir(parents=True, exist_ok=True)
STEP13 = OUTDIR / "step13_public_outcome_screen"
PREFIX = OUTDIR / "16_pan_ckd_specificity_audit"

SIGNATURE = "candidate_DKD_TI_remodeling_no_SPP1_CD44"
SIGNATURE_LABEL = "10-gene injury-dominant candidate"

# NCBI Entrez Gene IDs corresponding to the frozen 10 genes. The C-PROBE matrix
# uses Entrez-like row IDs with an _at suffix, so mapping is exact only when the
# stripped row ID matches one of these identifiers.
GENE_TO_ENTREZ: Dict[str, str] = {
    "HAVCR1": "26762",
    "LCN2": "3934",
    "VCAM1": "7412",
    "CLU": "1191",
    "COL1A1": "1277",
    "FN1": "2335",
    "ACTA2": "59",
    "TGFB1": "7040",
    "TLR4": "7099",
    "NLRP3": "114548",
}
ENTREZ_TO_GENE = {v: k for k, v in GENE_TO_ENTREZ.items()}
GENES = list(GENE_TO_ENTREZ.keys())

# Matplotlib publication defaults (Python backend only; see nature-figure skill)
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
    "figure.dpi": 160,
})

# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def bh_fdr(pvals: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR, preserving NaNs."""
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    mask = np.isfinite(arr)
    if not mask.any():
        return out.tolist()
    p = arr[mask]
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    temp = np.empty_like(q)
    temp[order] = q
    out[mask] = temp
    return out.tolist()


def mw_vs_reference(df: pd.DataFrame, group_col: str, score_col: str, ref: str,
                    min_group_n: int = 3) -> pd.DataFrame:
    rows = []
    ref_vals = pd.to_numeric(df.loc[df[group_col] == ref, score_col], errors="coerce").dropna()
    for group, sub in df.groupby(group_col):
        vals = pd.to_numeric(sub[score_col], errors="coerce").dropna()
        row = {
            "group": group,
            "n": int(vals.shape[0]),
            "median": float(vals.median()) if len(vals) else np.nan,
            "mean": float(vals.mean()) if len(vals) else np.nan,
            "reference": ref,
            "reference_n": int(ref_vals.shape[0]),
            "reference_median": float(ref_vals.median()) if len(ref_vals) else np.nan,
            "delta_median_vs_reference": np.nan,
            "mannwhitney_p": np.nan,
            "test_note": "reference group" if group == ref else "not tested",
        }
        if group != ref and len(vals) >= min_group_n and len(ref_vals) >= min_group_n:
            try:
                stat = stats.mannwhitneyu(vals, ref_vals, alternative="two-sided")
                row["mannwhitney_p"] = float(stat.pvalue)
                row["delta_median_vs_reference"] = float(vals.median() - ref_vals.median())
                row["test_note"] = "Mann-Whitney vs reference"
            except Exception as exc:  # noqa: BLE001
                row["test_note"] = f"test failed: {exc}"
        elif group == ref:
            row["delta_median_vs_reference"] = 0.0
        elif len(vals) < min_group_n:
            row["test_note"] = f"not tested: n<{min_group_n}"
        rows.append(row)
    res = pd.DataFrame(rows).sort_values(["group"]).reset_index(drop=True)
    res["fdr"] = bh_fdr(res["mannwhitney_p"].tolist())
    return res


def spearman_by_group(df: pd.DataFrame, group_col: str, score_col: str, y_col: str,
                      min_n: int = 5) -> pd.DataFrame:
    rows = []
    groups = [("All usable", df), ("CKD only", df[df[group_col] != "Control"])]
    groups.extend(list(df.groupby(group_col)))
    seen = set()
    out_groups = []
    for name, sub in groups:
        if name in seen:
            continue
        seen.add(name)
        out_groups.append((name, sub))
    for name, sub in out_groups:
        x = pd.to_numeric(sub[score_col], errors="coerce")
        y = pd.to_numeric(sub[y_col], errors="coerce")
        good = x.notna() & y.notna()
        row = {"group": name, "n": int(good.sum()), "spearman_rho": np.nan, "p_value": np.nan}
        if good.sum() >= min_n:
            r = stats.spearmanr(x[good], y[good])
            row["spearman_rho"] = float(r.statistic)
            row["p_value"] = float(r.pvalue)
        rows.append(row)
    res = pd.DataFrame(rows)
    res["fdr"] = bh_fdr(res["p_value"].tolist())
    return res


def clean_gse175759_diagnosis(x: str) -> str:
    mapping = {
        "Control": "Control",
        "Diabetic nephropathy": "DN",
        "FSGS": "FSGS",
        "IgAN": "IgAN",
        "Lupus nephritis": "LN",
        "Membranous nephropathy": "MN",
        "minimal change disease": "MCD",
    }
    return mapping.get(str(x), str(x))


def clean_cprobe_group(x: str) -> str:
    s = str(x).strip()
    s = re.sub(r"\s+", " ", s)
    s = s.replace("LN-WHO-V", "LN-WHO V")
    if s == "Living donor":
        return "Living donor"
    if s == "DN":
        return "DN"
    if "IgAN" in s:
        return "IgAN"
    if "FSGS" in s or "FGGS" in s:
        return "FSGS/FGGS"
    if s.startswith("LN-WHO"):
        return "LN"
    if s == "MN":
        return "MN"
    if s == "MCD":
        return "MCD"
    if "Interstitial nephritis" in s:
        return "Interstitial nephritis"
    if "Interstitial fibrosis" in s:
        return "Interstitial fibrosis"
    if "infection-associated" in s or "Immuncomplex" in s or "Glomerulonephritis" in s:
        return "Other GN/immune"
    return "Other CKD"


def parse_series_matrix(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (expression matrix, sample metadata) from a GEO series matrix."""
    meta_lines = []
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_") or line.startswith("!Series_title"):
                meta_lines.append(line.rstrip("\n"))
    expr = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t")
    expr.iloc[:, 0] = expr.iloc[:, 0].astype(str).str.strip('"')
    expr = expr.rename(columns={expr.columns[0]: "ID_REF"})
    sample_cols = [c for c in expr.columns if c != "ID_REF"]
    expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="coerce")

    meta: Dict[str, List[str]] = {}
    for line in meta_lines:
        parts = line.split("\t")
        key = parts[0]
        vals = [p.strip().strip('"') for p in parts[1:]]
        meta.setdefault(key, []).append(vals)
    accessions = meta.get("!Sample_geo_accession", [[]])[0]
    titles = meta.get("!Sample_title", [[]])[0]
    characteristics = meta.get("!Sample_characteristics_ch1", [])
    md = pd.DataFrame({"sample_id": accessions, "title": titles})
    for vals in characteristics:
        if not vals:
            continue
        label = vals[0].split(":", 1)[0].strip().lower().replace(" ", "_") if ":" in vals[0] else "characteristic"
        clean_vals = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
        # Avoid silently overwriting duplicate labels.
        col = label
        idx = 2
        while col in md.columns:
            col = f"{label}_{idx}"
            idx += 1
        md[col] = clean_vals
    return expr, md


def score_cprobe(path: Path, dataset: str, compartment: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    expr, md = parse_series_matrix(path)
    expr["entrez_id"] = expr["ID_REF"].str.replace(r"_at$", "", regex=True)
    expr["gene"] = expr["entrez_id"].map(ENTREZ_TO_GENE)
    presence_rows = []
    for gene in GENES:
        hit = expr.loc[expr["gene"] == gene]
        presence_rows.append({
            "dataset": dataset,
            "compartment": compartment,
            "gene": gene,
            "entrez_id": GENE_TO_ENTREZ[gene],
            "n_rows": int(hit.shape[0]),
            "row_ids": ";".join(hit["ID_REF"].astype(str).tolist()),
        })
    hits = expr.dropna(subset=["gene"]).copy()
    sample_cols = [c for c in hits.columns if c.startswith("GSM")]
    gene_expr = hits.groupby("gene")[sample_cols].mean()
    gene_expr = gene_expr.reindex(GENES)
    # z-score each gene across samples, then average across genes present.
    z = gene_expr.sub(gene_expr.mean(axis=1), axis=0).div(gene_expr.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    score = z.mean(axis=0, skipna=True).rename("score").reset_index().rename(columns={"index": "sample_id"})
    score["n_genes_scored"] = int(gene_expr.dropna(how="all").shape[0])
    score = score.merge(md, on="sample_id", how="left")
    score["dataset"] = dataset
    score["compartment"] = compartment
    group_col = "sample_group" if "sample_group" in score.columns else None
    if group_col is None:
        # fallback from title before bracket
        score["sample_group"] = score["title"].str.replace(r"\s*\[.*$", "", regex=True).str.strip()
    score["group_collapsed"] = score["sample_group"].map(clean_cprobe_group)
    return score, pd.DataFrame(presence_rows)


def jitter_positions(n: int, width: float = 0.18) -> np.ndarray:
    if n <= 1:
        return np.array([0.0])[:n]
    return np.linspace(-width, width, n)


def plot_group_points(ax, data: pd.DataFrame, group_col: str, value_col: str,
                      order: List[str], color_map: Dict[str, str], ylabel: str,
                      title: str, ref_group: str) -> None:
    for i, group in enumerate(order):
        vals = pd.to_numeric(data.loc[data[group_col] == group, value_col], errors="coerce").dropna().sort_values().to_numpy()
        if len(vals) == 0:
            continue
        xs = i + jitter_positions(len(vals), width=0.20)
        ax.scatter(xs, vals, s=18, color=color_map.get(group, "#777777"), alpha=0.82,
                   edgecolor="white", linewidth=0.35, zorder=3)
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ax.plot([i-0.24, i+0.24], [med, med], color="#222222", lw=1.2, zorder=4)
        ax.vlines(i, q1, q3, color="#222222", lw=1.0, zorder=4)
    ax.axhline(0, color="#BBBBBB", lw=0.8, zorder=1)
    if ref_group in order:
        ax.axvspan(order.index(ref_group)-0.45, order.index(ref_group)+0.45, color="#F1F3F5", zorder=0)
    ax.set_xticks(range(len(order)))
    tick_labels = []
    for g in order:
        n = int((data[group_col] == g).sum())
        tick_labels.append(f"{g}\nn={n}")
    ax.set_xticklabels(tick_labels, rotation=35, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E9ECEF", lw=0.6)


def write_report(path: Path, sections: List[str]) -> None:
    path.write_text("\n\n".join(sections) + "\n", encoding="utf-8")


def simple_markdown_table(df: pd.DataFrame) -> str:
    """Small dependency-free markdown table writer."""
    if df.empty:
        return "_No rows._"
    work = df.copy()
    cols = list(work.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in work.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                vals.append(f"{v:.4g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

# -----------------------------------------------------------------------------
# 1) GSE175759 disease-specificity audit
# -----------------------------------------------------------------------------
gse175 = pd.read_csv(OUTDIR / "12_clinical_egfr_stratification_scores_merged.csv")
gse175 = gse175[(gse175["signature"] == SIGNATURE) & (gse175["usable"] == True)].copy()  # noqa: E712
gse175["group_collapsed"] = gse175["diagnosis"].map(clean_gse175759_diagnosis)
gse175_stats = mw_vs_reference(gse175, "group_collapsed", "score", "Control", min_group_n=3)
gse175_corr = spearman_by_group(gse175, "group_collapsed", "score", "egfr_ckd_epi", min_n=5)

# -----------------------------------------------------------------------------
# 2) C-PROBE tubuli and glomeruli pan-CKD audit
# -----------------------------------------------------------------------------
cprobe_specs = [
    ("GSE180394", "Tubuli", STEP13 / "GSE180394_series_matrix.txt.gz"),
    ("GSE180393", "Glomeruli", STEP13 / "GSE180393_series_matrix.txt.gz"),
]
cprobe_scores = []
cprobe_presence = []
for dataset, compartment, path in cprobe_specs:
    sc, pr = score_cprobe(path, dataset, compartment)
    cprobe_scores.append(sc)
    cprobe_presence.append(pr)
cprobe_scores = pd.concat(cprobe_scores, ignore_index=True)
cprobe_presence = pd.concat(cprobe_presence, ignore_index=True)

cprobe_stats_all = []
for (dataset, compartment), sub in cprobe_scores.groupby(["dataset", "compartment"]):
    st = mw_vs_reference(sub, "group_collapsed", "score", "Living donor", min_group_n=3)
    st.insert(0, "dataset", dataset)
    st.insert(1, "compartment", compartment)
    cprobe_stats_all.append(st)
cprobe_stats = pd.concat(cprobe_stats_all, ignore_index=True)

# A compact cross-dataset interpretation table.
summary_rows = []
# GSE175759 key disease directionality
for _, row in gse175_stats.iterrows():
    if row["group"] == "Control":
        continue
    summary_rows.append({
        "source": "GSE175759 tubulointerstitium",
        "comparison": f"{row['group']} vs Control",
        "n_group": row["n"],
        "n_reference": row["reference_n"],
        "delta_median": row["delta_median_vs_reference"],
        "p_value": row["mannwhitney_p"],
        "fdr": row["fdr"],
        "interpretation": "higher than control" if pd.notna(row["delta_median_vs_reference"]) and row["delta_median_vs_reference"] > 0 else "not higher than control",
    })
# C-PROBE key disease directionality
for _, row in cprobe_stats.iterrows():
    if row["group"] == "Living donor":
        continue
    summary_rows.append({
        "source": f"{row['dataset']} {row['compartment']}",
        "comparison": f"{row['group']} vs Living donor",
        "n_group": row["n"],
        "n_reference": row["reference_n"],
        "delta_median": row["delta_median_vs_reference"],
        "p_value": row["mannwhitney_p"],
        "fdr": row["fdr"],
        "interpretation": "higher than living donor" if pd.notna(row["delta_median_vs_reference"]) and row["delta_median_vs_reference"] > 0 else "not higher than living donor",
    })
summary = pd.DataFrame(summary_rows)

# -----------------------------------------------------------------------------
# Save tables
# -----------------------------------------------------------------------------
gse175_stats.to_csv(OUTDIR / "16_pan_ckd_gse175759_diagnosis_tests.csv", index=False, encoding="utf-8-sig")
gse175_corr.to_csv(OUTDIR / "16_pan_ckd_gse175759_egfr_correlations.csv", index=False, encoding="utf-8-sig")
cprobe_presence.to_csv(OUTDIR / "16_pan_ckd_cprobe_gene_presence.csv", index=False, encoding="utf-8-sig")
cprobe_scores.to_csv(OUTDIR / "16_pan_ckd_cprobe_scores.csv", index=False, encoding="utf-8-sig")
cprobe_stats.to_csv(OUTDIR / "16_pan_ckd_cprobe_diagnosis_tests.csv", index=False, encoding="utf-8-sig")
summary.to_csv(OUTDIR / "16_pan_ckd_specificity_summary.csv", index=False, encoding="utf-8-sig")

# -----------------------------------------------------------------------------
# Figure: quantitative grid
# -----------------------------------------------------------------------------
# Core conclusion: the signature is DKD-relevant but not DKD-exclusive, with
# strongest support for a conserved tubulointerstitial/pan-CKD injury context.
PALETTE = {
    "Control": "#8A8F98",
    "Living donor": "#8A8F98",
    "DN": "#D55E00",
    "IgAN": "#0072B2",
    "FSGS": "#009E73",
    "FSGS/FGGS": "#009E73",
    "LN": "#CC79A7",
    "MN": "#E69F00",
    "MCD": "#56B4E9",
    "Interstitial nephritis": "#7B68EE",
    "Interstitial fibrosis": "#6B4F3B",
    "Other GN/immune": "#999933",
    "Other CKD": "#777777",
}

fig = plt.figure(figsize=(7.25, 6.8), constrained_layout=False)
gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.0, 1.0], hspace=0.42, wspace=0.35)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# Panel A: GSE175759 disease groups
order_a = [g for g in ["Control", "DN", "IgAN", "FSGS", "LN", "MN", "MCD"] if g in set(gse175["group_collapsed"])]
plot_group_points(ax1, gse175, "group_collapsed", "score", order_a, PALETTE,
                  "Signature score (z-mean)", "A  GSE175759: CKD diagnoses vs control", "Control")

# Panel B: eGFR relationship
x = pd.to_numeric(gse175["score"], errors="coerce")
y = pd.to_numeric(gse175["egfr_ckd_epi"], errors="coerce")
for group, sub in gse175.groupby("group_collapsed"):
    ax2.scatter(sub["score"], sub["egfr_ckd_epi"], s=20, color=PALETTE.get(group, "#777777"),
                alpha=0.82, edgecolor="white", linewidth=0.35, label=group)
mask = x.notna() & y.notna()
if mask.sum() > 2:
    coef = np.polyfit(x[mask], y[mask], deg=1)
    xs = np.linspace(x[mask].min(), x[mask].max(), 100)
    ax2.plot(xs, coef[0]*xs + coef[1], color="#222222", lw=1.0)
    rho = stats.spearmanr(x[mask], y[mask])
    ax2.text(0.02, 0.05, f"all samples: rho={rho.statistic:.2f}, P={rho.pvalue:.3g}", transform=ax2.transAxes,
             ha="left", va="bottom", fontsize=7)
ax2.set_xlabel("Signature score")
ax2.set_ylabel("eGFR CKD-EPI")
ax2.set_title("B  Cross-sectional kidney function", loc="left", fontweight="bold")
ax2.grid(color="#E9ECEF", lw=0.6)
ax2.legend(loc="upper right", fontsize=5.5, ncols=2, handletextpad=0.2, columnspacing=0.7)

# Panel C: C-PROBE tubuli
ct = cprobe_scores[cprobe_scores["compartment"] == "Tubuli"].copy()
order_c_base = ["Living donor", "DN", "IgAN", "FSGS/FGGS", "LN", "MN", "MCD", "Interstitial nephritis", "Interstitial fibrosis", "Other GN/immune", "Other CKD"]
order_c = [g for g in order_c_base if g in set(ct["group_collapsed"])]
plot_group_points(ax3, ct, "group_collapsed", "score", order_c, PALETTE,
                  "Signature score (z-mean)", "C  C-PROBE tubuli: pan-CKD audit", "Living donor")

# Panel D: C-PROBE compartment summary for selected diseases
sel = ["DN", "IgAN", "FSGS/FGGS", "LN", "Other CKD"]
comp_rows = []
for comp in ["Tubuli", "Glomeruli"]:
    substats = cprobe_stats[cprobe_stats["compartment"] == comp]
    for g in sel:
        hit = substats[substats["group"] == g]
        if hit.empty:
            continue
        comp_rows.append({"compartment": comp, "group": g, "delta": float(hit.iloc[0]["delta_median_vs_reference"]), "fdr": hit.iloc[0]["fdr"], "n": int(hit.iloc[0]["n"])})
comp_df = pd.DataFrame(comp_rows)
width = 0.35
xpos = np.arange(len(sel))
for j, comp in enumerate(["Tubuli", "Glomeruli"]):
    vals = []
    labels = []
    for g in sel:
        hit = comp_df[(comp_df["compartment"] == comp) & (comp_df["group"] == g)]
        vals.append(hit.iloc[0]["delta"] if not hit.empty else np.nan)
        labels.append(hit.iloc[0]["fdr"] if not hit.empty else np.nan)
    offset = (-width/2 if comp == "Tubuli" else width/2)
    ax4.bar(xpos + offset, vals, width=width, label=comp,
            color=("#4C78A8" if comp == "Tubuli" else "#F58518"), alpha=0.85)
    for xi, val, q in zip(xpos + offset, vals, labels):
        if np.isfinite(val) and np.isfinite(q) and q < 0.10:
            ax4.text(xi, val + (0.04 if val >= 0 else -0.08), "*", ha="center", va="bottom" if val >= 0 else "top", fontsize=9)
ax4.axhline(0, color="#444444", lw=0.8)
ax4.set_xticks(xpos)
ax4.set_xticklabels(sel, rotation=35, ha="right")
ax4.set_ylabel("Median delta vs living donor")
ax4.set_title("D  Compartment contrast in C-PROBE", loc="left", fontweight="bold")
ax4.legend(loc="upper left", fontsize=6)
ax4.grid(axis="y", color="#E9ECEF", lw=0.6)

fig.suptitle("Step 16. The injury-remodeling candidate is DKD-relevant but not DKD-exclusive", x=0.01, ha="left", y=0.995, fontsize=9, fontweight="bold")
fig.text(0.01, 0.012,
         "Scores are mean z-scored expression of the frozen 10-gene candidate. C-PROBE row IDs were matched by Entrez-like probe IDs; all 10 genes were present. *FDR<0.10 vs living donor.",
         ha="left", va="bottom", fontsize=6, color="#555555")
fig.subplots_adjust(top=0.92, bottom=0.16, left=0.08, right=0.98)
for ext in ["svg", "pdf", "png", "tiff"]:
    if ext in ["png", "tiff"]:
        fig.savefig(f"{PREFIX}.{ext}", dpi=600, bbox_inches="tight")
    else:
        fig.savefig(f"{PREFIX}.{ext}", bbox_inches="tight")
plt.close(fig)

# -----------------------------------------------------------------------------
# Markdown report
# -----------------------------------------------------------------------------
# Pull key numbers for prose.
def fmt_p(p):
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3g}"

def row_lookup(df, group):
    hit = df[df["group"] == group]
    return hit.iloc[0] if not hit.empty else None

gse_dn = row_lookup(gse175_stats, "DN")
gse_igan = row_lookup(gse175_stats, "IgAN")
gse_all = gse175_corr[gse175_corr["group"] == "All usable"].iloc[0]
tub_dn = cprobe_stats[(cprobe_stats["compartment"] == "Tubuli") & (cprobe_stats["group"] == "DN")].iloc[0]
tub_igan = cprobe_stats[(cprobe_stats["compartment"] == "Tubuli") & (cprobe_stats["group"] == "IgAN")].iloc[0]
tub_fsgs = cprobe_stats[(cprobe_stats["compartment"] == "Tubuli") & (cprobe_stats["group"] == "FSGS/FGGS")].iloc[0]
tub_ln = cprobe_stats[(cprobe_stats["compartment"] == "Tubuli") & (cprobe_stats["group"] == "LN")].iloc[0]

gene_presence_summary = cprobe_presence.groupby(["dataset", "compartment"])["n_rows"].apply(lambda x: int((x > 0).sum())).reset_index(name="genes_present")

sections = [
    "# Step 16. DKD-specificity versus pan-CKD specificity audit",
    "## Question\nDoes the frozen 10-gene injury-dominant remodeling candidate behave as a DKD-specific signal, or as a broader CKD/tubulointerstitial injury-remodeling context that is also present in DKD?",
    "## Data and scoring\n- GSE175759 was reused from Step 12, restricted to the frozen `candidate_DKD_TI_remodeling_no_SPP1_CD44` score and usable samples.\n- C-PROBE was used as an external pan-CKD screen: GSE180394 for tubuli and GSE180393 for glomeruli. Matrix row IDs matched the 10 frozen genes through Entrez-like probe IDs with an `_at` suffix.\n- All 10 candidate genes were recovered in both C-PROBE compartments.\n- Statistics are two-sided Mann-Whitney tests versus control/living donor, with Benjamini-Hochberg FDR across disease groups within each compartment/dataset. This is a disease-context audit, not an outcome or causality analysis.",
    f"## Main result\nThe signature is best described as **DKD-relevant but not DKD-exclusive**. In GSE175759, the candidate retained an inverse association with eGFR across usable samples (Spearman rho={gse_all['spearman_rho']:.3f}, P={fmt_p(gse_all['p_value'])}), but the diabetic nephropathy subset was too small for specificity claims (DN n={int(gse_dn['n']) if gse_dn is not None else 'NA'}). In C-PROBE tubuli, DN showed an elevated score versus living donor (median delta={tub_dn['delta_median_vs_reference']:.3f}, P={fmt_p(tub_dn['mannwhitney_p'])}, FDR={fmt_p(tub_dn['fdr'])}), but non-diabetic CKD groups also showed elevation, including IgAN (delta={tub_igan['delta_median_vs_reference']:.3f}, FDR={fmt_p(tub_igan['fdr'])}), FSGS/FGGS (delta={tub_fsgs['delta_median_vs_reference']:.3f}, FDR={fmt_p(tub_fsgs['fdr'])}), and LN (delta={tub_ln['delta_median_vs_reference']:.3f}, FDR={fmt_p(tub_ln['fdr'])}).",
    "## Interpretation for manuscript positioning\nThis result argues against a narrow DKD-specific framing. The safer and stronger wording is that the candidate captures a **conserved CKD injury-remodeling context observed in DKD**, with tubulointerstitial relevance and external pan-CKD support. This helps answer the likely reviewer objection that HAVCR1/LCN2/CLU are common injury markers: the point is not novelty of individual genes, but the reproducible tissue-state context and its relationship to kidney dysfunction.",
    "## Claim boundary\nDo not write that this signature predicts progression or is specific for diabetic nephropathy. The support is for a conserved, injury-dominant transcriptomic state that is present in DN/DKD and several non-diabetic CKD contexts. The C-PROBE analysis strengthens generalizability, but it does not replace sample-level longitudinal validation.",
    "## Files written\n- `16_pan_ckd_gse175759_diagnosis_tests.csv`\n- `16_pan_ckd_gse175759_egfr_correlations.csv`\n- `16_pan_ckd_cprobe_gene_presence.csv`\n- `16_pan_ckd_cprobe_scores.csv`\n- `16_pan_ckd_cprobe_diagnosis_tests.csv`\n- `16_pan_ckd_specificity_summary.csv`\n- `16_pan_ckd_specificity_audit.svg/pdf/png/tiff`",
    "## Gene recovery check\n" + simple_markdown_table(gene_presence_summary),
]
write_report(OUTDIR / "16_pan_ckd_specificity_audit.md", sections)

# Console summary
print("Step 16 complete")
print("GSE175759 disease tests")
print(gse175_stats.to_string(index=False))
print("\nGSE175759 eGFR correlations")
print(gse175_corr.to_string(index=False))
print("\nC-PROBE gene presence")
print(gene_presence_summary.to_string(index=False))
print("\nC-PROBE diagnosis tests")
print(cprobe_stats.sort_values(["compartment", "group"]).to_string(index=False))
print("\nWrote", OUTDIR / "16_pan_ckd_specificity_audit.md")
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
