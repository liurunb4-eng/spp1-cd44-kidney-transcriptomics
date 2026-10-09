
# -*- coding: utf-8 -*-
"""
Main Figure 2 assembly: kidney-function association and pan-CKD disease-context audit.

Core conclusion:
The frozen 10-gene injury-dominant candidate is associated with lower baseline
eGFR in a human tubulointerstitial CKD cohort and is elevated in C-PROBE tubuli
from diabetic nephropathy and several non-diabetic CKD groups. This supports a
DKD-relevant but not DKD-exclusive conserved CKD injury-remodeling context.

Archetype: quantitative grid.
Backend: Python/matplotlib only.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

OUTDIR = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUTDIR.mkdir(parents=True, exist_ok=True)
PREFIX = OUTDIR / "20_main_figure2_clinical_pan_ckd_audit"
SIG = "candidate_DKD_TI_remodeling_no_SPP1_CD44"

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

COL_CONTROL = "#8A8F98"
COL_CKD = "#0072B2"
COL_DN = "#D55E00"
COL_IgAN = "#4C78A8"
COL_FSGS = "#009E73"
COL_LN = "#CC79A7"
COL_OTHER = "#6C757D"
COL_LIGHT = "#E9ECEF"
COL_TUB = "#4C78A8"
COL_GLOM = "#F58518"

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
scores_all = pd.read_csv(OUTDIR / "12_clinical_egfr_stratification_scores_merged.csv")
gse = scores_all[(scores_all["signature"] == SIG) & (scores_all["usable"] == True)].copy()  # noqa: E712
gse["egfr_ckd_epi"] = pd.to_numeric(gse["egfr_ckd_epi"], errors="coerce")
gse["score"] = pd.to_numeric(gse["score"], errors="coerce")
# Compact disease grouping for scatter colors.
gse["plot_group"] = np.where(gse["diagnosis"].eq("Control"), "Control", "CKD")
gse.loc[gse["diagnosis"].eq("Diabetic nephropathy"), "plot_group"] = "DN"

tests = pd.read_csv(OUTDIR / "12_clinical_egfr_stratification_tests.csv")
corrs = pd.read_csv(OUTDIR / "12_clinical_egfr_stratification_correlations.csv")
tert_test = tests[(tests["signature"] == SIG) & (tests["contrast"] == "low_vs_high_egfr_tertile")].iloc[0]
all_corr = corrs[(corrs["signature"] == SIG) & (corrs["subgroup"] == "all_usable")].iloc[0]
ckd_corr = corrs[(corrs["signature"] == SIG) & (corrs["subgroup"] == "ckd_only")].iloc[0]

cprobe = pd.read_csv(OUTDIR / "16_pan_ckd_cprobe_scores.csv")
cprobe_stats = pd.read_csv(OUTDIR / "16_pan_ckd_cprobe_diagnosis_tests.csv")

# Source data exports.
gse.to_csv(OUTDIR / "20_main_figure2_panelAB_source_data.csv", index=False, encoding="utf-8-sig")
cprobe.to_csv(OUTDIR / "20_main_figure2_panelC_source_data.csv", index=False, encoding="utf-8-sig")
cprobe_stats.to_csv(OUTDIR / "20_main_figure2_panelD_source_data.csv", index=False, encoding="utf-8-sig")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def panel_label(ax, label):
    ax.text(-0.10, 1.06, label, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="top")


def jitter(n, width=0.16):
    if n <= 1:
        return np.array([0.0])[:n]
    return np.linspace(-width, width, n)


def median_iqr(ax, x, vals, color="#222222", width=0.28):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    ax.plot([x-width, x+width], [med, med], color=color, lw=1.25, zorder=4)
    ax.vlines(x, q1, q3, color=color, lw=1.0, zorder=4)


def p_label(p):
    if pd.isna(p):
        return "NA"
    p = float(p)
    return f"{p:.1e}" if p < 0.001 else f"{p:.3g}"

# -----------------------------------------------------------------------------
# Figure layout
# -----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.35, 6.75))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.0], hspace=0.42, wspace=0.34)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# Panel A: eGFR tertiles
panel_label(axA, "A")
order_tert = ["High eGFR", "Low eGFR"]
colors_tert = {"High eGFR": "#9AA1AA", "Low eGFR": COL_DN}
for i, grp in enumerate(order_tert):
    vals = gse.loc[gse["egfr_tertile"] == grp, "score"].dropna().sort_values().to_numpy()
    xs = i + jitter(len(vals), 0.18)
    axA.scatter(xs, vals, s=20, color=colors_tert[grp], edgecolor="white", linewidth=0.35, alpha=0.85, zorder=3)
    median_iqr(axA, i, vals)
axA.axhline(0, color="#BBBBBB", lw=0.8)
axA.set_xticks(range(len(order_tert)))
axA.set_xticklabels([f"{g}\nn={(gse['egfr_tertile']==g).sum()}" for g in order_tert])
axA.set_ylabel("Candidate score (z-mean)")
axA.set_title("Higher candidate score in lower-eGFR samples", loc="left", fontsize=8, fontweight="bold")
axA.text(0.02, 0.96, f"Low vs high: delta={tert_test['mean_delta']:.3f}\nFDR={tert_test['mannwhitney_fdr_within_contrast']:.3g}", transform=axA.transAxes, ha="left", va="top", fontsize=6.5)
axA.grid(axis="y", color=COL_LIGHT, lw=0.6)

# Panel B: continuous eGFR association
panel_label(axB, "B")
plot_colors = {"Control": COL_CONTROL, "CKD": COL_CKD, "DN": COL_DN}
for grp in ["Control", "CKD", "DN"]:
    sub = gse[gse["plot_group"] == grp]
    if sub.empty:
        continue
    axB.scatter(sub["score"], sub["egfr_ckd_epi"], s=22 if grp != "DN" else 30, color=plot_colors[grp], alpha=0.82, edgecolor="white", linewidth=0.35, label=f"{grp} (n={len(sub)})", zorder=4 if grp == "DN" else 3)
mask = gse["score"].notna() & gse["egfr_ckd_epi"].notna()
coef = np.polyfit(gse.loc[mask, "score"], gse.loc[mask, "egfr_ckd_epi"], deg=1)
xs = np.linspace(gse.loc[mask, "score"].min(), gse.loc[mask, "score"].max(), 100)
axB.plot(xs, coef[0]*xs + coef[1], color="#222222", lw=1.0)
axB.set_xlabel("Candidate score")
axB.set_ylabel("eGFR CKD-EPI")
axB.set_title("Continuous kidney-function association", loc="left", fontsize=8, fontweight="bold")
axB.text(0.02, 0.05, f"all samples: rho={all_corr['spearman_rho']:.3f}, FDR={all_corr['spearman_fdr_within_subgroup']:.3g}\nCKD only: rho={ckd_corr['spearman_rho']:.3f}, FDR={ckd_corr['spearman_fdr_within_subgroup']:.3g}", transform=axB.transAxes, ha="left", va="bottom", fontsize=6.5)
axB.legend(loc="upper right", fontsize=6, handletextpad=0.2)
axB.grid(color=COL_LIGHT, lw=0.6)

# Panel C: C-PROBE tubuli disease-context audit
panel_label(axC, "C")
tub = cprobe[cprobe["compartment"] == "Tubuli"].copy()
order_c = ["Living donor", "DN", "IgAN", "FSGS/FGGS", "LN", "Other CKD"]
colors_c = {"Living donor": COL_CONTROL, "DN": COL_DN, "IgAN": COL_IgAN, "FSGS/FGGS": COL_FSGS, "LN": COL_LN, "Other CKD": COL_OTHER}
for i, grp in enumerate(order_c):
    vals = tub.loc[tub["group_collapsed"] == grp, "score"].dropna().sort_values().to_numpy()
    xs = i + jitter(len(vals), 0.18)
    axC.scatter(xs, vals, s=22, color=colors_c[grp], edgecolor="white", linewidth=0.35, alpha=0.85, zorder=3)
    median_iqr(axC, i, vals)
axC.axhline(0, color="#BBBBBB", lw=0.8)
axC.axvspan(-0.45, 0.45, color="#F1F3F5", zorder=0)
axC.set_xticks(range(len(order_c)))
axC.set_xticklabels([f"{g}\nn={(tub['group_collapsed']==g).sum()}" for g in order_c], rotation=30, ha="right")
axC.set_ylabel("Candidate score (z-mean)")
axC.set_title("C-PROBE tubuli: DKD-relevant but not DKD-exclusive", loc="left", fontsize=8, fontweight="bold")
# Mark groups significant versus living donor; exact FDRs are in source data and legend.
for i, grp in enumerate(order_c[1:], start=1):
    hit = cprobe_stats[(cprobe_stats["compartment"] == "Tubuli") & (cprobe_stats["group"] == grp)]
    if not hit.empty and np.isfinite(hit.iloc[0]["fdr"]) and hit.iloc[0]["fdr"] < 0.05:
        vals = tub.loc[tub["group_collapsed"] == grp, "score"].dropna()
        y = vals.max() + 0.08 if len(vals) else 0.2
        axC.text(i, y, "*", ha="center", va="bottom", fontsize=9, color="#333333")
axC.text(0.02, 0.04, "* FDR<0.05 vs living donor", transform=axC.transAxes, ha="left", va="bottom", fontsize=6.2, color="#555555")
axC.grid(axis="y", color=COL_LIGHT, lw=0.6)

# Panel D: C-PROBE compartment contrast dot plot
panel_label(axD, "D")
order_d = ["DN", "IgAN", "FSGS/FGGS", "LN", "Other CKD", "Other GN/immune"]
ypos = np.arange(len(order_d))[::-1]
for comp_name, color, offset in [("Tubuli", COL_TUB, 0.12), ("Glomeruli", COL_GLOM, -0.12)]:
    xs_plot = []
    ys_plot = []
    labels = []
    for y, grp in zip(ypos, order_d):
        hit = cprobe_stats[(cprobe_stats["compartment"] == comp_name) & (cprobe_stats["group"] == grp)]
        if hit.empty or not np.isfinite(hit.iloc[0]["delta_median_vs_reference"]):
            continue
        row = hit.iloc[0]
        xs_plot.append(float(row["delta_median_vs_reference"]))
        ys_plot.append(y + offset)
        labels.append(row)
    axD.scatter(xs_plot, ys_plot, s=48, color=color, edgecolor="white", linewidth=0.5, label=comp_name, zorder=3)
    for xval, yval, row in zip(xs_plot, ys_plot, labels):
        if np.isfinite(row["fdr"]) and row["fdr"] < 0.05:
            axD.text(xval + 0.035, yval, "*", ha="left", va="center", fontsize=9, color=color)
axD.axvline(0, color="#BBBBBB", lw=0.8)
axD.set_yticks(ypos)
axD.set_yticklabels(order_d)
axD.set_xlabel("Median delta vs living donor")
axD.set_title("Compartment contrast confirms a pan-CKD context", loc="left", fontsize=8, fontweight="bold")
axD.legend(loc="lower right", fontsize=6)
axD.grid(axis="x", color=COL_LIGHT, lw=0.6)
axD.text(0.02, 0.04, "* FDR<0.05 vs living donor", transform=axD.transAxes, ha="left", va="bottom", fontsize=6.2, color="#555555")

fig.suptitle("Figure 2. Kidney-function association and pan-CKD disease-context audit", x=0.01, y=0.995, ha="left", fontsize=10, fontweight="bold")
fig.text(0.01, 0.010, "Candidate score is the mean z-scored expression of the frozen 10-gene injury-dominant candidate. C-PROBE analyses compare each disease group with living donor samples; rare groups are shown in source data.", ha="left", va="bottom", fontsize=6.2, color="#555555")
fig.subplots_adjust(top=0.93, bottom=0.14, left=0.08, right=0.98)

for ext in ["svg", "pdf", "png", "tiff"]:
    if ext in ["png", "tiff"]:
        fig.savefig(f"{PREFIX}.{ext}", dpi=600, bbox_inches="tight")
    else:
        fig.savefig(f"{PREFIX}.{ext}", bbox_inches="tight")
plt.close(fig)

legend = f"""# Main Figure 2. Kidney-function association and pan-CKD disease-context audit

## Core figure claim

The frozen 10-gene injury-dominant candidate is associated with lower baseline eGFR in GSE175759 and is elevated in C-PROBE tubuli from diabetic nephropathy and several non-diabetic CKD groups. This supports a DKD-relevant but not DKD-exclusive conserved CKD injury-remodeling context.

## Panel legend draft

**Figure 2 | Kidney-function association and pan-CKD disease-context audit.**
**A,** Candidate scores across GSE175759 eGFR tertiles. The low-eGFR tertile showed a higher candidate score than the high-eGFR tertile (mean-score delta={tert_test['mean_delta']:.3f}; Mann-Whitney FDR={tert_test['mannwhitney_fdr_within_contrast']:.3g}).
**B,** Continuous association between candidate score and eGFR in GSE175759. Across all usable samples, higher candidate score was associated with lower eGFR (Spearman rho={all_corr['spearman_rho']:.3f}; FDR={all_corr['spearman_fdr_within_subgroup']:.3g}); the CKD-only sensitivity analysis showed the same direction (rho={ckd_corr['spearman_rho']:.3f}; FDR={ckd_corr['spearman_fdr_within_subgroup']:.3g}).
**C,** C-PROBE tubuli disease-context audit. Diabetic nephropathy showed elevated candidate score relative to living donor samples, but non-diabetic CKD groups including IgA nephropathy, FSGS/FGGS, lupus nephritis and other CKD groups were also elevated.
**D,** C-PROBE compartment contrast. Points show median disease-group delta versus living donor for tubuli and glomeruli; asterisks indicate FDR<0.05 versus living donor.

## Claim boundary

This figure supports cross-sectional kidney-function association and pan-CKD disease-context generalizability. It does not establish individual-level longitudinal prediction, diagnostic specificity, or a DKD-specific biomarker.

## Files

- `20_main_figure2_clinical_pan_ckd_audit.svg`
- `20_main_figure2_clinical_pan_ckd_audit.pdf`
- `20_main_figure2_clinical_pan_ckd_audit.png`
- `20_main_figure2_clinical_pan_ckd_audit.tiff`
- `20_main_figure2_panelAB_source_data.csv`
- `20_main_figure2_panelC_source_data.csv`
- `20_main_figure2_panelD_source_data.csv`
"""
(OUTDIR / "20_main_figure2_clinical_pan_ckd_audit_legend.md").write_text(legend, encoding="utf-8")

print("Figure 2 assembly complete")
print(PREFIX.with_suffix(".png"))
print(OUTDIR / "20_main_figure2_clinical_pan_ckd_audit_legend.md")
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
