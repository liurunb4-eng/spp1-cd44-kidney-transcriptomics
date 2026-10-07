# -*- coding: utf-8 -*-
"""Main Figure 3: GSE137570 patient-level CKD progression and clinicopathologic support.

Core conclusion:
The frozen 10-gene injury-remodeling candidate is higher in patient-level
progressive CKD and aligns with lower GFR and greater tubulointerstitial
fibrosis in an independent CKD biopsy transcriptomic dataset.
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

OUTDIR = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUTDIR.mkdir(parents=True, exist_ok=True)
DATADIR = OUTDIR / "35_gse137570_ckd_progression_screen"
PREFIX = OUTDIR / "37_main_figure3_gse137570_patient_level_support"

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

COL_STABLE = "#6C757D"
COL_PROG = "#D55E00"
COL_BLUE = "#0072B2"
COL_GREEN = "#009E73"
COL_LIGHT = "#E9ECEF"
COL_TEXT = "#222222"

scores2 = pd.read_csv(DATADIR / "GSE137570_cohort2_candidate_scores.csv")
scores1 = pd.read_csv(DATADIR / "GSE137570_cohort1_candidate_scores.csv")
comp_tests = pd.read_csv(DATADIR / "GSE137570_cohort2_component_gene_progression_tests.csv")
cor = pd.read_csv(DATADIR / "GSE137570_cohort1_candidate_correlations.csv")

# Source data exports for final figure.
scores2.to_csv(OUTDIR / "37_main_figure3_panelA_source_data.csv", index=False, encoding="utf-8-sig")
scores1[["sample", "candidate_score", "GFR", "TIF"]].to_csv(OUTDIR / "37_main_figure3_panelBC_source_data.csv", index=False, encoding="utf-8-sig")
keep_order = [
    "candidate_10_gene",
    "injury_HAVCR1_LCN2_VCAM1",
    "fibrosis_ECM_COL1A1_FN1_ACTA2_TGFB1",
    "inflammation_TLR4_NLRP3",
    "CLU",
]
comp_plot = comp_tests[comp_tests["score_name"].isin(keep_order)].copy()
comp_plot["score_name"] = pd.Categorical(comp_plot["score_name"], categories=keep_order, ordered=True)
comp_plot = comp_plot.sort_values("score_name")
comp_plot.to_csv(OUTDIR / "37_main_figure3_panelD_source_data.csv", index=False, encoding="utf-8-sig")


def panel_label(ax, label):
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="top")


def mean_line(ax, x, y, width=0.22):
    ax.plot([x-width, x+width], [np.mean(y), np.mean(y)], color="black", lw=1.4, solid_capstyle="round", zorder=4)

fig = plt.figure(figsize=(7.4, 5.25))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.95], width_ratios=[1.05, 1.05, 1.05], hspace=0.58, wspace=0.45)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, :])

# A: progression groups.
panel_label(axA, "A")
order = [("non-progressive", "Non-progressive", COL_STABLE), ("progressive", "Progressive", COL_PROG)]
for i, (key, label, color) in enumerate(order):
    y = scores2.loc[scores2["progression_status"] == key, "candidate_score"].dropna().to_numpy()
    rng = np.random.default_rng(300 + i)
    x = np.full(len(y), i) + rng.uniform(-0.07, 0.07, len(y))
    axA.scatter(x, y, s=32, color=color, edgecolor="white", linewidth=0.6, zorder=3)
    mean_line(axA, i, y)
axA.set_xticks([0, 1])
axA.set_xticklabels(["Non-progressive", "Progressive"], rotation=0)
axA.set_ylabel("10-gene candidate score")
axA.set_title("Cohort 2 progression grouping", loc="left", fontsize=8, fontweight="bold")
axA.grid(axis="y", color=COL_LIGHT, lw=0.6)
main = comp_tests[comp_tests["score_name"] == "candidate_10_gene"].iloc[0]
axA.text(0.03, 0.96, f"n={int(main['n_nonprogressive'])} vs {int(main['n_progressive'])}\ndelta={main['delta_mean_progressive_minus_nonprogressive']:.2f}\nP={main['mannwhitney_p']:.3g}", transform=axA.transAxes, ha="left", va="top", fontsize=6.8)
axA.set_xlim(-0.35, 1.35)

# B: GFR.
panel_label(axB, "B")
sub = scores1[["GFR", "candidate_score"]].dropna()
axB.scatter(sub["GFR"], sub["candidate_score"], color=COL_BLUE, s=32, edgecolor="white", linewidth=0.6, zorder=3)
coef = np.polyfit(sub["GFR"], sub["candidate_score"], 1)
xs = np.linspace(sub["GFR"].min(), sub["GFR"].max(), 100)
axB.plot(xs, coef[0]*xs + coef[1], color="#333333", lw=1.0)
hit = cor[cor["variable"] == "GFR"].iloc[0]
axB.text(0.03, 0.96, f"n={int(hit['n'])}\nrho={hit['spearman_rho']:.2f}\nP={hit['p_value']:.3g}", transform=axB.transAxes, ha="left", va="top", fontsize=6.8)
axB.set_xlabel("GFR")
axB.set_ylabel("10-gene candidate score")
axB.set_title("Cohort 1 kidney function", loc="left", fontsize=8, fontweight="bold")
axB.grid(color=COL_LIGHT, lw=0.6)

# C: TIF.
panel_label(axC, "C")
sub = scores1[["TIF", "candidate_score"]].dropna()
axC.scatter(sub["TIF"], sub["candidate_score"], color=COL_GREEN, s=32, edgecolor="white", linewidth=0.6, zorder=3)
coef = np.polyfit(sub["TIF"], sub["candidate_score"], 1)
xs = np.linspace(sub["TIF"].min(), sub["TIF"].max(), 100)
axC.plot(xs, coef[0]*xs + coef[1], color="#333333", lw=1.0)
hit = cor[cor["variable"] == "TIF"].iloc[0]
axC.text(0.03, 0.96, f"n={int(hit['n'])}\nrho={hit['spearman_rho']:.2f}\nP={hit['p_value']:.3g}", transform=axC.transAxes, ha="left", va="top", fontsize=6.8)
axC.set_xlabel("Tubulointerstitial fibrosis (%)")
axC.set_ylabel("10-gene candidate score")
axC.set_title("Cohort 1 histology", loc="left", fontsize=8, fontweight="bold")
axC.grid(color=COL_LIGHT, lw=0.6)

# D: components.
panel_label(axD, "D")
display_labels = {
    "candidate_10_gene": "10-gene candidate",
    "injury_HAVCR1_LCN2_VCAM1": "Injury\nHAVCR1/LCN2/VCAM1",
    "fibrosis_ECM_COL1A1_FN1_ACTA2_TGFB1": "Fibrosis/ECM\nCOL1A1/FN1/ACTA2/TGFB1",
    "inflammation_TLR4_NLRP3": "Inflammation\nTLR4/NLRP3",
    "CLU": "CLU",
}
plot = comp_plot.copy()
x = np.arange(len(plot))
vals = plot["delta_mean_progressive_minus_nonprogressive"].to_numpy()
axD.bar(x, vals, color=[COL_PROG if v >= 0 else COL_STABLE for v in vals], width=0.68)
axD.axhline(0, color="#333333", lw=0.8)
for i, (_, r) in enumerate(plot.iterrows()):
    axD.text(i, r["delta_mean_progressive_minus_nonprogressive"] + 0.04, f"FDR={r['bh_fdr']:.3g}", ha="center", va="bottom", fontsize=6.4)
axD.set_xticks(x)
axD.set_xticklabels([display_labels[str(s)] for s in plot["score_name"]], rotation=0)
axD.set_ylabel("Mean score delta\n(progressive - non-progressive)")
axD.set_title("Cohort 2 component-level progression contrasts", loc="left", fontsize=8, fontweight="bold")
axD.grid(axis="y", color=COL_LIGHT, lw=0.6)
axD.set_ylim(0, max(vals) + 0.35)

fig.suptitle("Figure 3. Patient-level CKD progression and clinicopathologic support in GSE137570", x=0.01, y=0.995, ha="left", fontsize=10, fontweight="bold")
fig.text(0.01, 0.008, "Scores are mean gene-wise z scores. P values in A and D are two-sided Mann-Whitney tests; D shows BH FDR. B and C show Spearman correlations.", ha="left", va="bottom", fontsize=6.3, color="#555555")
fig.subplots_adjust(top=0.90, bottom=0.16, left=0.08, right=0.98)

for ext in ["svg", "pdf", "png", "tiff"]:
    out = f"{PREFIX}.{ext}"
    if ext in ["png", "tiff"]:
        fig.savefig(out, dpi=600, bbox_inches="tight")
    else:
        fig.savefig(out, bbox_inches="tight")
plt.close(fig)

legend = f"""# Main Figure 3. Patient-level CKD progression and clinicopathologic support in GSE137570

## Core figure claim

The 10-gene injury-remodeling candidate is higher in patient-level progressive CKD and aligns with lower kidney function and greater tubulointerstitial fibrosis in an independent CKD biopsy transcriptomic dataset.

## Panel legend draft

**Figure 3 | Patient-level CKD progression and clinicopathologic support in GSE137570.**
**A,** Candidate scores in GSE137570 Cohort 2 non-progressive and progressive CKD samples. Points represent individual biopsy samples; horizontal bars indicate group means. The candidate score was higher in progressive CKD (n = {int(main['n_progressive'])}) than non-progressive CKD (n = {int(main['n_nonprogressive'])}; mean-score delta = {main['delta_mean_progressive_minus_nonprogressive']:.2f}; two-sided Mann-Whitney P = {main['mannwhitney_p']:.4g}).
**B,** Association between candidate score and GFR in Cohort 1 (n = {int(cor[cor['variable']=='GFR'].iloc[0]['n'])}; Spearman rho = {cor[cor['variable']=='GFR'].iloc[0]['spearman_rho']:.3f}; P = {cor[cor['variable']=='GFR'].iloc[0]['p_value']:.4g}).
**C,** Association between candidate score and tubulointerstitial fibrosis in Cohort 1 (n = {int(cor[cor['variable']=='TIF'].iloc[0]['n'])}; Spearman rho = {cor[cor['variable']=='TIF'].iloc[0]['spearman_rho']:.3f}; P = {cor[cor['variable']=='TIF'].iloc[0]['p_value']:.4g}).
**D,** Component-level progression contrasts in Cohort 2. Bars show mean-score deltas for progressive minus non-progressive CKD; labels show Benjamini-Hochberg FDR from two-sided Mann-Whitney tests. Source data are provided as Source Data files.

## Claim boundary

This figure supports patient-level progression-group and clinicopathologic alignment for a conserved CKD injury-remodeling context. It does not establish an individualized prognostic model or time-to-event prediction.

## Files

- `37_main_figure3_gse137570_patient_level_support.svg`
- `37_main_figure3_gse137570_patient_level_support.pdf`
- `37_main_figure3_gse137570_patient_level_support.png`
- `37_main_figure3_gse137570_patient_level_support.tiff`
- `37_main_figure3_panelA_source_data.csv`
- `37_main_figure3_panelBC_source_data.csv`
- `37_main_figure3_panelD_source_data.csv`
"""
(OUTDIR / "37_main_figure3_gse137570_patient_level_support_legend.md").write_text(legend, encoding="utf-8")

print("WROTE", PREFIX.with_suffix(".png"))
print("WROTE", OUTDIR / "37_main_figure3_gse137570_patient_level_support_legend.md")
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
