# -*- coding: utf-8 -*-
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUTDIR.mkdir(parents=True, exist_ok=True)
PREFIX = OUTDIR / "31_main_figure1_signature_audit_reframed"

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

COL_NEW = "#D55E00"
COL_INJ = "#0072B2"
COL_GREEN = "#009E73"
COL_LIGHT = "#E9ECEF"
COL_NULL = "#C9D6E3"
COL_FIB = "#6C757D"
COL_INF = "#7A8793"

cand = pd.read_csv(OUTDIR / "06_frozen_signature_candidate_summary.csv")
nulls = pd.read_csv(OUTDIR / "07_candidate_random_control_nulls.csv.gz")
null_summary = pd.read_csv(OUTDIR / "07_candidate_random_control_summary.csv")
comp = pd.read_csv(OUTDIR / "07_candidate_component_sensitivity.csv")

label_map = {
    "candidate_DKD_TI_remodeling_no_SPP1_CD44": "Frozen 10-gene\ncandidate",
    "human_tubular_injury_without_SPP1_CD44": "Tubular injury\nreference",
    "human_tubular_injury_context": "Tubular injury\ncontext",
    "SPP1_CD44_anchor": "SPP1/CD44\ncell-context anchor",
    "curated_fibro_inflammatory_context": "Curated remodeling\nreference",
    "fibrosis_ecm_reference": "Fibrosis/ECM",
    "inflammation_reference": "Inflammation",
    "NicheNet_SPP1_receiver_target_set": "SPP1 receiver\ntarget set",
    "component_injury_HAVCR1_LCN2_VCAM1": "Injury\nHAVCR1/LCN2/VCAM1",
    "component_CLU": "CLU",
    "component_fibrosis_ECM": "Fibrosis/ECM",
    "component_inflammation": "Inflammation",
    "candidate_minus_injury": "Candidate\nminus injury",
    "candidate_minus_CLU": "Candidate\nminus CLU",
    "candidate_minus_fibrosis_ECM": "Candidate\nminus fibrosis/ECM",
    "candidate_minus_inflammation": "Candidate\nminus inflammation",
}

candidate_name = "candidate_DKD_TI_remodeling_no_SPP1_CD44"
key_sigs = [
    "candidate_DKD_TI_remodeling_no_SPP1_CD44",
    "human_tubular_injury_without_SPP1_CD44",
    "human_tubular_injury_context",
    "SPP1_CD44_anchor",
    "curated_fibro_inflammatory_context",
    "fibrosis_ecm_reference",
    "inflammation_reference",
    "NicheNet_SPP1_receiver_target_set",
]
plot_cand = cand[cand["signature"].isin(key_sigs)].copy()
plot_cand["label"] = plot_cand["signature"].map(label_map)
plot_cand[["signature", "label", "gse175759_rho", "gse175759_p", "gse30122_primary_tubule_delta", "gse30122_primary_tubule_p"]].to_csv(
    OUTDIR / "31_main_figure1_panelB_source_data.csv", index=False, encoding="utf-8-sig"
)
null_summary.to_csv(OUTDIR / "31_main_figure1_panelCD_random_control_summary.csv", index=False, encoding="utf-8-sig")

component_order = [
    "candidate_DKD_TI_remodeling_no_SPP1_CD44",
    "component_injury_HAVCR1_LCN2_VCAM1",
    "component_CLU",
    "component_fibrosis_ECM",
    "component_inflammation",
    "candidate_minus_injury",
    "candidate_minus_CLU",
    "candidate_minus_fibrosis_ECM",
    "candidate_minus_inflammation",
]
metric_specs = [
    ("GSE175759", "eGFR_spearman_all_samples", "-rho\n(lower eGFR)"),
    ("GSE30122", "primary_tubules_DKD_minus_control", "delta primary\ntubules"),
    ("GSE30122", "glomerulus_DKD_minus_control", "delta glomeruli"),
]
mat_rows = []
for sig in component_order:
    for ds, analysis, metric_label in metric_specs:
        hit = comp[(comp["signature"] == sig) & (comp["dataset"].str.strip() == ds) & (comp["analysis"] == analysis)]
        if hit.empty:
            val = np.nan; p = np.nan
        else:
            effect = float(hit.iloc[0]["effect"])
            p = float(hit.iloc[0]["p_value"])
            val = -effect if analysis == "eGFR_spearman_all_samples" else effect
        mat_rows.append({"signature": sig, "label": label_map.get(sig, sig), "metric": metric_label, "support_value": val, "p_value": p})
comp_plot = pd.DataFrame(mat_rows)
comp_plot.to_csv(OUTDIR / "31_main_figure1_panelE_component_source_data.csv", index=False, encoding="utf-8-sig")


def panel_label(ax, label):
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="top")


def draw_box(ax, xy, width, height, text, fc, ec):
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.02,rounding_size=0.025", linewidth=0.9, edgecolor=ec, facecolor=fc)
    ax.add_patch(patch)
    ax.text(x + width/2, y + height/2, text, ha="center", va="center", fontsize=6.4, color="#222222", linespacing=1.12)


def plot_null(ax, dataset, analysis, observed, empirical_p, direction, title, xlabel):
    sub = nulls[(nulls["dataset"].str.strip() == dataset) & (nulls["analysis"] == analysis)].copy()
    vals = sub["random_stat"].dropna().to_numpy()
    bins = np.linspace(vals.min(), vals.max(), 43)
    ax.hist(vals, bins=bins, color=COL_NULL, edgecolor="white", linewidth=0.35)
    ax.axvline(observed, color=COL_NEW, lw=1.8)
    ax.axvline(np.median(vals), color="#586069", lw=1.0, ls="--")
    if direction == "lower_or_equal":
        tail = vals[vals <= observed]
        label_x, label_ha = 0.09, "left"
    else:
        tail = vals[vals >= observed]
        label_x, label_ha = 0.92, "right"
    ax.hist(tail, bins=bins, color=COL_NEW, alpha=0.55, edgecolor="none")
    ax.text(label_x, 0.80, f"observed\n{observed:.3f}", transform=ax.transAxes, ha=label_ha, va="top", color=COL_NEW, fontsize=6.5)
    ax.text(0.02, 0.94, f"empirical P={empirical_p:.4g}\nn={len(vals):,} random sets", transform=ax.transAxes, ha="left", va="top", fontsize=6.5)
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Random sets")
    ax.grid(axis="y", color=COL_LIGHT, lw=0.6)

fig = plt.figure(figsize=(7.35, 7.55))
gs = fig.add_gridspec(3, 2, height_ratios=[0.95, 1.05, 1.15], width_ratios=[1.05, 1.0], hspace=0.58, wspace=0.38)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])
axD = fig.add_subplot(gs[2, 0])
axE = fig.add_subplot(gs[2, 1])

# A
axA.axis("off")
panel_label(axA, "A")
axA.set_title("Candidate definition and robustness audit", loc="left", fontsize=8.5, fontweight="bold", pad=2)
draw_box(axA, (0.02, 0.22), 0.22, 0.48, "Biological scope\ntubular injury/stress\nmatrix remodeling\ninflammatory sensors", "#F1F3F5", "#8A8F98")
draw_box(axA, (0.29, 0.22), 0.22, 0.48, "Frozen 10-gene candidate\nHAVCR1, LCN2, VCAM1, CLU\nCOL1A1, FN1, ACTA2\nTGFB1, TLR4, NLRP3", "#FFF3E8", COL_NEW)
draw_box(axA, (0.56, 0.22), 0.19, 0.48, "Robustness audit\neGFR rho=-0.227\ntubule delta=0.789\nemp. P=0.0002 / 0.0014", "#EAF4FB", COL_INJ)
draw_box(axA, (0.80, 0.22), 0.18, 0.48, "Interpretation\ninjury/CLU-dominant\ntubulointerstitial\nremodeling context", "#EAF6EF", COL_GREEN)
for start, end in [((0.245, 0.46), (0.285, 0.46)), ((0.515, 0.46), (0.555, 0.46)), ((0.755, 0.46), (0.795, 0.46))]:
    axA.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10, lw=1.0, color="#333333"))
axA.text(0.02, 0.08, "SPP1 and CD44 were reserved for cellular-context analyses and were not included in the candidate score.", ha="left", va="center", fontsize=6.5, color="#555555")

# B
panel_label(axB, "B")
axB.axvline(0, color="#CCCCCC", lw=0.8)
axB.axhline(0, color="#CCCCCC", lw=0.8)
axB.axvspan(-0.34, 0, ymin=0.5, ymax=1, color="#FFF5EB", alpha=0.9, zorder=0)
for _, r in plot_cand.iterrows():
    sig = r["signature"]
    if sig == candidate_name:
        color, size, z = COL_NEW, 74, 5
    elif "injury" in sig:
        color, size, z = COL_INJ, 52, 4
    elif sig == "fibrosis_ecm_reference":
        color, size, z = COL_FIB, 40, 3
    elif sig == "inflammation_reference":
        color, size, z = COL_INF, 40, 3
    else:
        color, size, z = "#9AA7B2", 35, 3
    axB.scatter(r["gse175759_rho"], r["gse30122_primary_tubule_delta"], s=size, color=color, edgecolor="white", linewidth=0.6, zorder=z)
label_offsets = {
    candidate_name: (0.012, 0.08),
    "human_tubular_injury_without_SPP1_CD44": (-0.052, -0.04),
    "fibrosis_ecm_reference": (0.010, 0.05),
    "inflammation_reference": (-0.030, -0.20),
    "SPP1_CD44_anchor": (0.012, 0.10),
}
for sig, (dx, dy) in label_offsets.items():
    r = plot_cand[plot_cand["signature"] == sig]
    if r.empty:
        continue
    r = r.iloc[0]
    axB.text(r["gse175759_rho"] + dx, r["gse30122_primary_tubule_delta"] + dy, label_map.get(sig, sig), fontsize=5.8, color="#222222")
axB.set_xlabel("GSE175759 Spearman rho vs eGFR\n(more negative = higher score with lower eGFR)")
axB.set_ylabel("GSE30122 DKD-control delta\nprimary tubules")
axB.set_title("Candidate and reference signatures across human contexts", loc="left", fontsize=8, fontweight="bold")
axB.set_xlim(-0.32, 0.02)
axB.set_ylim(-0.62, 1.88)
axB.grid(color=COL_LIGHT, lw=0.6)

# C/D
summary1 = null_summary[(null_summary["dataset"].str.strip() == "GSE175759") & (null_summary["analysis"] == "eGFR_spearman_all_samples")].iloc[0]
panel_label(axC, "C")
plot_null(axC, "GSE175759", "eGFR_spearman_all_samples", float(summary1["observed_stat"]), float(summary1["empirical_directional_p"]), str(summary1["direction_tested"]), "Expression-matched control: eGFR association", "Spearman rho vs eGFR")
summary2 = null_summary[(null_summary["dataset"].str.strip() == "GSE30122") & (null_summary["analysis"] == "primary_tubules_DKD_minus_control")].iloc[0]
panel_label(axD, "D")
plot_null(axD, "GSE30122", "primary_tubules_DKD_minus_control", float(summary2["observed_stat"]), float(summary2["empirical_directional_p"]), str(summary2["direction_tested"]), "Expression-matched control: DKD primary tubules", "DKD-control delta")

# E
panel_label(axE, "E")
mat = comp_plot.pivot(index="label", columns="metric", values="support_value")
row_labels = [label_map.get(s, s) for s in component_order]
col_labels = [m[2] for m in metric_specs]
mat = mat.reindex(row_labels)[col_labels]
values = mat.to_numpy(dtype=float)
max_abs = np.nanmax(np.abs(values))
im = axE.imshow(values, cmap="RdBu_r", vmin=-max_abs, vmax=max_abs, aspect="auto")
axE.set_xticks(np.arange(len(col_labels)))
axE.set_xticklabels(col_labels, rotation=0, ha="center")
axE.set_yticks(np.arange(len(row_labels)))
axE.set_yticklabels(row_labels, fontsize=5.8)
for i in range(values.shape[0]):
    for j in range(values.shape[1]):
        val = values[i, j]
        if np.isfinite(val):
            color = "white" if abs(val) > max_abs*0.55 else "#222222"
            axE.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color=color)
axE.set_title("Component sensitivity supports injury/CLU dominance", loc="left", fontsize=8, fontweight="bold")
for y in np.arange(-0.5, len(row_labels), 1):
    axE.axhline(y, color="white", lw=0.7)
for x in np.arange(-0.5, len(col_labels), 1):
    axE.axvline(x, color="white", lw=0.7)
cbar = fig.colorbar(im, ax=axE, fraction=0.046, pad=0.02)
cbar.set_label("Support for candidate direction", fontsize=6.5)
cbar.ax.tick_params(labelsize=6)

fig.suptitle("Figure 1. Signature definition and robustness audit", x=0.01, y=0.995, ha="left", fontsize=10, fontweight="bold")
fig.text(0.01, 0.008, "Scores are mean gene-wise z scores. Random controls are expression-matched 10-gene sets. Negative rho means higher score with lower eGFR; positive delta means higher score in DKD.", ha="left", va="bottom", fontsize=6.3, color="#555555")
fig.subplots_adjust(top=0.94, bottom=0.08, left=0.08, right=0.98)

for ext in ["svg", "pdf", "png", "tiff"]:
    if ext in ["png", "tiff"]:
        fig.savefig(f"{PREFIX}.{ext}", dpi=600, bbox_inches="tight")
    else:
        fig.savefig(f"{PREFIX}.{ext}", bbox_inches="tight")
plt.close(fig)

legend = f"""# Main Figure 1. Signature definition and robustness audit

## Core figure claim

A frozen 10-gene candidate captures an injury-dominant tubulointerstitial remodeling context that shows the expected kidney-function and DKD primary-tubule directions, exceeds expression-matched random controls, and is dominated by injury/CLU components.

## Panel legend draft

**Figure 1 | Signature definition and robustness audit identifies an injury-dominant remodeling candidate.**
**A,** Candidate definition and robustness workflow. The frozen candidate contains HAVCR1, LCN2, VCAM1, CLU, COL1A1, FN1, ACTA2, TGFB1, TLR4 and NLRP3. SPP1 and CD44 were reserved for cellular-context analyses and were not included in the candidate score.
**B,** Candidate and reference signatures were compared across a human tubulointerstitial kidney-function context (GSE175759 Spearman rho versus eGFR) and a DKD microdissected primary-tubule context (GSE30122 DKD-control delta). The frozen candidate occupied the expected direction of higher score with lower eGFR and higher score in DKD primary tubules.
**C,** Expression-matched random-control audit for the GSE175759 eGFR association. The observed candidate statistic fell beyond the lower tail of 5,000 matched random 10-gene sets (observed rho={float(summary1['observed_stat']):.3f}; empirical directional P={float(summary1['empirical_directional_p']):.4g}).
**D,** Expression-matched random-control audit for the GSE30122 DKD-control primary-tubule contrast. The observed candidate statistic exceeded the upper tail of 5,000 matched random 10-gene sets (observed delta={float(summary2['observed_stat']):.3f}; empirical directional P={float(summary2['empirical_directional_p']):.4g}).
**E,** Component sensitivity analysis. Values are oriented so that larger positive values indicate stronger support for the candidate direction. Injury markers and CLU contributed the strongest support, whereas the inflammatory component did not support the primary-tubule DKD direction.

## Claim boundary

This figure supports selection of a conserved injury-dominant remodeling candidate. It does not establish a DKD-specific biomarker or a validated prognostic score.

## Files

- `31_main_figure1_signature_audit_reframed.svg`
- `31_main_figure1_signature_audit_reframed.pdf`
- `31_main_figure1_signature_audit_reframed.png`
- `31_main_figure1_signature_audit_reframed.tiff`
- `31_main_figure1_panelB_source_data.csv`
- `31_main_figure1_panelCD_random_control_summary.csv`
- `31_main_figure1_panelE_component_source_data.csv`
"""
(OUTDIR / "31_main_figure1_signature_audit_reframed_legend.md").write_text(legend, encoding="utf-8")
print("WROTE", PREFIX.with_suffix(".png"))
print("WROTE", OUTDIR / "31_main_figure1_signature_audit_reframed_legend.md")
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
