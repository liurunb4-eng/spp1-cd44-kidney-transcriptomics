
# -*- coding: utf-8 -*-
"""
Main Figure 4 assembly: external DKD progression-context gene-level support.

Core conclusion:
External Acoba 2025 DKD progression differential-expression tables support parts
of the 10-gene injury-remodeling candidate at gene level, especially in
 tubulointerstitial CKD-stage advancement, but they do not constitute
sample-level prognostic validation.

Archetype: quantitative grid.
Backend: Python/matplotlib only.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

OUTDIR = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUTDIR.mkdir(parents=True, exist_ok=True)
PREFIX = OUTDIR / "45_main_figure4_acoba_progression_context_bhq"

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

COL_UP = "#C95A5A"
COL_DOWN = "#5B8CC0"
COL_MISS = "#DADDE2"
COL_LIGHT = "#E9ECEF"
COL_ACCENT = "#D55E00"
COL_GREY = "#6C757D"

GENES = ["COL1A1", "HAVCR1", "FN1", "ACTA2", "CLU", "VCAM1", "LCN2", "TGFB1", "TLR4", "NLRP3"]
ENDPOINTS = ["CKD_stage_advancement", "steep_eGFR_slope", "composite_outcome", "UACR_increase", "NRA"]
ENDPOINT_LABELS = {
    "CKD_stage_advancement": "CKD stage\nadvancement",
    "steep_eGFR_slope": "Steep eGFR\nslope",
    "composite_outcome": "Composite\noutcome",
    "UACR_increase": "UACR\nincrease",
    "NRA": "NRA",
}

align = pd.read_csv(OUTDIR / "13_acoba2025_signature_gene_progression_alignment.csv")
summary = pd.read_csv(OUTDIR / "13_acoba2025_signature_progression_summary.csv")
enrich = pd.read_csv(OUTDIR / "13_acoba2025_signature_progression_enrichment.csv")

# Use tubulointerstitium for main progression-context support.
ti = align[(align["signature"] == "candidate_10_gene") & (align["compartment"] == "tubulointerstitium")].copy()
ti["gene"] = pd.Categorical(ti["gene"], categories=GENES, ordered=True)
ti["endpoint"] = pd.Categorical(ti["endpoint"], categories=ENDPOINTS, ordered=True)

heat = ti.pivot(index="gene", columns="endpoint", values="log2FoldChange_rapid_vs_nonrapid").reindex(GENES)[ENDPOINTS]
nominal = ti.pivot(index="gene", columns="endpoint", values="nominal_p_lt_0_05").reindex(GENES)[ENDPOINTS]
found = ~heat.isna()

# Source data exports for the main figure.
ti.to_csv(OUTDIR / "45_main_figure4_panelAB_source_data.csv", index=False, encoding="utf-8-sig")
enrich_main = enrich[(enrich["signature"] == "candidate_10_gene") & (enrich["compartment"] == "tubulointerstitium")].copy()
# Benjamini-Hochberg correction across the candidate-gene tubulointerstitial enrichment family
# (5 progression endpoints x 2 nominal gene-set definitions).
pvals = enrich_main["fisher_p"].astype(float).to_numpy()
order = np.argsort(pvals)
qvals = np.empty_like(pvals, dtype=float)
prev = 1.0
m = len(pvals)
for rank_from_end, idx in enumerate(order[::-1], start=1):
    rank = m - rank_from_end + 1
    q = min(prev, pvals[idx] * m / rank)
    qvals[idx] = q
    prev = q
enrich_main["bh_q_displayed_family"] = qvals
enrich_main.to_csv(OUTDIR / "45_main_figure4_panelC_source_data_bhq.csv", index=False, encoding="utf-8-sig")
summary_main = summary[(summary["signature"] == "candidate_10_gene") & (summary["compartment"] == "tubulointerstitium")].copy()
summary_main.to_csv(OUTDIR / "45_main_figure4_panelD_source_data.csv", index=False, encoding="utf-8-sig")

# Helper functions.
def panel_label(ax, label):
    ax.text(-0.10, 1.05, label, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="top")


def p_label(p):
    if pd.isna(p):
        return "NA"
    p = float(p)
    return f"{p:.1e}" if p < 0.001 else f"{p:.3g}"

# Layout.
fig = plt.figure(figsize=(7.35, 6.85))
gs = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.08, 1.0], hspace=0.42, wspace=0.44)
axA = fig.add_subplot(gs[:, 0])
right = gs[:, 1].subgridspec(3, 1, height_ratios=[0.95, 0.85, 0.80], hspace=0.58)
axB = fig.add_subplot(right[0, 0])
axC = fig.add_subplot(right[1, 0])
axD = fig.add_subplot(right[2, 0])

# Panel A: heatmap across tubulointerstitial progression endpoints.
panel_label(axA, "A")
vals = heat.to_numpy(dtype=float)
masked = np.ma.masked_invalid(vals)
cmap = mpl.colormaps["RdBu_r"].copy()
cmap.set_bad(COL_MISS)
max_abs = np.nanmax(np.abs(vals))
im = axA.imshow(masked, cmap=cmap, vmin=-max_abs, vmax=max_abs, aspect="auto")
axA.set_xticks(np.arange(len(ENDPOINTS)))
axA.set_xticklabels([ENDPOINT_LABELS[e] for e in ENDPOINTS], rotation=35, ha="right")
axA.set_yticks(np.arange(len(GENES)))
axA.set_yticklabels(GENES)
axA.set_title("Candidate genes in Acoba 2025 tubulointerstitial progression-DE tables", loc="left", fontsize=8, fontweight="bold")
for i, gene in enumerate(GENES):
    for j, endpoint in enumerate(ENDPOINTS):
        v = heat.loc[gene, endpoint]
        if pd.isna(v):
            axA.text(j, i, "not\nfound", ha="center", va="center", fontsize=5.2, color="#555555")
        else:
            star = "*" if bool(nominal.loc[gene, endpoint]) else ""
            color = "white" if abs(float(v)) > max_abs*0.55 else "#222222"
            axA.text(j, i, f"{float(v):.2f}{star}", ha="center", va="center", fontsize=6.2, color=color)
for y in np.arange(-0.5, len(GENES), 1):
    axA.axhline(y, color="white", lw=0.7)
for x in np.arange(-0.5, len(ENDPOINTS), 1):
    axA.axvline(x, color="white", lw=0.7)
cbar = fig.colorbar(im, ax=axA, orientation="horizontal", fraction=0.035, pad=0.12)
cbar.set_label("log2FC in rapid progression", fontsize=6.5)
cbar.ax.tick_params(labelsize=6)
# Panel B: CKD stage advancement lollipop/dot for recovered genes.
panel_label(axB, "B")
ckd = ti[ti["endpoint"].astype(str) == "CKD_stage_advancement"].copy()
ckd["gene"] = pd.Categorical(ckd["gene"], categories=GENES, ordered=True)
ckd = ckd.sort_values("gene")
ypos = np.arange(len(GENES))[::-1]
for y, gene in zip(ypos, GENES):
    row = ckd[ckd["gene"].astype(str) == gene]
    if row.empty or pd.isna(row.iloc[0]["log2FoldChange_rapid_vs_nonrapid"]):
        axB.scatter(0, y, s=28, color=COL_MISS, edgecolor="white", linewidth=0.4, zorder=3)
        axB.text(0.05, y, "not found", va="center", ha="left", fontsize=5.7, color="#777777")
    else:
        v = float(row.iloc[0]["log2FoldChange_rapid_vs_nonrapid"])
        nom = bool(row.iloc[0]["nominal_p_lt_0_05"])
        axB.plot([0, v], [y, y], color="#CFCFCF", lw=1.0, zorder=1)
        axB.scatter(v, y, s=46, color=COL_UP if v >= 0 else COL_DOWN, edgecolor="white", linewidth=0.5, zorder=3)
        if nom:
            axB.text(v + 0.06, y, "*", va="center", ha="left", fontsize=9, color=COL_ACCENT)
axB.axvline(0, color="#BBBBBB", lw=0.8)
axB.set_yticks(ypos)
axB.set_yticklabels(GENES)
axB.set_xlabel("log2FC rapid vs non-rapid")
axB.set_title("CKD-stage advancement highlights injury/remodeling genes", loc="left", fontsize=8, fontweight="bold")
axB.grid(axis="x", color=COL_LIGHT, lw=0.6)
axB.set_xlim(-0.25, 1.55)

# Panel C: enrichment among progression-associated genes.
panel_label(axC, "C")
en_up = enrich_main[enrich_main["set_tested"] == "nominal_up_in_rapid"].copy()
en_up["endpoint"] = pd.Categorical(en_up["endpoint"], categories=ENDPOINTS, ordered=True)
en_up = en_up.sort_values("endpoint")
y = np.arange(len(en_up))[::-1]
colors = [COL_ACCENT if e == "CKD_stage_advancement" else "#B7C6D8" for e in en_up["endpoint"].astype(str)]
axC.barh(y, -np.log10(en_up["fisher_p"].astype(float).clip(lower=1e-300)), color=colors, edgecolor="white", height=0.65)
axC.axvline(-np.log10(0.05), color="#333333", lw=0.9, ls="--")
axC.set_yticks(y)
axC.set_yticklabels([ENDPOINT_LABELS[e].replace("\n", " ") for e in en_up["endpoint"].astype(str)])
axC.set_xlabel("-log10 Fisher P\n(overlap with nominal up genes)")
axC.set_title("Enrichment strongest for CKD-stage advancement", loc="left", fontsize=8, fontweight="bold")
for yi, (_, row) in zip(y, en_up.iterrows()):
    p = float(row["fisher_p"])
    orv = row["odds_ratio"]
    or_text = "inf" if np.isinf(orv) else f"{float(orv):.1f}"
    axC.text(-np.log10(max(p, 1e-300)) + 0.05, yi, f"OR {or_text}; P={p_label(p)}; q={p_label(row['bh_q_displayed_family'])}", va="center", ha="left", fontsize=5.6)
axC.grid(axis="x", color=COL_LIGHT, lw=0.6)

# Panel D: directionality among captured genes by endpoint.
panel_label(axD, "D")
s = summary_main.copy()
s["endpoint"] = pd.Categorical(s["endpoint"], categories=ENDPOINTS, ordered=True)
s = s.sort_values("endpoint")
y = np.arange(len(s))[::-1]
up = s["n_up_in_rapid"].astype(int).to_numpy()
down = -s["n_down_in_rapid"].astype(int).to_numpy()
axD.barh(y, up, color=COL_UP, edgecolor="white", height=0.65, label="Up in rapid")
axD.barh(y, down, color=COL_DOWN, edgecolor="white", height=0.65, label="Down in rapid")
axD.axvline(0, color="#333333", lw=0.8)
axD.set_yticks(y)
axD.set_yticklabels([ENDPOINT_LABELS[e].replace("\n", " ") for e in s["endpoint"].astype(str)])
axD.set_xlabel("Captured candidate genes")
axD.set_title("Most captured genes trend upward in rapid progression", loc="left", fontsize=8, fontweight="bold")
for yi, u, d in zip(y, up, -down):
    if u:
        axD.text(u + 0.12, yi, str(u), va="center", ha="left", fontsize=6.3)
    if d:
        axD.text(-d - 0.12, yi, str(d), va="center", ha="right", fontsize=6.3)
axD.legend(loc="lower right", fontsize=6)
axD.set_xlim(-6.2, 6.5)
axD.grid(axis="x", color=COL_LIGHT, lw=0.6)

fig.suptitle("Figure 4. External DKD progression data support the candidate at gene level", x=0.01, y=0.995, ha="left", fontsize=10, fontweight="bold")
fig.text(0.01, 0.010, "Acoba 2025 progression-DE tables were analysed at gene level. Enrichment q values use Benjamini-Hochberg correction across the candidate-gene tubulointerstitial enrichment family.", ha="left", va="bottom", fontsize=6.2, color="#555555")
fig.subplots_adjust(top=0.93, bottom=0.14, left=0.08, right=0.98)

for ext in ["svg", "pdf", "png", "tiff"]:
    if ext in ["png", "tiff"]:
        fig.savefig(f"{PREFIX}.{ext}", dpi=600, bbox_inches="tight")
    else:
        fig.savefig(f"{PREFIX}.{ext}", bbox_inches="tight")
plt.close(fig)

# Legend
ckd_en = en_up[en_up["endpoint"].astype(str) == "CKD_stage_advancement"].iloc[0]
legend = f"""# Main Figure 4. External DKD progression-context gene-level support

## Core figure claim

External Acoba 2025 DKD progression differential-expression tables support parts of the 10-gene injury-dominant remodeling candidate at gene level, especially in tubulointerstitial CKD-stage advancement. This is progression-context alignment, not individual-level prognostic validation.

## Panel legend draft

**Figure 4 | External DKD progression-associated differential-expression tables support the candidate at gene level.**
**A,** Candidate-gene log2 fold changes in Acoba 2025 tubulointerstitial DKD progression differential-expression tables. Asterisks mark nominal P<0.05; grey cells indicate genes not recovered from the parsed endpoint table.
**B,** Gene-level effects for the tubulointerstitial CKD-stage-advancement endpoint. COL1A1, HAVCR1 and FN1 were nominally upregulated in rapid progressors, while LCN2, TGFB1, TLR4 and NLRP3 were not recovered from the parsed endpoint table.
**C,** Fisher enrichment analysis testing whether recovered candidate genes overlapped nominally upregulated rapid-progression genes. Enrichment was strongest for tubulointerstitial CKD-stage advancement (odds ratio={float(ckd_en['odds_ratio']):.2f}; Fisher P={float(ckd_en['fisher_p']):.4g}; BH q={float(ckd_en['bh_q_displayed_family']):.4g} across the candidate-gene tubulointerstitial enrichment family).
**D,** Directionality among candidate genes captured in each tubulointerstitial endpoint table. CKD-stage advancement and composite outcome showed all captured candidate genes higher in rapid progressors.

## Claim boundary

This figure supports external DKD progression-context alignment at gene level, with the strongest enrichment nominally significant but attenuated after BH correction. It does not establish an individual-level prognostic score.

## Files

- `45_main_figure4_acoba_progression_context_bhq.svg`
- `45_main_figure4_acoba_progression_context_bhq.pdf`
- `45_main_figure4_acoba_progression_context_bhq.png`
- `45_main_figure4_acoba_progression_context_bhq.tiff`
- `45_main_figure4_panelAB_source_data.csv`
- `45_main_figure4_panelC_source_data_bhq.csv`
- `45_main_figure4_panelD_source_data.csv`
"""
(OUTDIR / "45_main_figure4_acoba_progression_context_bhq_legend.md").write_text(legend, encoding="utf-8")

print("Figure 4 assembly complete")
print(PREFIX.with_suffix(".png"))
print(OUTDIR / "45_main_figure4_acoba_progression_context_bhq_legend.md")



# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
