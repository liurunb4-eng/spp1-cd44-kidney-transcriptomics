# -*- coding: utf-8 -*-
"""
Step 12: clinical eGFR stratification audit for the relaunched injury-dominant DKD/TI signature.

Purpose:
- Recompute GSE175759 clinical cross-sectional eGFR association for the new primary candidate,
  not the old SPP1/CD44-centered program.
- Explicitly audit whether true longitudinal/prognostic analysis is possible from local metadata.
"""
from pathlib import Path
import os
import gzip
import re
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

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

REPO_ROOT = Path(__file__).resolve().parents[3]
BASE = Path(os.environ.get("AJP_DATA_ROOT", REPO_ROOT))
OUT = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUT.mkdir(parents=True, exist_ok=True)
SCORES = OUT / "06_frozen_signature_scores_gse175759.csv.gz"
META = BASE / "results" / "human_bulk_gse175759_clinical_correlation" / "gse175759_metadata_clean.csv"
CLIN_SCREEN = BASE / "data" / "public_datasets" / "Human_CKD_bulk_clinical_screen"

PRIMARY = "candidate_DKD_TI_remodeling_no_SPP1_CD44"
KEEP_SIGNATURES = [
    PRIMARY,
    "human_tubular_injury_without_SPP1_CD44",
    "human_tubular_injury_context",
    "SPP1_CD44_anchor",
    "old_NDT_composite",
]
LABELS = {
    PRIMARY: "10-gene injury-dominant\nDKD/TI candidate",
    "human_tubular_injury_without_SPP1_CD44": "Human tubular injury\nwithout SPP1/CD44",
    "human_tubular_injury_context": "Human tubular injury\ncontext",
    "SPP1_CD44_anchor": "SPP1/CD44\nanchor",
    "old_NDT_composite": "Old NDT\ncomposite",
}


def bh_fdr(pvals):
    p = np.asarray([np.nan if x is None else x for x in pvals], dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        return out
    idx = np.where(mask)[0]
    order = idx[np.argsort(p[mask])]
    ranked = p[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out[order] = adj
    return out


def safe_spearman(x, y):
    df = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(df) < 4 or df["x"].nunique() < 2 or df["y"].nunique() < 2:
        return len(df), np.nan, np.nan
    r, p = stats.spearmanr(df["x"], df["y"])
    return len(df), float(r), float(p)


def mannwhitney_delta(a, b):
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return len(a), len(b), np.nan, np.nan
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    return len(a), len(b), float(a.mean() - b.mean()), float(p)


def scan_longitudinal_fields():
    keywords = re.compile(r"(follow.?up|followup|esrd|end.?stage|dialysis|transplant|survival|death|event|outcome|endpoint|slope|progression|time.?to|doubling|longitudinal)", re.I)
    hits = []
    for p in CLIN_SCREEN.rglob("*"):
        if not p.is_file():
            continue
        name = p.name.lower()
        if not (name.endswith((".csv", ".tsv", ".txt")) or name.endswith((".txt.gz", ".csv.gz", ".tsv.gz"))):
            continue
        if p.stat().st_size > 50_000_000:
            continue
        rel = str(p.relative_to(CLIN_SCREEN))
        try:
            if name.endswith(".gz"):
                with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                    lines = []
                    for _ in range(300):
                        try:
                            lines.append(next(fh))
                        except StopIteration:
                            break
                matched = [ln.strip()[:220] for ln in lines if keywords.search(ln)]
            else:
                try:
                    df0 = pd.read_csv(p, nrows=5)
                except Exception:
                    df0 = pd.read_csv(p, sep="\t", nrows=5)
                matched = [str(c) for c in df0.columns if keywords.search(str(c))]
                rowtxt = df0.head(5).to_string()
                if keywords.search(rowtxt):
                    matched.append("row_text_contains_longitudinal_keyword")
            if matched:
                hits.append({"file": rel, "matched_terms_or_lines": " | ".join(matched[:8])})
        except Exception:
            pass
    return pd.DataFrame(hits)

# Load and merge
scores = pd.read_csv(SCORES)
meta = pd.read_csv(META)
meta = meta.rename(columns={"gsm": "sample_id"})
merged = scores.merge(meta, on="sample_id", how="left")
merged = merged[merged["signature"].isin(KEEP_SIGNATURES)].copy()
merged["usable"] = merged["egfr_ckd_epi"].notna() & (merged["technical_outlier"].astype(str).str.lower() == "no")
merged["disease_group"] = np.where(merged["diagnosis"].eq("Control"), "Control", "CKD")

usable_meta = meta[(meta["egfr_ckd_epi"].notna()) & (meta["technical_outlier"].astype(str).str.lower()=="no")].copy()
# tertiles from all usable samples, matching earlier logic
usable_meta["egfr_tertile"] = pd.qcut(usable_meta["egfr_ckd_epi"], 3, labels=["Low eGFR", "Mid eGFR", "High eGFR"])
tertile_map = usable_meta.set_index("sample_id")["egfr_tertile"].to_dict()
merged["egfr_tertile"] = merged["sample_id"].map(tertile_map)
merged["egfr_60_group"] = np.where(merged["egfr_ckd_epi"] < 60, "eGFR <60", "eGFR >=60")
merged.loc[~merged["usable"], "egfr_60_group"] = np.nan
merged["signature_label"] = merged["signature"].map(LABELS).fillna(merged["signature"])

# Correlation analyses by signature and subgroup
subgroups = {
    "all_usable": lambda d: d["usable"],
    "ckd_only": lambda d: d["usable"] & d["disease_group"].eq("CKD"),
    "igan_only": lambda d: d["usable"] & d["diagnosis"].eq("IgAN"),
    "non_control_non_igan": lambda d: d["usable"] & (~d["diagnosis"].isin(["Control", "IgAN"])),
    "diabetic_nephropathy_only_descriptive": lambda d: d["usable"] & d["diagnosis"].eq("Diabetic nephropathy"),
}
cor_rows = []
for sig, g in merged.groupby("signature"):
    for subgroup, fn in subgroups.items():
        dd = g.loc[fn(g)].copy()
        n, rho, p = safe_spearman(dd["score"], dd["egfr_ckd_epi"])
        cor_rows.append({
            "signature": sig,
            "signature_label": LABELS.get(sig, sig),
            "subgroup": subgroup,
            "n": n,
            "spearman_rho": rho,
            "spearman_p": p,
        })
cor = pd.DataFrame(cor_rows)
cor["spearman_fdr_within_subgroup"] = np.nan
for subgroup in cor["subgroup"].unique():
    idx = cor["subgroup"].eq(subgroup)
    cor.loc[idx, "spearman_fdr_within_subgroup"] = bh_fdr(cor.loc[idx, "spearman_p"].values)

# Low vs high eGFR tertile and eGFR<60 sensitivity
tert_rows = []
for sig, g in merged[merged["usable"]].groupby("signature"):
    low = g.loc[g["egfr_tertile"].eq("Low eGFR"), "score"]
    high = g.loc[g["egfr_tertile"].eq("High eGFR"), "score"]
    n_low, n_high, delta, p = mannwhitney_delta(low, high)
    tert_rows.append({
        "signature": sig, "signature_label": LABELS.get(sig, sig),
        "contrast": "low_vs_high_egfr_tertile", "n_low_or_case": n_low, "n_high_or_ref": n_high,
        "mean_delta": delta, "mannwhitney_p": p,
        "interpretation": "higher_in_low_eGFR" if np.isfinite(delta) and delta > 0 else "lower_or_no_increase_in_low_eGFR"
    })
    low60 = g.loc[g["egfr_60_group"].eq("eGFR <60"), "score"]
    hi60 = g.loc[g["egfr_60_group"].eq("eGFR >=60"), "score"]
    n_case, n_ref, delta60, p60 = mannwhitney_delta(low60, hi60)
    tert_rows.append({
        "signature": sig, "signature_label": LABELS.get(sig, sig),
        "contrast": "egfr_lt60_vs_ge60_descriptive", "n_low_or_case": n_case, "n_high_or_ref": n_ref,
        "mean_delta": delta60, "mannwhitney_p": p60,
        "interpretation": "higher_in_eGFR_lt60" if np.isfinite(delta60) and delta60 > 0 else "lower_or_no_increase_in_eGFR_lt60"
    })
tert = pd.DataFrame(tert_rows)
tert["mannwhitney_fdr_within_contrast"] = np.nan
for contrast in tert["contrast"].unique():
    idx = tert["contrast"].eq(contrast)
    tert.loc[idx, "mannwhitney_fdr_within_contrast"] = bh_fdr(tert.loc[idx, "mannwhitney_p"].values)

# Dataset inventory / longitudinal feasibility
long_hits = scan_longitudinal_fields()
# diagnosis composition
composition = usable_meta["diagnosis"].value_counts().rename_axis("diagnosis").reset_index(name="n")

# Source data for primary plot
primary_df = merged[(merged["signature"].eq(PRIMARY)) & (merged["usable"])].copy()
primary_df["egfr_tertile"] = pd.Categorical(primary_df["egfr_tertile"], ["Low eGFR", "Mid eGFR", "High eGFR"], ordered=True)
primary_df = primary_df.sort_values(["egfr_tertile", "diagnosis", "sample_id"])

# Write tables
prefix = OUT / "12_clinical_egfr_stratification"
merged.to_csv(OUT / "12_clinical_egfr_stratification_scores_merged.csv", index=False, encoding="utf-8-sig")
cor.to_csv(OUT / "12_clinical_egfr_stratification_correlations.csv", index=False, encoding="utf-8-sig")
tert.to_csv(OUT / "12_clinical_egfr_stratification_tests.csv", index=False, encoding="utf-8-sig")
composition.to_csv(OUT / "12_clinical_egfr_stratification_diagnosis_composition.csv", index=False, encoding="utf-8-sig")
long_hits.to_csv(OUT / "12_longitudinal_field_screen.csv", index=False, encoding="utf-8-sig")
primary_df.to_csv(OUT / "12_clinical_egfr_stratification_source_data.csv", index=False, encoding="utf-8-sig")

# Figure: quantitative grid
primary_cor_all = cor[(cor.signature==PRIMARY) & (cor.subgroup=="all_usable")].iloc[0]
primary_cor_ckd = cor[(cor.signature==PRIMARY) & (cor.subgroup=="ckd_only")].iloc[0]
primary_test_tert = tert[(tert.signature==PRIMARY) & (tert.contrast=="low_vs_high_egfr_tertile")].iloc[0]

fig = plt.figure(figsize=(8.2, 4.9), constrained_layout=False)
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.25], height_ratios=[1, 1], wspace=0.42, hspace=0.55)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[1,0])
ax4 = fig.add_subplot(gs[1,1])

palette = {"Low eGFR":"#c75d5d", "Mid eGFR":"#c9b36a", "High eGFR":"#5b8fc1"}
order=["Low eGFR", "Mid eGFR", "High eGFR"]
# violin + points
vals=[primary_df.loc[primary_df.egfr_tertile.eq(o), "score"].values for o in order]
parts=ax1.violinplot(vals, positions=np.arange(1,4), widths=0.75, showmeans=False, showmedians=True, showextrema=False)
for i, body in enumerate(parts['bodies']):
    body.set_facecolor(palette[order[i]])
    body.set_alpha(0.45)
    body.set_edgecolor('none')
parts['cmedians'].set_color('#333333')
parts['cmedians'].set_linewidth(1.1)
rng=np.random.default_rng(20261005)
for i,o in enumerate(order, start=1):
    y=primary_df.loc[primary_df.egfr_tertile.eq(o), "score"].values
    x=i+rng.normal(0,0.045,size=len(y))
    ax1.scatter(x,y,s=12,color=palette[o],edgecolor='white',linewidth=0.35,alpha=0.9,zorder=3)
ax1.set_xticks([1,2,3], ["Low", "Mid", "High"])
ax1.set_xlabel("eGFR tertile")
ax1.set_ylabel("10-gene candidate score")
ax1.set_title("A  eGFR tertiles", loc='left', fontweight='bold', fontsize=8.5)
ax1.text(0.02, 0.98, f"Low vs high: Δ={primary_test_tert.mean_delta:.2f}\nFDR={primary_test_tert.mannwhitney_fdr_within_contrast:.3g}", transform=ax1.transAxes, va='top', ha='left', fontsize=7)

# scatter
colors = primary_df["disease_group"].map({"CKD":"#6d6d6d", "Control":"#9bc3dc"})
ax2.scatter(primary_df["egfr_ckd_epi"], primary_df["score"], s=18, c=colors, edgecolor='white', linewidth=0.35, alpha=0.9)
# line fit for display only
x=primary_df["egfr_ckd_epi"].values; y=primary_df["score"].values
m,b=np.polyfit(x,y,1)
xx=np.linspace(np.nanmin(x),np.nanmax(x),100)
ax2.plot(xx,m*xx+b,color="#333333",lw=1.0)
ax2.set_xlabel("eGFR (CKD-EPI)")
ax2.set_ylabel("10-gene candidate score")
ax2.set_title("B  Continuous eGFR association", loc='left', fontweight='bold', fontsize=8.5)
ax2.text(0.02,0.98, f"All: ρ={primary_cor_all.spearman_rho:.2f}, FDR={primary_cor_all.spearman_fdr_within_subgroup:.3g}\nCKD only: ρ={primary_cor_ckd.spearman_rho:.2f}, FDR={primary_cor_ckd.spearman_fdr_within_subgroup:.3g}", transform=ax2.transAxes, va='top', ha='left', fontsize=7)
ax2.scatter([],[],s=18,c="#6d6d6d",label="CKD")
ax2.scatter([],[],s=18,c="#9bc3dc",label="Control")
ax2.legend(loc='lower left', borderaxespad=0.2, handletextpad=0.3, fontsize=7)

# sensitivity rho horizontal bar for primary
sens = cor[cor.signature.eq(PRIMARY)].copy()
sens_order=["all_usable","ckd_only","igan_only","non_control_non_igan","diabetic_nephropathy_only_descriptive"]
sens=sens.set_index('subgroup').loc[sens_order].reset_index()
sens_labels=["All (n=%d)"%sens.loc[0,'n'], "CKD only (n=%d)"%sens.loc[1,'n'], "IgAN only (n=%d)"%sens.loc[2,'n'], "non-IgAN CKD (n=%d)"%sens.loc[3,'n'], "DN only (n=%d)"%sens.loc[4,'n']]
ypos=np.arange(len(sens))[::-1]
barcolors=["#5b8fc1" if np.isfinite(r) and r<0 else "#d6d6d6" for r in sens.spearman_rho]
ax3.barh(ypos, sens.spearman_rho.fillna(0), color=barcolors, height=0.65)
ax3.axvline(0,color="#333333",lw=0.7)
ax3.set_yticks(ypos, sens_labels)
ax3.set_xlabel("Spearman ρ with eGFR")
ax3.set_xlim(-0.55,0.18)
ax3.set_title("C  Subgroup boundary", loc='left', fontweight='bold', fontsize=8.5)
for yv,row in zip(ypos, sens.itertuples()):
    if np.isfinite(row.spearman_fdr_within_subgroup):
        ax3.text(row.spearman_rho-0.02 if row.spearman_rho<0 else row.spearman_rho+0.02, yv, f"FDR {row.spearman_fdr_within_subgroup:.2g}", va='center', ha='right' if row.spearman_rho<0 else 'left', fontsize=6.5)
    else:
        ax3.text(0.03, yv, "n too small", va='center', ha='left', fontsize=6.5)

# comparison across selected signatures all usable
comp = cor[(cor.subgroup.eq("all_usable")) & (cor.signature.isin(KEEP_SIGNATURES))].copy()
comp = comp.set_index('signature').loc[KEEP_SIGNATURES].reset_index()
comp_labels=["10-gene candidate", "Tubular injury no SPP1/CD44", "Tubular injury context", "SPP1/CD44 anchor", "Old NDT composite"]
ypos=np.arange(len(comp))[::-1]
ax4.barh(ypos, comp.spearman_rho, color=["#c75d5d" if s==PRIMARY else "#b8c6d8" for s in comp.signature], height=0.65)
ax4.axvline(0,color="#333333",lw=0.7)
ax4.set_yticks(ypos, comp_labels)
ax4.set_xlabel("Spearman ρ with eGFR")
ax4.set_xlim(-0.38,0.05)
ax4.set_title("D  Context among retained scores", loc='left', fontweight='bold', fontsize=8.5)
for yv,row in zip(ypos, comp.itertuples()):
    ax4.text(row.spearman_rho-0.012 if row.spearman_rho<0 else row.spearman_rho+0.01, yv, f"{row.spearman_rho:.2f}", va='center', ha='right' if row.spearman_rho<0 else 'left', fontsize=6.5)

fig.suptitle("Cross-sectional clinical relevance of the injury-dominant DKD/TI candidate in GSE175759", x=0.02, y=0.995, ha='left', fontsize=9, fontweight='bold')
fig.text(0.02, 0.01, "GSE175759 provides baseline eGFR but no local longitudinal endpoint; do not write this as survival or progression prediction.", ha='left', va='bottom', fontsize=6.5)
fig.subplots_adjust(top=0.90, bottom=0.13, left=0.13, right=0.98)

fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
plt.close(fig)
# Markdown report
primary_lines = []
primary_lines.append("# Step 12. Clinical eGFR stratification audit for the relaunched signature\n")
primary_lines.append("Date: 2026-10-05\n")
primary_lines.append("## Purpose\n")
primary_lines.append("Recalculate the human clinical-extension analysis for the current primary 10-gene injury-dominant DKD/TI candidate, because the older GSE175759 clinical figure was built around the previous SPP1/CD44-centered program.\n")
primary_lines.append("## Figure contract\n")
primary_lines.append("Core conclusion: GSE175759 supports a cross-sectional association between the new injury-dominant candidate and lower eGFR, but it does not provide a longitudinal or survival endpoint.\n")
primary_lines.append("Archetype: quantitative grid.\n")
primary_lines.append("Evidence chain: eGFR tertiles, continuous eGFR correlation, subgroup boundary analysis, and comparison with retained scores.\n")
primary_lines.append("## Dataset boundary\n")
primary_lines.append(f"Usable samples with baseline eGFR and no technical outlier: n={len(usable_meta)}.\n")
primary_lines.append("Diagnosis composition:\n\n")
for row in composition.itertuples(index=False):
    primary_lines.append(f"- {row.diagnosis}: n={row.n}\n")
primary_lines.append("\nNo local field screening hit supported ESRD, dialysis, transplant, survival, follow-up duration, eGFR slope, time-to-event, or other longitudinal outcomes in the current Human_CKD_bulk_clinical_screen files. Therefore Cox/Kaplan-Meier/eGFR-slope analysis is not justified from the currently available local files.\n")
primary_lines.append("\n## Primary 10-gene candidate results\n")
primary_lines.append(f"- Low-vs-high eGFR tertile mean-score delta = {primary_test_tert.mean_delta:.3f}; Mann-Whitney p = {primary_test_tert.mannwhitney_p:.4g}; FDR within contrast = {primary_test_tert.mannwhitney_fdr_within_contrast:.4g}.\n")
primary_lines.append(f"- Continuous eGFR association in all usable samples: Spearman rho = {primary_cor_all.spearman_rho:.3f}; p = {primary_cor_all.spearman_p:.4g}; FDR = {primary_cor_all.spearman_fdr_within_subgroup:.4g}.\n")
primary_lines.append(f"- CKD-only sensitivity: Spearman rho = {primary_cor_ckd.spearman_rho:.3f}; p = {primary_cor_ckd.spearman_p:.4g}; FDR = {primary_cor_ckd.spearman_fdr_within_subgroup:.4g}.\n")
for subgroup in ["igan_only","non_control_non_igan","diabetic_nephropathy_only_descriptive"]:
    row=cor[(cor.signature==PRIMARY)&(cor.subgroup==subgroup)].iloc[0]
    primary_lines.append(f"- {subgroup}: n = {int(row.n)}; rho = {row.spearman_rho if np.isfinite(row.spearman_rho) else 'NA'}; FDR = {row.spearman_fdr_within_subgroup if np.isfinite(row.spearman_fdr_within_subgroup) else 'NA'}.\n")
primary_lines.append("\n## Manuscript-safe interpretation\n")
primary_lines.append("This analysis supports clinical relevance as a cross-sectional eGFR association in a heterogeneous human tubulointerstitial CKD cohort. It should be written as `associated with lower eGFR` or `marks reduced kidney function`, not as `predicts progression`, `drives eGFR decline`, or `stratifies ESRD risk`.\n")
primary_lines.append("\nThe diabetic nephropathy subset in GSE175759 is too small for disease-specific inference. The result remains a pan-CKD human clinical-extension layer that complements the DKD compartment evidence from GSE30122 and the single-cell/snRNA validation.\n")
primary_lines.append("\n## Outputs\n")
for fn in [
    "12_clinical_egfr_stratification.svg",
    "12_clinical_egfr_stratification.pdf",
    "12_clinical_egfr_stratification.png",
    "12_clinical_egfr_stratification.tiff",
    "12_clinical_egfr_stratification_source_data.csv",
    "12_clinical_egfr_stratification_correlations.csv",
    "12_clinical_egfr_stratification_tests.csv",
    "12_longitudinal_field_screen.csv",
]:
    primary_lines.append(f"- `{fn}`\n")

(OUT / "12_clinical_egfr_stratification.md").write_text("".join(primary_lines), encoding="utf-8")

print("DONE")
print((OUT / "12_clinical_egfr_stratification.md"))
print(primary_test_tert.to_dict())
print(primary_cor_all.to_dict())
print(primary_cor_ckd.to_dict())
print('longitudinal_hits', len(long_hits))

# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
