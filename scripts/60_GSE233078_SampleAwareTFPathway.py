#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    }
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
IN_DIR = PROJECT_ROOT / "results" / "public_data_gse233078_spp1_cd44_deg_enrichment"
OUT_DIR = PROJECT_ROOT / "results" / "gse233078_sample_aware_tf_pathway"


PROGRAMS = {
    "DTL_HIF_STAT_AP1_injury": {
        "focus_group": "DTL",
        "theme": "tubular hypoxia/inflammatory stress",
        "genes": ["Hif1a", "Stat3", "Jun", "Fos", "Jund", "Atf3", "Egr1", "Spp1", "Havcr1", "Lcn2", "Vcam1", "Clu", "Axl"],
    },
    "DTL_SPP1_CD44_damage_context": {
        "focus_group": "DTL",
        "theme": "tubular Spp1/Cd44-associated damage context",
        "genes": ["Spp1", "Cd44", "Havcr1", "Lcn2", "Vcam1", "Clu", "Axl", "C3", "Ccl2", "Cxcl10", "Tgfb1"],
    },
    "MonoMacro_NFKB_AP1_CD44_inflammation": {
        "focus_group": "Mono_Macro",
        "theme": "myeloid NF-kB/AP-1 inflammatory remodeling",
        "genes": ["Nfkb1", "Rela", "Relb", "Jun", "Fos", "Jund", "Stat3", "Cd44", "Tnf", "Il1b", "Ccl2", "Ccl3", "Ccl4", "Nlrp3", "Tlr4"],
    },
    "MonoMacro_TREM2_APOE_phagocytic_state": {
        "focus_group": "Mono_Macro",
        "theme": "Trem2/Apoe phagocytic macrophage-like state",
        "genes": ["Trem2", "Apoe", "Tyrobp", "Aif1", "Lyz2", "Csf1r", "Itgam", "C1qa", "C1qb", "C1qc", "Lpl", "Spp1", "Cd44"],
    },
    "MonoMacro_SPP1_CD44_fibroinflammatory_context": {
        "focus_group": "Mono_Macro",
        "theme": "Spp1/Cd44 fibro-inflammatory context",
        "genes": ["Spp1", "Cd44", "Tgfb1", "Fn1", "Col1a1", "Ccl2", "Cxcl10", "Nlrp3", "Tlr4", "Mmp9", "Mmp12"],
    },
}

CONTRASTS = [
    ("obese_vs_lean", "lean", "obese", "disease_induction"),
    ("enalapril_vs_obese", "obese", "obese_enalapril", "treatment_reversal"),
]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return 0.0
    pooled = (((len(a) - 1) * np.var(a, ddof=1)) + ((len(b) - 1) * np.var(b, ddof=1))) / max(len(a) + len(b) - 2, 1)
    if not np.isfinite(pooled) or pooled <= 0:
        return 0.0
    return float((np.mean(b) - np.mean(a)) / math.sqrt(pooled))


def welch_p(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return 1.0
    p = ttest_ind(b, a, equal_var=False, nan_policy="omit").pvalue
    return float(p) if np.isfinite(p) else 1.0


def mann_p(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 1 or len(b) < 1:
        return 1.0
    try:
        p = mannwhitneyu(b, a, alternative="two-sided").pvalue
    except ValueError:
        p = 1.0
    return float(p) if np.isfinite(p) else 1.0


def load_logcpm() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(IN_DIR / "pseudobulk_counts_focus_groups.csv")
    meta = pd.read_csv(IN_DIR / "pseudobulk_group_metadata.csv")
    counts = counts.set_index("gene")
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib = meta.set_index("group_id")["library_size"].astype(float)
    common = [col for col in counts.columns if col in lib.index]
    counts = counts[common]
    lib = lib.loc[common]
    logcpm = np.log2(counts.div(lib, axis=1) * 1_000_000 + 1.0)
    return logcpm, meta


def score_programs(logcpm: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    presence_rows = []
    for program, spec in PROGRAMS.items():
        focus_group = spec["focus_group"]
        requested = spec["genes"]
        present = [gene for gene in requested if gene in logcpm.index]
        presence_rows.append(
            {
                "program": program,
                "focus_group": focus_group,
                "theme": spec["theme"],
                "n_requested": len(requested),
                "n_present": len(present),
                "present_genes": ";".join(present),
                "missing_genes": ";".join([gene for gene in requested if gene not in logcpm.index]),
            }
        )
        if not present:
            continue
        meta_sub = meta[meta["focus_group"] == focus_group].copy()
        group_ids = [gid for gid in meta_sub["group_id"] if gid in logcpm.columns]
        mat = logcpm.loc[present, group_ids]
        z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
        scores = z.mean(axis=0)
        for group_id, score in scores.items():
            m = meta_sub[meta_sub["group_id"] == group_id].iloc[0]
            rows.append(
                {
                    "program": program,
                    "theme": spec["theme"],
                    "focus_group": focus_group,
                    "group_id": group_id,
                    "sample_id": m["sample_id"],
                    "condition": m["condition"],
                    "n_cells": int(m["n_cells"]),
                    "library_size": float(m["library_size"]),
                    "score": float(score),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(presence_rows)


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for program, sub in scores.groupby("program", sort=False):
        for contrast, ref_group, test_group, role in CONTRASTS:
            ref = sub[sub["condition"] == ref_group]["score"].to_numpy(dtype=float)
            test = sub[sub["condition"] == test_group]["score"].to_numpy(dtype=float)
            delta = float(np.mean(test) - np.mean(ref)) if len(ref) and len(test) else 0.0
            rows.append(
                {
                    "program": program,
                    "theme": sub["theme"].iloc[0],
                    "focus_group": sub["focus_group"].iloc[0],
                    "contrast": contrast,
                    "contrast_role": role,
                    "reference_group": ref_group,
                    "test_group": test_group,
                    "n_reference": len(ref),
                    "n_test": len(test),
                    "reference_mean": float(np.mean(ref)) if len(ref) else 0.0,
                    "test_mean": float(np.mean(test)) if len(test) else 0.0,
                    "delta": delta,
                    "cohen_d": cohen_d(ref, test),
                    "welch_p": welch_p(ref, test),
                    "mannwhitney_p": mann_p(ref, test),
                    "direction": "up" if delta > 0.15 else "down" if delta < -0.15 else "flat",
                }
            )
    out = pd.DataFrame(rows)
    out["bh_fdr_welch"] = multipletests(out["welch_p"].fillna(1.0), method="fdr_bh")[1]
    out["bh_fdr_mannwhitney"] = multipletests(out["mannwhitney_p"].fillna(1.0), method="fdr_bh")[1]
    return out


def save_boxplot(scores: pd.DataFrame, out_path: Path) -> None:
    programs = list(PROGRAMS.keys())
    fig, axes = plt.subplots(len(programs), 1, figsize=(8.5, 11), sharex=True)
    conds = ["lean", "obese", "obese_enalapril"]
    colors = ["#d8e8d2", "#d9a6a6", "#a8c7e6"]
    for ax, program in zip(axes, programs):
        sub = scores[scores["program"] == program]
        values = [sub[sub["condition"] == cond]["score"].to_numpy(dtype=float) for cond in conds]
        ax.boxplot(values, tick_labels=conds, patch_artist=True, medianprops={"color": "black"})
        for patch, color in zip(ax.patches, colors):
            patch.set_facecolor(color)
        for i, vals in enumerate(values, start=1):
            ax.scatter(np.full(len(vals), i), vals, s=34, color="#2f4858", alpha=0.8, zorder=3)
        ax.axhline(0, lw=0.8, color="#999999", alpha=0.7)
        ax.set_ylabel("z-score")
        ax.set_title(f"{program} ({PROGRAMS[program]['focus_group']})", fontsize=10)
    fig.suptitle("Sample-aware TF/pathway activity in GSE233078 focus cell groups", fontsize=14, y=0.995)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_delta_heatmap(summary: pd.DataFrame, out_path: Path) -> None:
    mat = summary.pivot_table(index="program", columns="contrast", values="delta", aggfunc="first").fillna(0.0)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    vmax = max(0.5, float(np.max(np.abs(mat.values))))
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=25, ha="right", fontsize=9)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels([idx.replace("_", " ") for idx in mat.index], fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("TF/pathway program score shifts")
    fig.colorbar(im, ax=ax, label="Mean score delta")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logcpm, meta = load_logcpm()
    scores, presence = score_programs(logcpm, meta)
    summary = summarize(scores)
    scores.to_csv(OUT_DIR / "sample_aware_tf_pathway_scores.csv", index=False, encoding="utf-8-sig")
    presence.to_csv(OUT_DIR / "sample_aware_tf_pathway_gene_presence.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT_DIR / "sample_aware_tf_pathway_contrast_summary.csv", index=False, encoding="utf-8-sig")
    save_boxplot(scores, OUT_DIR / "FigureBR_sample_aware_tf_pathway_boxplots.png")
    save_delta_heatmap(summary, OUT_DIR / "FigureBS_sample_aware_tf_pathway_delta_heatmap.png")

    lines = [
        "GSE233078 sample-aware TF/pathway activity supplement",
        "",
        "Purpose:",
        "This layer complements the previous regulator atlas by scoring curated TF/pathway programs at the sample-pseudobulk level in DTL and Mono/Macro.",
        "Scores are based on logCPM gene-level z-scores within each focus group, then averaged per program.",
        "",
        "Key contrast readouts:",
    ]
    for _, row in summary.sort_values(["program", "contrast"]).iterrows():
        lines.append(
            f"- {row['program']} | {row['contrast']}: delta={row['delta']:.3f}, "
            f"d={row['cohen_d']:.3f}, Welch p={row['welch_p']:.4g}, MW p={row['mannwhitney_p']:.4g}, direction={row['direction']}"
        )
    lines.extend(
        [
            "",
            "Writing boundary:",
            "1. This is curated TF/pathway activity scoring, not full pySCENIC regulon inference.",
            "2. Because sample size is small, this should be written as mechanistic support rather than definitive regulatory proof.",
            "3. Programs showing obese-up and enalapril-down behavior are suitable for Figure 4 support or Supplementary Figure.",
        ]
    )
    write_text(OUT_DIR / "status.txt", "completed\n")


if __name__ == "__main__":
    main()


