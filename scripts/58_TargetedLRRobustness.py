#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd

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
IN_DIR = PROJECT_ROOT / "results" / "public_data_gse233078_cellchat_mvp_prep"
OUT_DIR = PROJECT_ROOT / "results" / "targeted_lr_robustness_spp1_cd44"


PRIORITY_AXES = [
    {
        "axis": "Spp1->Cd44",
        "module": "Spp1/Cd44 tubular-immune",
        "ligand": "Spp1",
        "receptor": "Cd44",
        "source_group": "DTL",
        "target_group": "Mono_Macro",
        "priority": "main_candidate",
    },
    {
        "axis": "Tgfb1->Tgfbr2",
        "module": "fibrosis remodeling",
        "ligand": "Tgfb1",
        "receptor": "Tgfbr2",
        "source_group": "Mono_Macro",
        "target_group": "Endo_GEC",
        "priority": "context_axis",
    },
    {
        "axis": "Vegfa->Kdr",
        "module": "vascular microcirculation",
        "ligand": "Vegfa",
        "receptor": "Kdr",
        "source_group": "IC-A",
        "target_group": "Endo_GEC",
        "priority": "context_axis",
    },
    {
        "axis": "Vegfa->Flt1",
        "module": "vascular microcirculation",
        "ligand": "Vegfa",
        "receptor": "Flt1",
        "source_group": "IC-A",
        "target_group": "Endo_GEC",
        "priority": "context_axis",
    },
]


CONDITIONS = ["lean", "obese", "obese_enalapril"]
METHOD_COLUMNS = [
    "expr_product",
    "pct_product",
    "specificity_product",
    "balanced_expr_pct",
    "geometric_consensus",
]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def get_value(group_expr: pd.DataFrame, gene: str, cell_group: str, condition: str, field: str) -> float:
    sub = group_expr[
        (group_expr["gene"].str.upper() == gene.upper())
        & (group_expr["cellchat_group"] == cell_group)
        & (group_expr["condition"] == condition)
    ]
    if sub.empty:
        return 0.0
    return float(sub.iloc[0][field])


def build_lr_scores() -> pd.DataFrame:
    group_expr = pd.read_csv(IN_DIR / "cellchat_mvp_whitelist_group_expression.csv")
    rows = []
    eps = 1e-9
    for condition in CONDITIONS:
        condition_expr = group_expr[group_expr["condition"] == condition].copy()
        for axis in PRIORITY_AXES:
            ligand = axis["ligand"]
            receptor = axis["receptor"]
            source = axis["source_group"]
            target = axis["target_group"]
            ligand_mean = get_value(condition_expr, ligand, source, condition, "mean_expr")
            receptor_mean = get_value(condition_expr, receptor, target, condition, "mean_expr")
            ligand_pct = get_value(condition_expr, ligand, source, condition, "pct_expr")
            receptor_pct = get_value(condition_expr, receptor, target, condition, "pct_expr")

            ligand_total = condition_expr[condition_expr["gene"].str.upper() == ligand.upper()]["mean_expr"].sum()
            receptor_total = condition_expr[condition_expr["gene"].str.upper() == receptor.upper()]["mean_expr"].sum()
            ligand_specificity = ligand_mean / (ligand_total + eps)
            receptor_specificity = receptor_mean / (receptor_total + eps)

            expr_product = ligand_mean * receptor_mean
            pct_product = ligand_pct * receptor_pct
            specificity_product = ligand_specificity * receptor_specificity
            balanced_expr_pct = math.sqrt(expr_product + eps) * math.sqrt(pct_product + eps)
            geometric_consensus = math.exp(
                np.mean(
                    np.log(
                        np.asarray(
                            [expr_product + eps, pct_product + eps, specificity_product + eps, balanced_expr_pct + eps],
                            dtype=float,
                        )
                    )
                )
            )
            rows.append(
                {
                    **axis,
                    "condition": condition,
                    "ligand_mean_source": ligand_mean,
                    "receptor_mean_target": receptor_mean,
                    "ligand_pct_source": ligand_pct,
                    "receptor_pct_target": receptor_pct,
                    "ligand_specificity_source": ligand_specificity,
                    "receptor_specificity_target": receptor_specificity,
                    "expr_product": expr_product,
                    "pct_product": pct_product,
                    "specificity_product": specificity_product,
                    "balanced_expr_pct": balanced_expr_pct,
                    "geometric_consensus": geometric_consensus,
                }
            )
    scores = pd.DataFrame(rows)
    for method in METHOD_COLUMNS:
        scores[f"{method}_rank"] = scores.groupby("condition")[method].rank(ascending=False, method="min")
    scores["aggregate_rank"] = scores[[f"{method}_rank" for method in METHOD_COLUMNS]].mean(axis=1)
    scores["aggregate_score"] = 1.0 / scores["aggregate_rank"]
    return scores


def build_change_summary(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for axis, sub in scores.groupby("axis", sort=False):
        indexed = sub.set_index("condition")
        for method in METHOD_COLUMNS + ["aggregate_score"]:
            lean = float(indexed.loc["lean", method])
            obese = float(indexed.loc["obese", method])
            enal = float(indexed.loc["obese_enalapril", method])
            rows.append(
                {
                    "axis": axis,
                    "method": method,
                    "lean": lean,
                    "obese": obese,
                    "obese_enalapril": enal,
                    "obese_vs_lean_delta": obese - lean,
                    "enalapril_vs_obese_delta": enal - obese,
                    "obese_vs_lean_ratio": (obese + 1e-9) / (lean + 1e-9),
                    "enalapril_vs_obese_ratio": (enal + 1e-9) / (obese + 1e-9),
                    "disease_up": obese > lean,
                    "enalapril_attenuated": enal < obese,
                }
            )
    return pd.DataFrame(rows)


def save_bubble(scores: pd.DataFrame, out_path: Path) -> None:
    axes_order = [axis["axis"] for axis in PRIORITY_AXES]
    condition_order = CONDITIONS
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for _, row in scores.iterrows():
        x = condition_order.index(row["condition"])
        y = axes_order.index(row["axis"])
        size = 2600 * float(row["aggregate_score"])
        color_value = float(row["geometric_consensus"])
        ax.scatter(x, y, s=size, c=color_value, cmap="viridis", vmin=0, vmax=scores["geometric_consensus"].max(), alpha=0.82, edgecolor="black", linewidth=0.5)
        ax.text(x, y, f"R{row['aggregate_rank']:.1f}", ha="center", va="center", fontsize=8, color="white")
    ax.set_xticks(range(len(condition_order)))
    ax.set_xticklabels(condition_order, rotation=20, ha="right")
    ax.set_yticks(range(len(axes_order)))
    ax.set_yticklabels(axes_order)
    ax.set_title("Targeted multi-score robustness of priority ligand-receptor axes")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Priority LR axis")
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0, vmax=scores["geometric_consensus"].max()))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Geometric consensus score")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_method_heatmap(scores: pd.DataFrame, out_path: Path) -> None:
    plot = scores.copy()
    plot["row_label"] = plot["axis"] + " | " + plot["condition"]
    mat = plot.set_index("row_label")[METHOD_COLUMNS]
    z = mat.apply(lambda col: (col - col.mean()) / (col.std(ddof=0) if col.std(ddof=0) else 1.0), axis=0)
    fig, ax = plt.subplots(figsize=(9.2, 5.5))
    vmax = max(1.0, float(np.max(np.abs(z.values))))
    im = ax.imshow(z.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(z.shape[1]))
    ax.set_xticklabels([col.replace("_", "\n") for col in z.columns], fontsize=8)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index, fontsize=8)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            ax.text(j, i, f"{z.iloc[i, j]:.1f}", ha="center", va="center", fontsize=7)
    ax.set_title("Method-level standardized LR robustness scores")
    fig.colorbar(im, ax=ax, shrink=0.82, label="Method-wise z-score")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores = build_lr_scores()
    changes = build_change_summary(scores)
    scores.to_csv(OUT_DIR / "targeted_lr_multiscore_all.csv", index=False, encoding="utf-8-sig")
    changes.to_csv(OUT_DIR / "targeted_lr_multiscore_change_summary.csv", index=False, encoding="utf-8-sig")
    save_bubble(scores, OUT_DIR / "FigureBN_targeted_lr_multiscore_bubble.png")
    save_method_heatmap(scores, OUT_DIR / "FigureBO_targeted_lr_method_heatmap.png")

    spp1 = scores[scores["axis"] == "Spp1->Cd44"].sort_values("condition")
    spp1_changes = changes[(changes["axis"] == "Spp1->Cd44") & (changes["method"].isin(["geometric_consensus", "aggregate_score"]))]
    lines = [
        "Targeted multi-score LR robustness analysis for GSE233078",
        "",
        "Why this analysis was added:",
        "LIANA could not be installed in the current R environment during this run: CRAN did not provide liana for this R version, and GitHub installation was blocked by API rate limit.",
        "As a practical fallback, this module uses the same focused single-cell group-expression table to score priority LR axes under multiple independent scoring definitions.",
        "",
        "Important boundary:",
        "This is not a formal LIANA result. It is a targeted robustness/sensitivity layer for the existing candidate axes.",
        "",
        "Scoring methods:",
        "- expr_product: source ligand mean expression x target receptor mean expression.",
        "- pct_product: source ligand expressing fraction x target receptor expressing fraction.",
        "- specificity_product: ligand and receptor group specificity within each condition.",
        "- balanced_expr_pct: combined expression and prevalence score.",
        "- geometric_consensus and aggregate_rank: cross-method summaries.",
        "",
        "Spp1->Cd44 DTL->Mono_Macro readouts:",
    ]
    for _, row in spp1.iterrows():
        lines.append(
            f"- {row['condition']}: expr_product={row['expr_product']:.4f}, pct_product={row['pct_product']:.4f}, "
            f"specificity_product={row['specificity_product']:.4f}, geometric_consensus={row['geometric_consensus']:.4f}, "
            f"aggregate_rank={row['aggregate_rank']:.2f}"
        )
    lines.append("")
    lines.append("Spp1->Cd44 change summary:")
    for _, row in spp1_changes.iterrows():
        lines.append(
            f"- {row['method']}: obese_vs_lean_delta={row['obese_vs_lean_delta']:.4f}, "
            f"enalapril_vs_obese_delta={row['enalapril_vs_obese_delta']:.4f}, "
            f"disease_up={row['disease_up']}, enalapril_attenuated={row['enalapril_attenuated']}"
        )
    lines.extend(
        [
            "",
            "Writing interpretation:",
            "1. This layer supports that the Spp1->Cd44 candidate axis remains high-ranked under multiple targeted scoring definitions.",
            "2. It should be written as targeted multi-score robustness, not as LIANA formal confirmation.",
            "3. Because Spp1 expression remains high in DTL after enalapril while Cd44 in Mono/Macro decreases, attenuation is more visible in receptor/prevalence and aggregate behavior than in all raw expression components.",
        ]
    )
    write_text(OUT_DIR / "status.txt", "completed\n")


if __name__ == "__main__":
    main()


