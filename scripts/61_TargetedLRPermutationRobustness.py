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
OUT_DIR = PROJECT_ROOT / "results" / "targeted_lr_permutation_robustness_spp1_cd44"

N_PERM = 1000
RANDOM_SEED = 20260426

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


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def read_small_mtx(path: Path) -> np.ndarray:
    """Read the tiny coordinate MatrixMarket file without relying on scipy's parser."""
    with path.open("r", encoding="utf-8") as handle:
        banner = handle.readline().strip()
        if not banner.startswith("%%MatrixMarket"):
            raise ValueError(f"Not a MatrixMarket file: {path}")
        line = handle.readline().strip()
        while line.startswith("%"):
            line = handle.readline().strip()
        n_rows, n_cols, _n_nonzero = [int(x) for x in line.split()]
        mat = np.zeros((n_rows, n_cols), dtype=float)
        for raw in handle:
            if not raw.strip():
                continue
            row, col, value = raw.split()
            mat[int(row) - 1, int(col) - 1] = float(value)
    return mat


def bh_fdr(p_values: pd.Series) -> pd.Series:
    values = p_values.astype(float).to_numpy()
    n = len(values)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        adjusted[order[i]] = val
        prev = val
    return pd.Series(adjusted, index=p_values.index)


def score_axis(
    gene_matrix: np.ndarray,
    gene_index: dict[str, int],
    groups: np.ndarray,
    axis: dict[str, str],
) -> dict[str, float]:
    ligand = gene_matrix[gene_index[axis["ligand"].upper()], :]
    receptor = gene_matrix[gene_index[axis["receptor"].upper()], :]
    source_mask = groups == axis["source_group"]
    target_mask = groups == axis["target_group"]
    if source_mask.sum() == 0 or target_mask.sum() == 0:
        return {
            "ligand_mean_source": 0.0,
            "receptor_mean_target": 0.0,
            "ligand_pct_source": 0.0,
            "receptor_pct_target": 0.0,
            "expr_product": 0.0,
            "pct_product": 0.0,
            "balanced_expr_pct": 0.0,
        }
    ligand_mean = float(ligand[source_mask].mean())
    receptor_mean = float(receptor[target_mask].mean())
    ligand_pct = float((ligand[source_mask] > 0).mean())
    receptor_pct = float((receptor[target_mask] > 0).mean())
    expr_product = ligand_mean * receptor_mean
    pct_product = ligand_pct * receptor_pct
    balanced = math.sqrt(expr_product + 1e-12) * math.sqrt(pct_product + 1e-12)
    return {
        "ligand_mean_source": ligand_mean,
        "receptor_mean_target": receptor_mean,
        "ligand_pct_source": ligand_pct,
        "receptor_pct_target": receptor_pct,
        "expr_product": expr_product,
        "pct_product": pct_product,
        "balanced_expr_pct": balanced,
    }


def build_permutation_results() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(IN_DIR / "cellchat_mvp_focus_metadata.csv")
    genes = pd.read_csv(IN_DIR / "cellchat_mvp_focus_whitelist_genes.tsv", sep="\t", header=None)[0].astype(str).tolist()
    matrix = read_small_mtx(IN_DIR / "cellchat_mvp_focus_whitelist_counts.mtx")
    if matrix.shape[1] != len(meta):
        raise ValueError(f"Matrix columns ({matrix.shape[1]}) do not match metadata rows ({len(meta)}).")
    gene_index = {gene.upper(): idx for idx, gene in enumerate(genes)}
    rng = np.random.default_rng(RANDOM_SEED)

    observed_rows = []
    null_rows = []
    for condition in CONDITIONS:
        condition_mask = meta["condition"].to_numpy() == condition
        condition_matrix = matrix[:, condition_mask]
        condition_groups = meta.loc[condition_mask, "cellchat_group"].astype(str).to_numpy()
        permuted_groups = [rng.permutation(condition_groups) for _ in range(N_PERM)]
        for axis in PRIORITY_AXES:
            observed = score_axis(condition_matrix, gene_index, condition_groups, axis)
            null_scores = np.asarray(
                [
                    score_axis(condition_matrix, gene_index, perm_groups, axis)["balanced_expr_pct"]
                    for perm_groups in permuted_groups
                ],
                dtype=float,
            )
            null_mean = float(null_scores.mean())
            null_sd = float(null_scores.std(ddof=1))
            empirical_p = float((1 + np.sum(null_scores >= observed["balanced_expr_pct"])) / (N_PERM + 1))
            z_score = float((observed["balanced_expr_pct"] - null_mean) / null_sd) if null_sd > 0 else np.nan
            observed_rows.append(
                {
                    **axis,
                    "condition": condition,
                    **observed,
                    "null_mean_balanced_expr_pct": null_mean,
                    "null_sd_balanced_expr_pct": null_sd,
                    "permutation_z": z_score,
                    "empirical_p": empirical_p,
                    "n_permutations": N_PERM,
                    "source_cell_count": int(np.sum(condition_groups == axis["source_group"])),
                    "target_cell_count": int(np.sum(condition_groups == axis["target_group"])),
                }
            )
            for i, value in enumerate(null_scores):
                null_rows.append(
                    {
                        "axis": axis["axis"],
                        "condition": condition,
                        "permutation_id": i + 1,
                        "null_balanced_expr_pct": float(value),
                    }
                )
    observed_df = pd.DataFrame(observed_rows)
    observed_df["empirical_fdr"] = bh_fdr(observed_df["empirical_p"])
    null_df = pd.DataFrame(null_rows)
    return observed_df, null_df


def save_bubble(observed: pd.DataFrame, out_path: Path) -> None:
    axes_order = [axis["axis"] for axis in PRIORITY_AXES]
    fig, ax = plt.subplots(figsize=(8.8, 4.9))
    max_z = max(1.0, np.nanmax(observed["permutation_z"].to_numpy()))
    for _, row in observed.iterrows():
        x = CONDITIONS.index(row["condition"])
        y = axes_order.index(row["axis"])
        size = 140 + 180 * min(-math.log10(max(float(row["empirical_p"]), 1e-12)), 6)
        color_value = float(row["permutation_z"])
        ax.scatter(
            x,
            y,
            s=size,
            c=color_value,
            cmap="magma",
            vmin=0,
            vmax=max_z,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.45,
        )
        ax.text(x, y, f"p={row['empirical_p']:.3f}", ha="center", va="center", fontsize=7, color="white")
    ax.set_xticks(range(len(CONDITIONS)))
    ax.set_xticklabels(CONDITIONS, rotation=20, ha="right")
    ax.set_yticks(range(len(axes_order)))
    ax.set_yticklabels(axes_order)
    ax.set_title("Within-condition permutation robustness of priority LR axes")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Priority LR axis")
    sm = plt.cm.ScalarMappable(cmap="magma", norm=plt.Normalize(vmin=0, vmax=max_z))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Permutation z-score")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_spp1_null_plot(observed: pd.DataFrame, nulls: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3), sharey=True)
    for ax, condition in zip(axes, CONDITIONS):
        null = nulls[(nulls["axis"] == "Spp1->Cd44") & (nulls["condition"] == condition)]["null_balanced_expr_pct"].to_numpy()
        obs = float(
            observed[(observed["axis"] == "Spp1->Cd44") & (observed["condition"] == condition)]["balanced_expr_pct"].iloc[0]
        )
        p = float(observed[(observed["axis"] == "Spp1->Cd44") & (observed["condition"] == condition)]["empirical_p"].iloc[0])
        ax.hist(null, bins=40, color="#B9C6D3", edgecolor="white")
        ax.axvline(obs, color="#B22727", linewidth=2.2)
        ax.set_title(f"{condition}\nobs={obs:.3f}, p={p:.3f}", fontsize=10)
        ax.set_xlabel("Null balanced score")
    axes[0].set_ylabel("Permutation count")
    fig.suptitle("Spp1->Cd44 DTL-to-Mono/Macro observed score versus condition-matched null")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_change_plot(observed: pd.DataFrame, out_path: Path) -> None:
    plot = observed.pivot(index="axis", columns="condition", values="balanced_expr_pct").loc[
        [axis["axis"] for axis in PRIORITY_AXES],
        CONDITIONS,
    ]
    changes = pd.DataFrame(
        {
            "obese_vs_lean": plot["obese"] - plot["lean"],
            "enalapril_vs_obese": plot["obese_enalapril"] - plot["obese"],
        }
    )
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    vmax = max(0.5, float(np.max(np.abs(changes.to_numpy()))))
    im = ax.imshow(changes.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(changes.shape[1]))
    ax.set_xticklabels(changes.columns, rotation=20, ha="right")
    ax.set_yticks(range(changes.shape[0]))
    ax.set_yticklabels(changes.index)
    for i in range(changes.shape[0]):
        for j in range(changes.shape[1]):
            ax.text(j, i, f"{changes.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Observed LR balanced-score condition deltas")
    fig.colorbar(im, ax=ax, label="Delta balanced score")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    observed, nulls = build_permutation_results()
    observed.to_csv(OUT_DIR / "targeted_lr_permutation_observed.csv", index=False, encoding="utf-8-sig")
    nulls.to_csv(OUT_DIR / "targeted_lr_permutation_null_scores.csv", index=False, encoding="utf-8-sig")
    save_bubble(observed, OUT_DIR / "FigureBT_targeted_lr_permutation_bubble.png")
    save_spp1_null_plot(observed, nulls, OUT_DIR / "FigureBU_spp1_cd44_permutation_null.png")
    save_change_plot(observed, OUT_DIR / "FigureBV_targeted_lr_permutation_condition_deltas.png")

    spp1 = observed[observed["axis"] == "Spp1->Cd44"].set_index("condition").loc[CONDITIONS]
    lines = [
        "GSE233078 targeted LR permutation robustness analysis",
        "",
        "Purpose:",
        "This analysis adds a condition-matched permutation layer to the focused Spp1/Cd44 communication evidence.",
        "It asks whether a priority ligand-receptor source-target assignment is stronger than expected after shuffling cell-group labels within the same condition.",
        "",
        "Boundary:",
        "This is still not formal LIANA / CellPhoneDB output. It is a transparent permutation robustness test based on the exported focused single-cell matrix.",
        "",
        f"Permutation design: {N_PERM} shuffles per condition and priority LR axis; empirical p = (1 + null >= observed) / (N + 1).",
        "",
        "Spp1->Cd44 DTL->Mono_Macro results:",
    ]
    for condition, row in spp1.iterrows():
        lines.append(
            f"- {condition}: balanced_expr_pct={row['balanced_expr_pct']:.4f}, "
            f"null_mean={row['null_mean_balanced_expr_pct']:.4f}, z={row['permutation_z']:.2f}, "
            f"empirical_p={row['empirical_p']:.4f}, empirical_fdr={row['empirical_fdr']:.4f}, "
            f"source_n={int(row['source_cell_count'])}, target_n={int(row['target_cell_count'])}"
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "1. If obese shows a high observed score with a low empirical p value, it supports that the DTL-to-Mono/Macro Spp1-Cd44 assignment is stronger than random cell-group relabeling within obese kidneys.",
            "2. If lean or obese_enalapril are also significant, the result should be interpreted as a condition-specific strength change rather than an on/off switch.",
            "3. The safest writing remains: the Spp1-Cd44 axis is prioritized by focused LR scoring and supported by permutation robustness, but still requires experimental validation for causality.",
        ]
    )
    write_text(OUT_DIR / "status.txt", "completed\n")


if __name__ == "__main__":
    main()


