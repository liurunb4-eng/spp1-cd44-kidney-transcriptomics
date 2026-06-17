#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv
import math
import shutil

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.csgraph import dijkstra
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_ORDER = ["lean", "obese", "obese_enalapril"]
FOCUS_ORDER = ["DTL", "Mono_Macro"]
SCORE_COLS = ["spp1_cd44_axis", "tubular_stress", "myeloid_activation", "fibro_inflammatory_context"]


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def zscore(values: pd.Series) -> pd.Series:
    mu = float(values.mean())
    sd = float(values.std(ddof=0))
    if sd == 0 or math.isnan(sd):
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - mu) / sd


def add_state_axis(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in SCORE_COLS:
        out[f"{col}_z"] = out.groupby("focus_group")[col].transform(zscore)

    is_dtl = out["focus_group"] == "DTL"
    out["state_axis_z"] = np.nan
    out.loc[is_dtl, "state_axis_z"] = out.loc[
        is_dtl,
        ["spp1_cd44_axis_z", "tubular_stress_z", "fibro_inflammatory_context_z"],
    ].mean(axis=1)
    out.loc[~is_dtl, "state_axis_z"] = out.loc[
        ~is_dtl,
        ["spp1_cd44_axis_z", "myeloid_activation_z", "fibro_inflammatory_context_z"],
    ].mean(axis=1)
    return out


def choose_root_cells(group_df: pd.DataFrame, quantile: float = 0.20) -> np.ndarray:
    lean_mask = group_df["condition"] == "lean"
    if lean_mask.any():
        lean_scores = group_df.loc[lean_mask, "state_axis_z"]
        cutoff = float(lean_scores.quantile(quantile))
        roots = np.where(lean_mask.to_numpy() & (group_df["state_axis_z"].to_numpy() <= cutoff))[0]
        if len(roots) > 0:
            return roots
    fallback_cutoff = float(group_df["state_axis_z"].quantile(0.10))
    return np.where(group_df["state_axis_z"].to_numpy() <= fallback_cutoff)[0]


def graph_pseudotime(group_df: pd.DataFrame, n_neighbors: int = 25) -> pd.DataFrame:
    if len(group_df) < 3:
        out = group_df.copy()
        out["graph_pseudotime"] = 0.0
        out["trajectory_root"] = out["condition"] == "lean"
        return out

    feature_cols = ["UMAP_1", "UMAP_2"] + SCORE_COLS
    features = group_df[feature_cols].astype(float).to_numpy()
    scaled = StandardScaler().fit_transform(features)
    k = min(n_neighbors + 1, len(group_df))
    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nbrs.fit(scaled)
    distances, indices = nbrs.kneighbors(scaled)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for idx in range(len(group_df)):
        for distance, neighbor in zip(distances[idx, 1:], indices[idx, 1:]):
            rows.extend([idx, int(neighbor)])
            cols.extend([int(neighbor), idx])
            data.extend([float(distance) + 1e-6, float(distance) + 1e-6])
    graph = csr_matrix((data, (rows, cols)), shape=(len(group_df), len(group_df)))

    roots = choose_root_cells(group_df)
    super_root = len(group_df)
    graph_with_root = lil_matrix((len(group_df) + 1, len(group_df) + 1))
    graph_with_root[: len(group_df), : len(group_df)] = graph
    for root in roots:
        graph_with_root[super_root, int(root)] = 1e-9
        graph_with_root[int(root), super_root] = 1e-9

    distances_from_root = dijkstra(graph_with_root.tocsr(), directed=False, indices=super_root)[: len(group_df)]
    finite = np.isfinite(distances_from_root)
    if finite.any():
        max_finite = float(np.max(distances_from_root[finite]))
        distances_from_root[~finite] = max_finite
        pseudotime = distances_from_root / max_finite if max_finite > 0 else np.zeros(len(group_df))
    else:
        pseudotime = np.zeros(len(group_df))

    out = group_df.copy()
    out["graph_pseudotime"] = pseudotime
    out["trajectory_root"] = False
    out.iloc[roots, out.columns.get_loc("trajectory_root")] = True
    return out


def build_cell_pseudotime(df: pd.DataFrame) -> pd.DataFrame:
    state_df = add_state_axis(df)
    pieces = []
    for focus_group in FOCUS_ORDER:
        group_df = state_df[state_df["focus_group"] == focus_group].copy()
        pieces.append(graph_pseudotime(group_df))
    return pd.concat(pieces, ignore_index=True)


def build_sample_summary(cell_df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        cell_df.groupby(["sample_id", "condition", "focus_group"], as_index=False)
        .agg(
            n_cells=("cell_id", "count"),
            graph_pseudotime_mean=("graph_pseudotime", "mean"),
            state_axis_z_mean=("state_axis_z", "mean"),
            root_fraction=("trajectory_root", "mean"),
        )
        .sort_values(["focus_group", "condition", "sample_id"])
    )
    return grouped


def build_condition_summary(sample_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for focus_group in FOCUS_ORDER:
        for condition in CONDITION_ORDER:
            subset = sample_df[(sample_df["focus_group"] == focus_group) & (sample_df["condition"] == condition)]
            rows.append(
                {
                    "focus_group": focus_group,
                    "condition": condition,
                    "n_samples": int(len(subset)),
                    "total_cells": int(subset["n_cells"].sum()) if len(subset) else 0,
                    "graph_pseudotime_mean": float(subset["graph_pseudotime_mean"].mean()) if len(subset) else 0.0,
                    "state_axis_z_mean": float(subset["state_axis_z_mean"].mean()) if len(subset) else 0.0,
                    "root_fraction_mean": float(subset["root_fraction"].mean()) if len(subset) else 0.0,
                }
            )
    return pd.DataFrame(rows)


def classify_trajectory(lean: float, obese: float, enalapril: float) -> str:
    if obese > lean and enalapril < obese:
        if enalapril <= lean:
            return "obese_late_with_full_return"
        return "obese_late_with_partial_return"
    if obese > lean:
        return "obese_late_without_return"
    if obese < lean and enalapril > obese:
        return "obese_early_with_recovery"
    return "mixed_or_flat"


def build_trajectory_calls(condition_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for focus_group in FOCUS_ORDER:
        values = {
            row["condition"]: float(row["graph_pseudotime_mean"])
            for _, row in condition_df[condition_df["focus_group"] == focus_group].iterrows()
        }
        lean = values.get("lean", 0.0)
        obese = values.get("obese", 0.0)
        enalapril = values.get("obese_enalapril", 0.0)
        rows.append(
            {
                "focus_group": focus_group,
                "lean_graph_pseudotime": lean,
                "obese_graph_pseudotime": obese,
                "obese_enalapril_graph_pseudotime": enalapril,
                "obese_vs_lean_delta": obese - lean,
                "enalapril_vs_obese_delta": enalapril - obese,
                "trajectory_call": classify_trajectory(lean, obese, enalapril),
            }
        )
    return pd.DataFrame(rows)


def save_pseudotime_map(cell_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, focus_group in zip(axes, FOCUS_ORDER):
        subset = cell_df[cell_df["focus_group"] == focus_group]
        sc = ax.scatter(
            subset["UMAP_1"],
            subset["UMAP_2"],
            c=subset["graph_pseudotime"],
            s=4,
            cmap="viridis",
            linewidths=0,
            alpha=0.85,
        )
        roots = subset[subset["trajectory_root"]]
        ax.scatter(roots["UMAP_1"], roots["UMAP_2"], s=8, c="#f97316", linewidths=0, alpha=0.85, label="root")
        ax.set_title(f"{focus_group} graph pseudotime")
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="pseudotime")
    fig.suptitle("GSE233078 DTL / Mono_Macro graph-based pseudotime", fontsize=14, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def save_condition_plot(sample_df: pd.DataFrame, out_path: Path) -> None:
    color_map = {"lean": "#16a34a", "obese": "#dc2626", "obese_enalapril": "#2563eb"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, constrained_layout=True)
    rng = np.random.default_rng(20260420)
    for ax, focus_group in zip(axes, FOCUS_ORDER):
        subset = sample_df[sample_df["focus_group"] == focus_group]
        for idx, condition in enumerate(CONDITION_ORDER):
            values = subset.loc[subset["condition"] == condition, "graph_pseudotime_mean"].astype(float).to_numpy()
            if len(values) == 0:
                continue
            jitter = rng.normal(0, 0.035, size=len(values))
            ax.scatter(
                np.full(len(values), idx) + jitter,
                values,
                s=52,
                color=color_map[condition],
                edgecolor="white",
                linewidth=0.8,
                alpha=0.9,
                label=condition if focus_group == FOCUS_ORDER[0] else None,
            )
            ax.plot([idx - 0.18, idx + 0.18], [values.mean(), values.mean()], color="#111827", linewidth=2.2)
        ax.set_title(f"{focus_group} sample-aware pseudotime")
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels(CONDITION_ORDER, rotation=20, ha="right")
        ax.set_ylabel("Mean graph pseudotime")
        ax.set_ylim(-0.03, 1.03)
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.8)
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle("Disease-associated trajectory shift and enalapril return", fontsize=14, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def format_float(value: float) -> str:
    return f"{float(value):.6f}"


def build_summary(condition_df: pd.DataFrame, calls_df: pd.DataFrame) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def build_readable_result(summary_text: str) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_to_quick(project_root: Path, out_dir: Path) -> None:
    return None

def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    in_path = project_root / "results/public_data_gse233078_dtl_monomacro_state_shift/dtl_monomacro_cell_state_scores.csv"
    out_dir = project_root / "results/public_data_gse233078_dtl_monomacro_graph_pseudotime"
    out_dir.mkdir(parents=True, exist_ok=True)

    cell_input = pd.read_csv(in_path)
    cell_df = build_cell_pseudotime(cell_input)
    sample_df = build_sample_summary(cell_df)
    condition_df = build_condition_summary(sample_df)
    calls_df = build_trajectory_calls(condition_df)

    cell_df.to_csv(out_dir / "dtl_monomacro_cell_graph_pseudotime.csv", index=False, encoding="utf-8-sig")
    sample_df.to_csv(out_dir / "dtl_monomacro_sample_graph_pseudotime.csv", index=False, encoding="utf-8-sig")
    condition_df.to_csv(out_dir / "dtl_monomacro_condition_graph_pseudotime.csv", index=False, encoding="utf-8-sig")
    calls_df.to_csv(out_dir / "dtl_monomacro_graph_pseudotime_calls.csv", index=False, encoding="utf-8-sig")

    save_pseudotime_map(cell_df, out_dir / "FigureAD_gse233078_dtl_monomacro_graph_pseudotime_map.svg")
    save_condition_plot(sample_df, out_dir / "FigureAE_gse233078_dtl_monomacro_graph_pseudotime_condition.svg")

    summary = build_summary(condition_df, calls_df)
    write_text(out_dir / "gse233078_dtl_monomacro_graph_pseudotime_summary.txt", summary)
    write_text(out_dir / "status.txt", "GSE233078 DTL / Mono_Macro graph pseudotime completed.\n")
    print(f"Graph pseudotime outputs written to: {out_dir}")


if __name__ == "__main__":
    main()


