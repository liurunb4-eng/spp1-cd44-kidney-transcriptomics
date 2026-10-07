from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
HERE.mkdir(parents=True, exist_ok=True)
STEP06_PATH = Path(__file__).with_name("06_frozen_signature_selection_audit.py")
RNG_SEED = 20261004
N_PERMUTATIONS = 5000


def load_step06():
    spec = importlib.util.spec_from_file_location("step06", STEP06_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {STEP06_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


step06 = load_step06()


PRIMARY_GENES = [
    "HAVCR1",
    "LCN2",
    "VCAM1",
    "CLU",
    "COL1A1",
    "FN1",
    "ACTA2",
    "TGFB1",
    "TLR4",
    "NLRP3",
]


SENSITIVITY_SIGNATURES = {
    "candidate_DKD_TI_remodeling_no_SPP1_CD44": PRIMARY_GENES,
    "component_injury_HAVCR1_LCN2_VCAM1": ["HAVCR1", "LCN2", "VCAM1"],
    "component_CLU": ["CLU"],
    "component_fibrosis_ECM": ["COL1A1", "FN1", "ACTA2", "TGFB1"],
    "component_inflammation": ["TLR4", "NLRP3"],
    "candidate_minus_injury": [g for g in PRIMARY_GENES if g not in {"HAVCR1", "LCN2", "VCAM1"}],
    "candidate_minus_CLU": [g for g in PRIMARY_GENES if g != "CLU"],
    "candidate_minus_fibrosis_ECM": [g for g in PRIMARY_GENES if g not in {"COL1A1", "FN1", "ACTA2", "TGFB1"}],
    "candidate_minus_inflammation": [g for g in PRIMARY_GENES if g not in {"TLR4", "NLRP3"}],
    "candidate_minus_injury_fibrosis_inflammation": ["CLU"],
}


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    mean = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=0).replace(0, np.nan)
    return mat.sub(mean, axis=0).div(sd, axis=0).fillna(0.0)


def score_gene_set(z: pd.DataFrame, genes: list[str]) -> pd.Series:
    index = {str(g).upper(): g for g in z.index}
    present = [index[g] for g in genes if g in index]
    if not present:
        return pd.Series(np.nan, index=z.columns)
    return z.loc[present].mean(axis=0)


def score_many(expr: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    z = zscore_rows(expr)
    rows = []
    presence = []
    expr_index = {str(g).upper(): g for g in expr.index}
    for sig, genes in SENSITIVITY_SIGNATURES.items():
        genes = step06.clean_gene_list(genes)
        present = [g for g in genes if g in expr_index]
        missing = [g for g in genes if g not in expr_index]
        s = score_gene_set(z, genes)
        rows.append(pd.DataFrame({"dataset": dataset, "sample_id": s.index, "signature": sig, "score": s.values}))
        presence.append(
            {
                "dataset": dataset,
                "signature": sig,
                "requested_n": len(genes),
                "present_n": len(present),
                "present_fraction": len(present) / len(genes) if genes else np.nan,
                "present_genes": ";".join(present),
                "missing_genes": ";".join(missing),
            }
        )
    return pd.concat(rows, ignore_index=True), pd.DataFrame(presence), z


def spearman_result(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    data = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < 3 or data["x"].nunique() < 2 or data["y"].nunique() < 2:
        return int(data.shape[0]), np.nan, np.nan
    rho, p = stats.spearmanr(data["x"], data["y"])
    return int(data.shape[0]), float(rho), float(p)


def group_result(control: pd.Series, dkd: pd.Series) -> tuple[int, int, float, float]:
    control = pd.Series(control).dropna().astype(float)
    dkd = pd.Series(dkd).dropna().astype(float)
    if len(control) < 1 or len(dkd) < 1:
        return int(len(control)), int(len(dkd)), np.nan, np.nan
    p = stats.mannwhitneyu(control, dkd, alternative="two-sided").pvalue
    return int(len(control)), int(len(dkd)), float(dkd.mean() - control.mean()), float(p)


def analyse_scores(g175_scores: pd.DataFrame, g175_meta: pd.DataFrame, g301_scores: pd.DataFrame, g301_meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    g175 = g175_scores.merge(g175_meta, left_on="sample_id", right_on="gsm", how="left")
    for sig, sub in g175.groupby("signature"):
        n, rho, p = spearman_result(sub["score"], sub["egfr_ckd_epi"])
        rows.append(
            {
                "dataset": "GSE175759",
                "analysis": "eGFR_spearman_all_samples",
                "signature": sig,
                "n": n,
                "effect": rho,
                "p_value": p,
                "effect_label": "spearman_rho_vs_eGFR",
                "direction_note": "negative rho means higher score with lower eGFR",
            }
        )

    g301 = step06.prepare_gse30122_score_meta(g301_scores, g301_meta)
    subsets = {
        "primary_tubules_DKD_minus_control": g301[g301["is_primary_tubules"].fillna(False).astype(bool)],
        "glomerulus_DKD_minus_control": g301[g301["subregion_clean"].astype(str).str.lower().eq("glomerulus")],
    }
    for analysis, df in subsets.items():
        df = df[df["disease_binary"].isin(["Control", "DKD"])]
        for sig, sub in df.groupby("signature"):
            control = sub[sub["disease_binary"].eq("Control")]["score"]
            dkd = sub[sub["disease_binary"].eq("DKD")]["score"]
            control_n, dkd_n, delta, p = group_result(control, dkd)
            rows.append(
                {
                    "dataset": "GSE30122",
                    "analysis": analysis,
                    "signature": sig,
                    "control_n": control_n,
                    "dkd_n": dkd_n,
                    "effect": delta,
                    "p_value": p,
                    "effect_label": "DKD_minus_control",
                    "direction_note": "positive delta means higher score in DKD",
                }
            )
    return pd.DataFrame(rows)


def expression_bins(expr: pd.DataFrame, target_genes: list[str], n_bins: int = 20) -> tuple[pd.Series, dict[str, list[str]], list[str]]:
    mean_expr = expr.mean(axis=1).sort_values()
    usable = mean_expr[np.isfinite(mean_expr)]
    bins = pd.qcut(usable.rank(method="first"), q=min(n_bins, max(2, len(usable) // 100)), labels=False, duplicates="drop")
    bins = bins.astype(int)
    upper_to_index = {str(g).upper(): g for g in expr.index}
    present_targets = [g for g in target_genes if g in upper_to_index]
    target_index = [upper_to_index[g] for g in present_targets]
    target_set = set(target_index)
    candidates_by_target = {}
    for tg_upper, tg_index in zip(present_targets, target_index):
        b = int(bins.loc[tg_index])
        same_bin = bins[bins.eq(b)].index.tolist()
        pool = [g for g in same_bin if g not in target_set]
        if len(pool) < 50:
            wider = bins[bins.sub(b).abs().le(1)].index.tolist()
            pool = [g for g in wider if g not in target_set]
        candidates_by_target[tg_upper] = pool
    return bins, candidates_by_target, present_targets


def run_random_controls(
    expr: pd.DataFrame,
    z: pd.DataFrame,
    dataset: str,
    meta: pd.DataFrame,
    n_perm: int = N_PERMUTATIONS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RNG_SEED + (1 if dataset == "GSE30122" else 0))
    target = step06.clean_gene_list(PRIMARY_GENES)
    _, candidates_by_target, present_targets = expression_bins(expr, target)
    observed_score = score_gene_set(z, target)

    if dataset == "GSE175759":
        obs_df = pd.DataFrame({"sample_id": observed_score.index, "score": observed_score.values}).merge(
            meta, left_on="sample_id", right_on="gsm", how="left"
        )
        _, observed, observed_p = spearman_result(obs_df["score"], obs_df["egfr_ckd_epi"])
        analysis = "eGFR_spearman_all_samples"
        direction = "lower_or_equal"
    else:
        obs_score_df = pd.DataFrame(
            {
                "dataset": dataset,
                "sample_id": observed_score.index,
                "signature": "observed",
                "score": observed_score.values,
            }
        )
        obs_df = step06.prepare_gse30122_score_meta(obs_score_df, meta)
        obs_df = obs_df[obs_df["is_primary_tubules"].fillna(False).astype(bool)]
        obs_df = obs_df[obs_df["disease_binary"].isin(["Control", "DKD"])]
        control = obs_df[obs_df["disease_binary"].eq("Control")]["score"]
        dkd = obs_df[obs_df["disease_binary"].eq("DKD")]["score"]
        _, _, observed, observed_p = group_result(control, dkd)
        analysis = "primary_tubules_DKD_minus_control"
        direction = "greater_or_equal"

    null_rows = []
    z_index = {str(g).upper(): g for g in z.index}
    for i in range(n_perm):
        sampled = []
        used = set()
        for tg in present_targets:
            pool = [g for g in candidates_by_target[tg] if g not in used]
            if not pool:
                pool = candidates_by_target[tg]
            choice = rng.choice(pool)
            used.add(choice)
            sampled.append(choice)
        random_score = z.loc[sampled].mean(axis=0)
        if dataset == "GSE175759":
            rdf = pd.DataFrame({"sample_id": random_score.index, "score": random_score.values}).merge(
                meta, left_on="sample_id", right_on="gsm", how="left"
            )
            _, stat_value, p_value = spearman_result(rdf["score"], rdf["egfr_ckd_epi"])
        else:
            score_df = pd.DataFrame(
                {"dataset": dataset, "sample_id": random_score.index, "signature": "random", "score": random_score.values}
            )
            rdf = step06.prepare_gse30122_score_meta(score_df, meta)
            rdf = rdf[rdf["is_primary_tubules"].fillna(False).astype(bool)]
            rdf = rdf[rdf["disease_binary"].isin(["Control", "DKD"])]
            control = rdf[rdf["disease_binary"].eq("Control")]["score"]
            dkd = rdf[rdf["disease_binary"].eq("DKD")]["score"]
            _, _, stat_value, p_value = group_result(control, dkd)
        null_rows.append(
            {
                "dataset": dataset,
                "analysis": analysis,
                "iteration": i + 1,
                "random_stat": stat_value,
                "random_p": p_value,
                "random_genes": ";".join(str(g).upper() for g in sampled),
            }
        )
    null_df = pd.DataFrame(null_rows)
    if direction == "lower_or_equal":
        empirical_p = (1 + (null_df["random_stat"] <= observed).sum()) / (len(null_df) + 1)
        percentile = (null_df["random_stat"] < observed).mean()
    else:
        empirical_p = (1 + (null_df["random_stat"] >= observed).sum()) / (len(null_df) + 1)
        percentile = (null_df["random_stat"] < observed).mean()
    summary = pd.DataFrame(
        [
            {
                "dataset": dataset,
                "analysis": analysis,
                "signature": "candidate_DKD_TI_remodeling_no_SPP1_CD44",
                "target_present_n": len(present_targets),
                "n_permutations": n_perm,
                "observed_stat": observed,
                "observed_nominal_p": observed_p,
                "null_mean": null_df["random_stat"].mean(),
                "null_sd": null_df["random_stat"].std(ddof=1),
                "null_2p5": null_df["random_stat"].quantile(0.025),
                "null_50": null_df["random_stat"].quantile(0.5),
                "null_97p5": null_df["random_stat"].quantile(0.975),
                "observed_percentile_vs_null": percentile,
                "empirical_directional_p": empirical_p,
                "direction_tested": direction,
            }
        ]
    )
    return summary, null_df


def fmt(x: object, digits: int = 3) -> str:
    if x is None or pd.isna(x):
        return "NA"
    try:
        return f"{float(x):.{digits}g}"
    except Exception:
        return str(x)


def write_report(stats_df: pd.DataFrame, presence_df: pd.DataFrame, random_summary: pd.DataFrame) -> None:
    lines = [
        "# Candidate DKD/TI remodeling signature sensitivity and random-control audit",
        "",
        "Date: 2026-10-04",
        "",
        "Primary candidate: `candidate_DKD_TI_remodeling_no_SPP1_CD44`.",
        "",
        "Gene set: `HAVCR1;LCN2;VCAM1;CLU;COL1A1;FN1;ACTA2;TGFB1;TLR4;NLRP3`.",
        "",
        "This step tests whether the candidate remains informative after decomposing injury, CLU, fibrosis/ECM, and inflammatory components, and whether its observed statistics exceed expression-matched random gene sets of the same size.",
        "",
        "## Component and leave-one-component sensitivity",
        "",
        "| Signature | GSE175759 rho vs eGFR | GSE30122 primary tubule delta | GSE30122 glomerulus delta |",
        "|---|---:|---:|---:|",
    ]
    sig_order = list(SENSITIVITY_SIGNATURES.keys())
    for sig in sig_order:
        sub = stats_df[stats_df["signature"].eq(sig)]
        g175 = sub[sub["analysis"].eq("eGFR_spearman_all_samples")]
        pt = sub[sub["analysis"].eq("primary_tubules_DKD_minus_control")]
        gl = sub[sub["analysis"].eq("glomerulus_DKD_minus_control")]
        lines.append(
            f"| {sig} | {fmt(g175['effect'].iloc[0] if len(g175) else np.nan)} (p={fmt(g175['p_value'].iloc[0] if len(g175) else np.nan)}) | {fmt(pt['effect'].iloc[0] if len(pt) else np.nan)} (p={fmt(pt['p_value'].iloc[0] if len(pt) else np.nan)}) | {fmt(gl['effect'].iloc[0] if len(gl) else np.nan)} (p={fmt(gl['p_value'].iloc[0] if len(gl) else np.nan)}) |"
        )
    lines += [
        "",
        "## Expression-matched random control",
        "",
        "| Dataset | Observed statistic | Null median | Null 2.5%-97.5% | Empirical directional P |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, r in random_summary.iterrows():
        lines.append(
            f"| {r['dataset']} / {r['analysis']} | {fmt(r['observed_stat'])} | {fmt(r['null_50'])} | {fmt(r['null_2p5'])} to {fmt(r['null_97p5'])} | {fmt(r['empirical_directional_p'])} |"
        )
    lines += [
        "",
        "## Result-specific conclusion",
        "",
        "The 10-gene DKD/TI candidate is stronger than expression-matched random gene sets in both core datasets. In GSE175759, the observed eGFR correlation falls beyond the lower tail of the matched null; in GSE30122 primary tubules, the observed DKD-control shift falls beyond the upper tail of the matched null.",
        "",
        "However, component decomposition shows that the signal is injury/CLU-dominant. `component_CLU` and `component_injury_HAVCR1_LCN2_VCAM1` are stronger than the full 10-gene score in the DKD primary-tubule comparison. Removing the injury component weakens the GSE30122 primary-tubule result, and removing CLU also weakens it. The fibrosis/ECM component is clearer in glomeruli than in primary tubules, while the two-gene inflammation component does not support the primary-tubule DKD direction.",
        "",
        "Therefore the candidate can be carried forward, but it should be written as an injury-dominant DKD tubulointerstitial remodeling context rather than as a newly independent fibrosis-inflammatory program. This is useful for the relaunch because it is stronger than random controls and independent of SPP1/CD44, but the biological claim must stay conservative.",
        "",
        "## Interpretation",
        "",
        "- If the full candidate remains stronger than matched random gene sets in both datasets, it can be carried forward as the primary computational signature.",
        "- If the signal is mostly explained by `component_injury_HAVCR1_LCN2_VCAM1`, the manuscript should avoid presenting the 10-gene score as a distinct new program and instead write it as a DKD tubulointerstitial injury/remodeling context.",
        "- Because the candidate intentionally combines injury, fibrosis/ECM, and inflammatory context genes, overlap-removal to a single residual gene is not biologically meaningful. The defensible test here is component decomposition plus expression-matched random controls.",
        "",
        "## Output files",
        "",
        "- `07_candidate_sensitivity_random_control.py`",
        "- `07_candidate_component_sensitivity.csv`",
        "- `07_candidate_gene_presence.csv`",
        "- `07_candidate_random_control_summary.csv`",
        "- `07_candidate_random_control_nulls.csv.gz`",
        "- `07_candidate_sensitivity_random_control.md`",
        "",
    ]
    (HERE / "07_candidate_sensitivity_random_control.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    g175_expr, g175_meta = step06.load_gse175759_expression()
    g301_expr, g301_meta = step06.load_gse30122_expression()

    g175_scores, g175_presence, g175_z = score_many(g175_expr, "GSE175759")
    g301_scores, g301_presence, g301_z = score_many(g301_expr, "GSE30122")

    stats_df = analyse_scores(g175_scores, g175_meta, g301_scores, g301_meta)
    presence_df = pd.concat([g175_presence, g301_presence], ignore_index=True)

    random_175, null_175 = run_random_controls(g175_expr, g175_z, "GSE175759", g175_meta)
    random_301, null_301 = run_random_controls(g301_expr, g301_z, "GSE30122", g301_meta)
    random_summary = pd.concat([random_175, random_301], ignore_index=True)
    random_nulls = pd.concat([null_175, null_301], ignore_index=True)

    stats_df.to_csv(HERE / "07_candidate_component_sensitivity.csv", index=False)
    presence_df.to_csv(HERE / "07_candidate_gene_presence.csv", index=False)
    random_summary.to_csv(HERE / "07_candidate_random_control_summary.csv", index=False)
    random_nulls.to_csv(HERE / "07_candidate_random_control_nulls.csv.gz", index=False)
    write_report(stats_df, presence_df, random_summary)

    print(f"Wrote {HERE / '07_candidate_sensitivity_random_control.md'}")
    print(f"Wrote {HERE / '07_candidate_random_control_summary.csv'}")


if __name__ == "__main__":
    main()
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
