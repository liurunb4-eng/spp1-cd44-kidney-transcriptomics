from __future__ import annotations

import gzip
import os
from io import StringIO
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

try:
    import statsmodels.formula.api as smf
except Exception:  # pragma: no cover - optional dependency boundary
    smf = None


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(os.environ.get("AJP_DATA_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

S2_PATH = ROOT / "results" / "manuscript_supplementary_tables" / "Supplementary_Table_S2_signature_module_definitions.csv"
GSE175759_COUNTS = ROOT / "results" / "human_bulk_gse175759_clinical_correlation" / "gse175759_merged_counts_by_ensembl.csv.gz"
GSE175759_MAP = ROOT / "results" / "human_bulk_gse175759_clinical_correlation" / "gse175759_ensembl_to_symbol_map.csv"
GSE175759_META = ROOT / "results" / "human_bulk_gse175759_clinical_correlation" / "gse175759_metadata_clean.csv"
GSE175759_EXISTING = ROOT / "results" / "human_bulk_gse175759_clinical_correlation" / "gse175759_priority_signature_scores_wide.csv"

GSE30122_DIR = ROOT / "data" / "public_datasets" / "DKD_bulk_candidate" / "GSE30122"
GSE30122_SERIES = GSE30122_DIR / "GSE30122_series_matrix.txt.gz"
GSE30122_GPL = GSE30122_DIR / "GPL571.soft.gz"
GSE30122_EXISTING = ROOT / "results" / "human_dkd_bulk_gse30122_signature_validation" / "gse30122_signature_scores_long.csv"

NICHENET_GENES = ROOT / "results" / "human_nichenet_targetset_validation" / "nichenet_spp1_top100_targetset_genes.csv"
GSE175759_NICHENET = ROOT / "results" / "human_nichenet_targetset_validation" / "gse175759_nichenet_targetset_scores.csv"
GSE30122_NICHENET = ROOT / "results" / "human_nichenet_targetset_validation" / "gse30122_nichenet_targetset_scores.csv"


def clean_gene_list(value: str | Iterable[str]) -> list[str]:
    if isinstance(value, str):
        raw = value.replace(",", ";").split(";")
    else:
        raw = list(value)
    genes = []
    for g in raw:
        s = str(g).strip().upper()
        if s and s not in genes:
            genes.append(s)
    return genes


def load_candidate_signatures() -> dict[str, dict[str, object]]:
    s2 = pd.read_csv(S2_PATH)
    by_name = {
        row["signature_or_module"]: clean_gene_list(row["gene_list"])
        for _, row in s2.iterrows()
    }
    old_main = by_name["Spp1_Cd44_tubular_immune_program"]
    curated = by_name["curated_fibro_inflammatory_context"]
    fibrosis = by_name["fibrosis_ecm"]
    inflammation = by_name["inflammation"]
    dkd_ti_no_anchor = clean_gene_list(
        "HAVCR1;LCN2;VCAM1;CLU;COL1A1;FN1;ACTA2;TGFB1;TLR4;NLRP3"
    )
    nn = pd.read_csv(NICHENET_GENES)
    nichenet = clean_gene_list(nn["gene"].tolist())

    signatures: dict[str, dict[str, object]] = {
        "old_NDT_composite": {
            "genes": old_main,
            "class": "old rat-derived composite",
            "role": "negative-control candidate; was the NDT main program",
        },
        "old_NDT_composite_without_SPP1_CD44": {
            "genes": [g for g in old_main if g not in {"SPP1", "CD44"}],
            "class": "old rat-derived composite sensitivity",
            "role": "tests whether old composite signal is only the two headline genes",
        },
        "SPP1_CD44_anchor": {
            "genes": ["SPP1", "CD44"],
            "class": "two-gene anchor",
            "role": "biological anchor; should not be the sole manuscript signature",
        },
        "human_tubular_injury_context": {
            "genes": ["SPP1", "CD44", "HAVCR1", "LCN2", "VCAM1"],
            "class": "curated human injury context",
            "role": "existing injury-context score used in prior validation",
        },
        "human_tubular_injury_without_SPP1_CD44": {
            "genes": ["HAVCR1", "LCN2", "VCAM1"],
            "class": "headline-gene removed sensitivity",
            "role": "tests whether tubular injury context persists after removing SPP1/CD44",
        },
        "curated_fibro_inflammatory_context": {
            "genes": curated,
            "class": "curated remodeling context",
            "role": "existing fibro-inflammatory context score",
        },
        "curated_fibro_inflammatory_without_SPP1_CD44": {
            "genes": [g for g in curated if g not in {"SPP1", "CD44"}],
            "class": "headline-gene removed sensitivity",
            "role": "tests remodeling context without SPP1/CD44",
        },
        "candidate_DKD_TI_remodeling_no_SPP1_CD44": {
            "genes": dkd_ti_no_anchor,
            "class": "exploratory frozen-candidate screen",
            "role": "new DKD/TI remodeling candidate excluding the two headline genes",
        },
        "NicheNet_SPP1_receiver_target_set": {
            "genes": nichenet,
            "class": "receiver-target program",
            "role": "SPP1-prioritized predicted receiver-target layer, not direct ligand-receptor proof",
        },
        "fibrosis_ecm_reference": {
            "genes": fibrosis,
            "class": "overlap reference",
            "role": "reference module for artifact/overlap checking",
        },
        "inflammation_reference": {
            "genes": inflammation,
            "class": "overlap reference",
            "role": "reference module for artifact/overlap checking",
        },
    }
    return signatures


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    mean = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=0).replace(0, np.nan)
    z = mat.sub(mean, axis=0).div(sd, axis=0)
    return z.fillna(0.0)


def score_signatures(expr: pd.DataFrame, signatures: dict[str, dict[str, object]], dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_index = {str(g).upper(): g for g in expr.index}
    score_frames = []
    presence_rows = []
    z = zscore_rows(expr)
    for sig, meta in signatures.items():
        genes = clean_gene_list(meta["genes"])
        present_upper = [g for g in genes if g in expr_index]
        present_index = [expr_index[g] for g in present_upper]
        missing = [g for g in genes if g not in expr_index]
        if present_index:
            scores = z.loc[present_index].mean(axis=0)
        else:
            scores = pd.Series(np.nan, index=expr.columns)
        score_frames.append(
            pd.DataFrame(
                {
                    "dataset": dataset,
                    "sample_id": scores.index,
                    "signature": sig,
                    "score": scores.values,
                }
            )
        )
        presence_rows.append(
            {
                "dataset": dataset,
                "signature": sig,
                "signature_class": meta["class"],
                "role": meta["role"],
                "requested_n": len(genes),
                "present_n": len(present_upper),
                "present_fraction": len(present_upper) / len(genes) if genes else np.nan,
                "present_genes": ";".join(present_upper),
                "missing_genes": ";".join(missing),
            }
        )
    return pd.concat(score_frames, ignore_index=True), pd.DataFrame(presence_rows)


def load_gse175759_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(GSE175759_COUNTS)
    mapping = pd.read_csv(GSE175759_MAP)
    counts["ensembl"] = counts["ensembl_version"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    merged = counts.merge(mapping, on="ensembl", how="inner")
    sample_cols = [c for c in merged.columns if c.startswith("GSM")]
    gene_expr = (
        merged.assign(symbol=merged["symbol"].astype(str).str.upper())
        .groupby("symbol", as_index=True)[sample_cols]
        .mean()
    )
    gene_expr = np.log1p(gene_expr)
    meta = pd.read_csv(GSE175759_META)
    return gene_expr, meta


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))

    sample_rows = [line for line in lines if line.startswith("!Sample_")]
    row_names = [row.split("\t", 1)[0].replace("!Sample_", "") for row in sample_rows]
    counts: dict[str, int] = {}
    unique_names = []
    for name in row_names:
        counts[name] = counts.get(name, 0) + 1
        if row_names.count(name) > 1:
            unique_names.append(f"{name}_{counts[name]}")
        else:
            unique_names.append(name)
    meta = None
    for name, row in zip(unique_names, sample_rows):
        parts = [x.strip('"') for x in row.split("\t")[1:]]
        if meta is None:
            meta = pd.DataFrame({"sample_index": range(1, len(parts) + 1)})
        meta[name] = parts
    assert meta is not None
    if "geo_accession" in meta.columns:
        meta = meta.rename(columns={"geo_accession": "gsm"})
    for col, prefix, out in [
        ("characteristics_ch1_1", "tissue: ", "tissue"),
        ("characteristics_ch1_2", "tissue subregion: ", "tissue_subregion"),
        ("characteristics_ch1_3", "disease state: ", "disease_state"),
        ("characteristics_ch1_4", "individual: ", "individual"),
    ]:
        if col in meta.columns:
            meta[out] = meta[col].astype(str).str.replace(f"^{prefix}", "", regex=True)
    table_text = "\n".join(lines[begin + 1 : end])
    expr = pd.read_csv(StringIO(table_text), sep="\t")
    expr = expr.rename(columns={expr.columns[0]: "probe_id"})
    expr["probe_id"] = expr["probe_id"].astype(str).str.strip('"')
    sample_cols = [c.strip('"') for c in expr.columns[1:]]
    expr.columns = ["probe_id"] + sample_cols
    expr = expr.set_index("probe_id")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr, meta


def parse_gpl571_symbols(path: Path) -> pd.DataFrame:
    magic = path.read_bytes()[:2]
    if magic == b"\x1f\x8b":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            lines = fh.read().splitlines()
    else:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_end"))
    gpl = pd.read_csv(StringIO("\n".join(lines[begin + 1 : end])), sep="\t", dtype=str)
    symbol_col = "Gene Symbol"
    if symbol_col not in gpl.columns:
        raise ValueError("GPL571 table does not contain 'Gene Symbol'")
    rows = []
    for _, row in gpl[["ID", symbol_col]].dropna().iterrows():
        raw = str(row[symbol_col]).replace(" /// ", ";").replace("///", ";")
        for sym in raw.split(";"):
            s = sym.strip().upper()
            if s:
                rows.append((str(row["ID"]), s))
    return pd.DataFrame(rows, columns=["probe_id", "symbol"]).drop_duplicates()


def load_gse30122_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    probe_expr, meta = parse_series_matrix(GSE30122_SERIES)
    probe_map = parse_gpl571_symbols(GSE30122_GPL)
    merged = probe_map.merge(probe_expr.reset_index(), on="probe_id", how="inner")
    sample_cols = [c for c in probe_expr.columns]
    gene_expr = merged.groupby("symbol", as_index=True)[sample_cols].mean()
    return gene_expr, meta


def spearman_result(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    data = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < 3 or data["x"].nunique() < 2 or data["y"].nunique() < 2:
        return int(data.shape[0]), np.nan, np.nan
    rho, p = stats.spearmanr(data["x"], data["y"])
    return int(data.shape[0]), float(rho), float(p)


def cohen_d(control: pd.Series, dkd: pd.Series) -> float:
    control = pd.Series(control).dropna().astype(float)
    dkd = pd.Series(dkd).dropna().astype(float)
    if len(control) < 2 or len(dkd) < 2:
        return np.nan
    pooled = np.sqrt(((len(control) - 1) * control.var(ddof=1) + (len(dkd) - 1) * dkd.var(ddof=1)) / (len(control) + len(dkd) - 2))
    if pooled == 0 or not np.isfinite(pooled):
        return np.nan
    return float((dkd.mean() - control.mean()) / pooled)


def group_result(control: pd.Series, dkd: pd.Series) -> tuple[int, int, float, float, float]:
    control = pd.Series(control).dropna().astype(float)
    dkd = pd.Series(dkd).dropna().astype(float)
    if len(control) < 1 or len(dkd) < 1:
        return int(len(control)), int(len(dkd)), np.nan, np.nan, np.nan
    p = stats.mannwhitneyu(control, dkd, alternative="two-sided").pvalue
    return int(len(control)), int(len(dkd)), float(dkd.mean() - control.mean()), float(p), cohen_d(control, dkd)


def analyse_gse175759(scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    merged = scores.merge(meta, left_on="sample_id", right_on="gsm", how="left")
    rows = []
    for sig, sub in merged.groupby("signature"):
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
        disease = sub[sub["diagnosis"].astype(str).str.lower().ne("control")]
        n2, rho2, p2 = spearman_result(disease["score"], disease["egfr_ckd_epi"])
        rows.append(
            {
                "dataset": "GSE175759",
                "analysis": "eGFR_spearman_disease_only",
                "signature": sig,
                "n": n2,
                "effect": rho2,
                "p_value": p2,
                "effect_label": "spearman_rho_vs_eGFR",
                "direction_note": "negative rho means higher score with lower eGFR; excludes controls if labelled",
            }
        )
    return pd.DataFrame(rows)


def prepare_gse30122_score_meta(scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    meta["disease_binary"] = np.select(
        [
            meta.get("disease_state", pd.Series(index=meta.index, dtype=str)).astype(str).eq("diabetic kidney disease (DKD)"),
            meta.get("disease_state", pd.Series(index=meta.index, dtype=str)).astype(str).str.lower().eq("control"),
        ],
        ["DKD", "Control"],
        default=None,
    )
    meta["subregion_clean"] = np.select(
        [
            meta.get("tissue_subregion", pd.Series(index=meta.index, dtype=str)).astype(str).isin(["tubules", "tubulus"]),
            meta.get("tissue_subregion", pd.Series(index=meta.index, dtype=str)).astype(str).str.lower().eq("glomerulus"),
        ],
        ["tubular", "glomerulus"],
        default=meta.get("tissue_subregion", pd.Series(index=meta.index, dtype=str)).astype(str),
    )
    meta["is_primary_tubules"] = meta.get("tissue_subregion", pd.Series(index=meta.index, dtype=str)).astype(str).eq("tubules")
    meta["is_tubular_any"] = meta.get("tissue_subregion", pd.Series(index=meta.index, dtype=str)).astype(str).isin(["tubules", "tubulus"])
    if "individual" not in meta.columns:
        meta["individual"] = np.nan
    return scores.merge(meta, left_on="sample_id", right_on="gsm", how="left")


def analyse_gse30122(scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    merged = prepare_gse30122_score_meta(scores, meta)
    rows = []
    subsets = {
        "primary_tubules_DKD_minus_control": merged[merged["is_primary_tubules"].fillna(False).astype(bool)],
        "all_tubular_DKD_minus_control": merged[merged["is_tubular_any"].fillna(False).astype(bool)],
        "glomerulus_DKD_minus_control": merged[merged["subregion_clean"].astype(str).str.lower().eq("glomerulus")],
    }
    for analysis, df in subsets.items():
        df = df[df["disease_binary"].isin(["Control", "DKD"])]
        for sig, sub in df.groupby("signature"):
            control = sub[sub["disease_binary"].eq("Control")]["score"]
            dkd = sub[sub["disease_binary"].eq("DKD")]["score"]
            control_n, dkd_n, delta, p, d = group_result(control, dkd)
            rows.append(
                {
                    "dataset": "GSE30122",
                    "analysis": analysis,
                    "signature": sig,
                    "control_n": control_n,
                    "dkd_n": dkd_n,
                    "effect": delta,
                    "p_value": p,
                    "cohen_d": d,
                    "effect_label": "DKD_minus_control",
                    "direction_note": "positive delta means higher score in DKD",
                }
            )
    if smf is not None:
        model_df = merged[
            merged["disease_binary"].isin(["Control", "DKD"])
            & merged["subregion_clean"].isin(["tubular", "glomerulus"])
        ].copy()
        model_df["disease_code"] = model_df["disease_binary"].map({"Control": 0, "DKD": 1})
        model_df["tubular_code"] = model_df["subregion_clean"].map({"glomerulus": 0, "tubular": 1})
        for sig, sub in model_df.groupby("signature"):
            sub = sub.dropna(subset=["score", "disease_code", "tubular_code"])
            if sub.shape[0] < 8:
                continue
            try:
                fit = smf.ols("score ~ disease_code * tubular_code", data=sub).fit()
                term = "disease_code:tubular_code"
                rows.append(
                    {
                        "dataset": "GSE30122",
                        "analysis": "disease_by_tubular_compartment_interaction",
                        "signature": sig,
                        "n": int(fit.nobs),
                        "effect": float(fit.params.get(term, np.nan)),
                        "p_value": float(fit.pvalues.get(term, np.nan)),
                        "effect_label": "interaction_beta",
                        "direction_note": "positive beta means stronger DKD-control shift in tubules than glomeruli; unclustered first-pass OLS",
                    }
                )
            except Exception:
                pass
    return pd.DataFrame(rows)


def add_existing_score_comparison() -> pd.DataFrame:
    rows = []
    if GSE175759_EXISTING.exists():
        wide = pd.read_csv(GSE175759_EXISTING)
        mapping = {
            "Spp1_Cd44_tubular_immune_program": "old_NDT_composite",
            "SPP1_CD44_anchor_score": "SPP1_CD44_anchor",
            "curated_fibro_inflammatory_context": "curated_fibro_inflammatory_context",
            "human_tubular_injury_context": "human_tubular_injury_context",
        }
        for col, sig in mapping.items():
            if col in wide.columns:
                n, rho, p = spearman_result(wide[col], wide["egfr_ckd_epi"])
                rows.append(
                    {
                        "dataset": "GSE175759",
                        "analysis": "existing_pipeline_eGFR_spearman_all_samples",
                        "signature": sig,
                        "n": n,
                        "effect": rho,
                        "p_value": p,
                        "effect_label": "spearman_rho_vs_eGFR",
                        "direction_note": "existing previous pipeline score; negative rho means higher score with lower eGFR",
                    }
                )
    if GSE175759_NICHENET.exists():
        nn = pd.read_csv(GSE175759_NICHENET)
        for score_type, sub in nn.groupby("score_type"):
            n, rho, p = spearman_result(sub["score"], sub["egfr_ckd_epi"])
            rows.append(
                {
                    "dataset": "GSE175759",
                    "analysis": f"existing_pipeline_NicheNet_{score_type}_eGFR_spearman",
                    "signature": "NicheNet_SPP1_receiver_target_set",
                    "n": n,
                    "effect": rho,
                    "p_value": p,
                    "effect_label": "spearman_rho_vs_eGFR",
                    "direction_note": "existing previous pipeline score; negative rho means higher score with lower eGFR",
                }
            )
    if GSE30122_EXISTING.exists():
        long = pd.read_csv(GSE30122_EXISTING)
        mapping = {
            "Spp1_Cd44_tubular_immune_program": "old_NDT_composite",
            "SPP1_CD44_anchor_score": "SPP1_CD44_anchor",
            "curated_fibro_inflammatory_context": "curated_fibro_inflammatory_context",
            "human_tubular_injury_context": "human_tubular_injury_context",
        }
        primary = long[long["is_primary_tubules"].fillna(False).astype(bool)]
        glom = long[long["subregion_clean"].astype(str).str.lower().eq("glomerulus")]
        for analysis, df in [("existing_pipeline_primary_tubules_DKD_minus_control", primary), ("existing_pipeline_glomerulus_DKD_minus_control", glom)]:
            for old_sig, sig in mapping.items():
                sub = df[df["signature"].eq(old_sig)]
                control = sub[sub["disease_binary"].astype(str).str.lower().eq("control")]["score"]
                dkd = sub[sub["disease_binary"].astype(str).str.upper().eq("DKD")]["score"]
                control_n, dkd_n, delta, p, d = group_result(control, dkd)
                rows.append(
                    {
                        "dataset": "GSE30122",
                        "analysis": analysis,
                        "signature": sig,
                        "control_n": control_n,
                        "dkd_n": dkd_n,
                        "effect": delta,
                        "p_value": p,
                        "cohen_d": d,
                        "effect_label": "DKD_minus_control",
                        "direction_note": "existing previous pipeline score; positive delta means higher score in DKD",
                    }
                )
    if GSE30122_NICHENET.exists():
        nn = pd.read_csv(GSE30122_NICHENET)
        subsets = {
            "existing_pipeline_NicheNet_primary_tubules_DKD_minus_control": nn[nn["is_primary_tubules"].fillna(False).astype(bool)],
            "existing_pipeline_NicheNet_glomerulus_DKD_minus_control": nn[nn["tissue_subregion"].astype(str).str.lower().eq("glomerulus")],
        }
        for analysis, df in subsets.items():
            for score_type, sub in df.groupby("score_type"):
                control = sub[sub["disease_binary"].astype(str).str.lower().eq("control")]["score"]
                dkd = sub[sub["disease_binary"].astype(str).str.upper().eq("DKD")]["score"]
                control_n, dkd_n, delta, p, d = group_result(control, dkd)
                rows.append(
                    {
                        "dataset": "GSE30122",
                        "analysis": f"{analysis}_{score_type}",
                        "signature": "NicheNet_SPP1_receiver_target_set",
                        "control_n": control_n,
                        "dkd_n": dkd_n,
                        "effect": delta,
                        "p_value": p,
                        "cohen_d": d,
                        "effect_label": "DKD_minus_control",
                        "direction_note": "existing previous pipeline score; positive delta means higher score in DKD",
                    }
                )
    return pd.DataFrame(rows)


def overlap_table(signatures: dict[str, dict[str, object]]) -> pd.DataFrame:
    references = {
        "headline_SPP1_CD44": {"SPP1", "CD44"},
        "fibrosis_ecm_reference": set(clean_gene_list(signatures["fibrosis_ecm_reference"]["genes"])),
        "inflammation_reference": set(clean_gene_list(signatures["inflammation_reference"]["genes"])),
        "human_tubular_injury_context": set(clean_gene_list(signatures["human_tubular_injury_context"]["genes"])),
    }
    rows = []
    for sig, meta in signatures.items():
        genes = set(clean_gene_list(meta["genes"]))
        for ref, ref_genes in references.items():
            overlap = sorted(genes & ref_genes)
            rows.append(
                {
                    "signature": sig,
                    "reference": ref,
                    "signature_n": len(genes),
                    "reference_n": len(ref_genes),
                    "overlap_n": len(overlap),
                    "overlap_genes": ";".join(overlap),
                }
            )
    return pd.DataFrame(rows)


def rank_candidates(stats_df: pd.DataFrame, presence_df: pd.DataFrame) -> pd.DataFrame:
    pivot_rows = []
    for sig in sorted(stats_df["signature"].dropna().unique()):
        sub = stats_df[stats_df["signature"].eq(sig)]
        get = lambda analysis: sub[sub["analysis"].eq(analysis)]["effect"].dropna()
        getp = lambda analysis: sub[sub["analysis"].eq(analysis)]["p_value"].dropna()
        g175 = get("eGFR_spearman_all_samples")
        g175p = getp("eGFR_spearman_all_samples")
        pt = get("primary_tubules_DKD_minus_control")
        ptp = getp("primary_tubules_DKD_minus_control")
        gl = get("glomerulus_DKD_minus_control")
        glp = getp("glomerulus_DKD_minus_control")
        inter = get("disease_by_tubular_compartment_interaction")
        interp = getp("disease_by_tubular_compartment_interaction")
        present = presence_df[presence_df["signature"].eq(sig)]
        min_present = present["present_fraction"].min() if not present.empty else np.nan
        coherent = (
            len(g175) > 0
            and len(pt) > 0
            and float(g175.iloc[0]) < 0
            and float(pt.iloc[0]) > 0
        )
        pivot_rows.append(
            {
                "signature": sig,
                "min_present_fraction_across_GSE175759_GSE30122": min_present,
                "gse175759_rho": float(g175.iloc[0]) if len(g175) else np.nan,
                "gse175759_p": float(g175p.iloc[0]) if len(g175p) else np.nan,
                "gse30122_primary_tubule_delta": float(pt.iloc[0]) if len(pt) else np.nan,
                "gse30122_primary_tubule_p": float(ptp.iloc[0]) if len(ptp) else np.nan,
                "gse30122_glomerulus_delta": float(gl.iloc[0]) if len(gl) else np.nan,
                "gse30122_glomerulus_p": float(glp.iloc[0]) if len(glp) else np.nan,
                "gse30122_interaction_beta": float(inter.iloc[0]) if len(inter) else np.nan,
                "gse30122_interaction_p": float(interp.iloc[0]) if len(interp) else np.nan,
                "coherent_direction_first_pass": coherent,
            }
        )
    return pd.DataFrame(pivot_rows).sort_values(
        ["coherent_direction_first_pass", "gse30122_primary_tubule_delta", "gse175759_rho"],
        ascending=[False, False, True],
    )


def fmt(x: object, digits: int = 3) -> str:
    if x is None or pd.isna(x):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    try:
        return f"{float(x):.{digits}g}"
    except Exception:
        return str(x)


def write_markdown(summary_df: pd.DataFrame, stats_df: pd.DataFrame, presence_df: pd.DataFrame, overlap_df: pd.DataFrame) -> None:
    lines = [
        "# Frozen-signature selection audit",
        "",
        "Date: 2026-10-04",
        "",
        "Purpose: choose the bioinformatics signature that should be frozen before the next manuscript rebuild. This responds directly to the NDT desk-rejection concern that the previous manuscript read as descriptive rather than hypothesis-tested.",
        "",
        "Method: candidate gene sets were scored from local expression matrices using mean gene-wise z scores. GSE175759 counts were collapsed from Ensembl IDs to gene symbols and log1p-transformed before scoring. GSE30122 was parsed from the GEO series matrix and GPL571 annotation, with multiple probes averaged per gene. Existing previous-pipeline scores are retained as a comparison layer in the CSV.",
        "",
        "Interpretation rule: a main signature should show the expected direction in both human kidney-function context (higher score with lower eGFR in GSE175759) and human DKD compartment context (higher in GSE30122 DKD primary tubules). A two-gene anchor can support the story but should not become the whole manuscript claim.",
        "",
        "## First-pass candidate ranking",
        "",
        "| Candidate | Min detected fraction | GSE175759 rho vs eGFR | GSE30122 primary tubule delta | GSE30122 glomerulus delta | Interaction beta | First-pass direction |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for _, r in summary_df.iterrows():
        lines.append(
            f"| {r['signature']} | {fmt(r['min_present_fraction_across_GSE175759_GSE30122'])} | {fmt(r['gse175759_rho'])} (p={fmt(r['gse175759_p'])}) | {fmt(r['gse30122_primary_tubule_delta'])} (p={fmt(r['gse30122_primary_tubule_p'])}) | {fmt(r['gse30122_glomerulus_delta'])} (p={fmt(r['gse30122_glomerulus_p'])}) | {fmt(r['gse30122_interaction_beta'])} (p={fmt(r['gse30122_interaction_p'])}) | {'coherent' if r['coherent_direction_first_pass'] else 'not coherent'} |"
        )
    lines += [
        "",
        "## Practical conclusion",
        "",
    ]
    top = summary_df[summary_df["coherent_direction_first_pass"]].head(4)
    if not top.empty:
        lines.append("The strongest first-pass candidates are:")
        for _, r in top.iterrows():
            lines.append(
                f"- `{r['signature']}`: rho={fmt(r['gse175759_rho'])} in GSE175759 and primary-tubule DKD-control delta={fmt(r['gse30122_primary_tubule_delta'])} in GSE30122."
            )
    else:
        lines.append("No candidate met the minimal coherent-direction rule in this first-pass audit.")
    lines += [
        "",
        "The old NDT composite should not be frozen as the new main signature unless later sensitivity analyses rescue it, because the previous audit already showed weak or absent GSE30122 primary-tubule DKD elevation. The safer manuscript direction is to treat it as historical/cross-species context and move the main claim toward a DKD tubulointerstitial remodeling program, preferably one that remains positive after removing SPP1 and CD44.",
        "",
        "My working recommendation is to carry `candidate_DKD_TI_remodeling_no_SPP1_CD44` into the next sensitivity round as the primary candidate, with `human_tubular_injury_without_SPP1_CD44` as a compact positive-control injury layer. The reason is pragmatic: the three-gene injury score is statistically strong but too close to a generic injury marker set to carry the paper by itself; the 10-gene DKD/TI candidate is weaker but more biologically complete and already excludes SPP1/CD44.",
        "",
        "A key writing boundary is that GSE30122 does not prove tubule-specificity for these candidates. Several candidates also rise in glomerular samples and the first-pass disease-by-compartment interaction is not significant for the leading injury/remodeling candidates. Therefore GSE30122 should be written as human DKD compartment support, while cell-level localization must come from single-cell/snRNA analyses.",
        "",
        "## Gene detection boundary",
        "",
        "Low detection fractions make a signature unsuitable for primary claims even if one statistic looks favorable. See `06_frozen_signature_gene_presence.csv` for full present/missing gene lists.",
        "",
        "## Overlap boundary",
        "",
        "The overlap table is saved as `06_frozen_signature_overlap_audit.csv`. The next step should rerun overlap-removed tests only after one main signature is selected.",
        "",
        "## Output files",
        "",
        "- `06_frozen_signature_selection_audit.py`",
        "- `06_frozen_signature_gene_presence.csv`",
        "- `06_frozen_signature_selection_audit.csv`",
        "- `06_frozen_signature_candidate_summary.csv`",
        "- `06_frozen_signature_overlap_audit.csv`",
        "- `06_frozen_signature_scores_gse175759.csv.gz`",
        "- `06_frozen_signature_scores_gse30122.csv.gz`",
        "",
    ]
    (OUT_DIR / "06_frozen_signature_selection_audit.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    signatures = load_candidate_signatures()
    g175_expr, g175_meta = load_gse175759_expression()
    g175_scores, g175_presence = score_signatures(g175_expr, signatures, "GSE175759")
    g301_expr, g301_meta = load_gse30122_expression()
    g301_scores, g301_presence = score_signatures(g301_expr, signatures, "GSE30122")

    g175_stats = analyse_gse175759(g175_scores, g175_meta)
    g301_stats = analyse_gse30122(g301_scores, g301_meta)
    existing = add_existing_score_comparison()

    stats_df = pd.concat([g175_stats, g301_stats, existing], ignore_index=True, sort=False)
    presence_df = pd.concat([g175_presence, g301_presence], ignore_index=True, sort=False)
    overlap_df = overlap_table(signatures)
    summary_df = rank_candidates(stats_df[~stats_df["analysis"].astype(str).str.startswith("existing_pipeline")], presence_df)

    g175_scores.to_csv(OUT_DIR / "06_frozen_signature_scores_gse175759.csv.gz", index=False)
    g301_scores.to_csv(OUT_DIR / "06_frozen_signature_scores_gse30122.csv.gz", index=False)
    presence_df.to_csv(OUT_DIR / "06_frozen_signature_gene_presence.csv", index=False)
    overlap_df.to_csv(OUT_DIR / "06_frozen_signature_overlap_audit.csv", index=False)
    stats_df.to_csv(OUT_DIR / "06_frozen_signature_selection_audit.csv", index=False)
    summary_df.to_csv(OUT_DIR / "06_frozen_signature_candidate_summary.csv", index=False)
    write_markdown(summary_df, stats_df, presence_df, overlap_df)

    print(f"Wrote {OUT_DIR / '06_frozen_signature_selection_audit.md'}")
    print(f"Wrote {OUT_DIR / '06_frozen_signature_candidate_summary.csv'}")


if __name__ == "__main__":
    main()
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
