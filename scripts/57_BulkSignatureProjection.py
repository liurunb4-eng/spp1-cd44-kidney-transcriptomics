#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import gzip
import json
import math
import re
import time
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, ttest_ind
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
OUT_DIR = PROJECT_ROOT / "results" / "bulk_signature_projection_spp1_cd44_program"
DEG_PATH = PROJECT_ROOT / "results" / "public_data_gse233078_spp1_cd44_deg_enrichment" / "spp1_cd44_axis_deg_all.csv"
GSE216376_DIR = PROJECT_ROOT / "data" / "public_datasets" / "GSE216376" / "processed"
GSE183841_DIR = PROJECT_ROOT / "data" / "public_datasets" / "GSE183841" / "processed"


MODULE_GENES = {
    "fibrosis_ecm": ["Tgfb1", "Col1a1", "Fn1", "Acta2"],
    "inflammation": ["Tlr4", "Nlrp3"],
    "vascular_microcirculation": ["Nos3", "Vegfa", "Hif1a", "Akt1"],
}


@dataclass
class Signature:
    name: str
    focus_group: str
    genes: list[str]
    construction: str


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def clean_gene(gene: object) -> str:
    return str(gene).strip()


def is_likely_gene_symbol(gene: str) -> bool:
    if not gene or "." in gene:
        return False
    if gene.upper().startswith(("LOC", "AABR", "GM", "RGD")):
        return False
    return bool(re.match(r"^[A-Za-z][A-Za-z0-9-]*$", gene))


def dedupe_keep_order(genes: list[str]) -> list[str]:
    seen = set()
    out = []
    for gene in genes:
        key = gene.upper()
        if key not in seen:
            seen.add(key)
            out.append(gene)
    return out


def build_signatures(deg_path: Path, top_n: int = 30) -> list[Signature]:
    deg = pd.read_csv(deg_path)
    deg["gene"] = deg["gene"].map(clean_gene)
    deg = deg[deg["gene"].map(is_likely_gene_symbol)].copy()
    signatures: list[Signature] = []

    for focus_group in ["DTL", "Mono_Macro"]:
        pivot = deg[deg["focus_group"] == focus_group].pivot_table(
            index="gene",
            columns="contrast",
            values=["logFC", "signed_score", "welch_p"],
            aggfunc="first",
        )
        reversible_rows = []
        for gene in pivot.index:
            try:
                disease_logfc = float(pivot.loc[gene, ("logFC", "obese_vs_lean")])
                treatment_logfc = float(pivot.loc[gene, ("logFC", "enalapril_vs_obese")])
                disease_score = float(pivot.loc[gene, ("signed_score", "obese_vs_lean")])
                treatment_score = float(pivot.loc[gene, ("signed_score", "enalapril_vs_obese")])
            except Exception:
                continue
            if disease_logfc > 0.25 and treatment_logfc < -0.25:
                reversible_rows.append(
                    {
                        "gene": gene,
                        "disease_logfc": disease_logfc,
                        "treatment_logfc": treatment_logfc,
                        "rank_score": disease_score - treatment_score,
                    }
                )
        reversible = pd.DataFrame(reversible_rows).sort_values("rank_score", ascending=False)
        if len(reversible) >= 10:
            genes = reversible["gene"].head(top_n).tolist()
            construction = "obese_vs_lean logFC>0.25 and enalapril_vs_obese logFC<-0.25; ranked by disease_score - treatment_score"
        else:
            fallback = deg[
                (deg["focus_group"] == focus_group)
                & (deg["contrast"] == "obese_vs_lean")
                & (deg["logFC"] > 0.25)
            ].sort_values("signed_score", ascending=False)
            genes = fallback["gene"].head(top_n).tolist()
            construction = "fallback obese_vs_lean logFC>0.25; ranked by signed_score because reversible overlap had <10 genes"
        signatures.append(
            Signature(
                name=f"{focus_group}_obese_up_enalapril_down_signature",
                focus_group=focus_group,
                genes=dedupe_keep_order(genes),
                construction=construction,
            )
        )

    combined = dedupe_keep_order(signatures[0].genes[:20] + signatures[1].genes[:20] + ["Spp1", "Cd44"])
    signatures.append(
        Signature(
            name="Spp1_Cd44_tubular_immune_program",
            focus_group="DTL_plus_Mono_Macro",
            genes=combined,
            construction="top reversible DTL genes + top reversible Mono/Macro genes + Spp1/Cd44 anchors",
        )
    )
    signatures.append(
        Signature(
            name="curated_fibro_inflammatory_context",
            focus_group="curated_context",
            genes=dedupe_keep_order(sum(MODULE_GENES.values(), []) + ["Spp1", "Cd44", "Havcr1", "Lcn2", "Vcam1"]),
            construction="curated fibrosis/inflammation/vascular/tubular injury context genes",
        )
    )
    return signatures


def read_gse216376_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    gene_col = df.columns[0]
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df = df.drop_duplicates(subset=gene_col).set_index(gene_col)
    df.index = df.index.str.upper()
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def read_gse183841_tpm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    gene_col = df.columns[0]
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df = df.drop_duplicates(subset=gene_col).set_index(gene_col)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def query_rat_ensembl_for_symbol(symbol: str) -> str:
    query = urlencode({"q": f"symbol:{symbol}", "species": "rat", "fields": "ensembl.gene,symbol", "size": "5"})
    url = f"https://mygene.info/v3/query?{query}"
    try:
        with urlopen(url, timeout=15) as handle:
            payload = json.loads(handle.read().decode("utf-8"))
    except Exception:
        return ""
    hits = payload.get("hits", [])
    for hit in hits:
        hit_symbol = str(hit.get("symbol", ""))
        if hit_symbol.upper() != symbol.upper():
            continue
        ens = hit.get("ensembl")
        if isinstance(ens, list) and ens:
            gene_id = ens[0].get("gene", "")
        elif isinstance(ens, dict):
            gene_id = ens.get("gene", "")
        else:
            gene_id = ""
        if gene_id:
            return str(gene_id)
    return ""


def load_or_build_rat_mapping(symbols: list[str], cache_path: Path) -> dict[str, str]:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    mapping: dict[str, str] = {}
    if cache_path.exists():
        cached = pd.read_csv(cache_path)
        mapping.update({str(row["symbol"]): str(row["rat_ensembl_id"]) for _, row in cached.iterrows() if str(row["rat_ensembl_id"]) != "nan"})

    missing = [symbol for symbol in dedupe_keep_order(symbols) if symbol not in mapping]
    for symbol in missing:
        mapping[symbol] = query_rat_ensembl_for_symbol(symbol)
        time.sleep(0.08)

    rows = [{"symbol": symbol, "rat_ensembl_id": mapping.get(symbol, "")} for symbol in sorted(mapping)]
    pd.DataFrame(rows).to_csv(cache_path, index=False, encoding="utf-8-sig")
    return mapping


def signature_score(expr: pd.DataFrame, genes: list[str], species: str, rat_mapping: dict[str, str] | None = None) -> tuple[pd.Series, list[str]]:
    if species == "rat":
        assert rat_mapping is not None
        mapped_ids = [rat_mapping.get(gene, "") for gene in genes]
        present_keys = [gene_id for gene_id in mapped_ids if gene_id and gene_id in expr.index]
        present_labels = [gene for gene in genes if rat_mapping.get(gene, "") in expr.index]
    else:
        present_keys = [gene.upper() for gene in genes if gene.upper() in expr.index]
        present_labels = [gene for gene in genes if gene.upper() in expr.index]
    if not present_keys:
        return pd.Series(0.0, index=expr.columns), []
    mat = np.log2(expr.loc[present_keys].astype(float) + 1.0)
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return z.mean(axis=0), present_labels


def welch_p(a: pd.Series, b: pd.Series) -> float:
    if len(a) < 2 or len(b) < 2:
        return 1.0
    p = ttest_ind(b, a, equal_var=False, nan_policy="omit").pvalue
    return float(p) if np.isfinite(p) else 1.0


def cohen_d(a: pd.Series, b: pd.Series) -> float:
    if len(a) < 2 or len(b) < 2:
        return 0.0
    pooled = (((len(a) - 1) * a.var(ddof=1)) + ((len(b) - 1) * b.var(ddof=1))) / max(len(a) + len(b) - 2, 1)
    if pooled <= 0 or not np.isfinite(pooled):
        return 0.0
    return float((b.mean() - a.mean()) / math.sqrt(pooled))


def score_gse216376(signatures: list[Signature]) -> tuple[pd.DataFrame, pd.DataFrame]:
    model_files = {
        "adenine": GSE216376_DIR / "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz",
        "uuo": GSE216376_DIR / "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz",
    }
    score_rows = []
    presence_rows = []
    for model_arm, path in model_files.items():
        expr = read_gse216376_matrix(path)
        column_groups = {}
        for sample_id in expr.columns:
            sid = str(sample_id).lower()
            if "control" in sid:
                column_groups[sample_id] = "control"
            elif "adenine" in sid:
                column_groups[sample_id] = "adenine"
            elif "sham" in sid:
                column_groups[sample_id] = "sham"
            elif "uuo" in sid:
                column_groups[sample_id] = "uuo"
        for signature in signatures:
            scores, present = signature_score(expr, signature.genes, species="mouse")
            presence_rows.append(
                {
                    "dataset": "GSE216376",
                    "species": "mouse",
                    "model_arm": model_arm,
                    "signature": signature.name,
                    "n_signature_genes": len(signature.genes),
                    "n_present_genes": len(present),
                    "present_genes": ";".join(present),
                    "construction": signature.construction,
                }
            )
            for sample_id, group in column_groups.items():
                if sample_id not in scores.index:
                    continue
                score_rows.append(
                    {
                        "dataset": "GSE216376",
                        "species": "mouse",
                        "model_arm": model_arm,
                        "sample_id": sample_id,
                        "group": group,
                        "signature": signature.name,
                        "score": float(scores[sample_id]),
                    }
                )
    return pd.DataFrame(score_rows), pd.DataFrame(presence_rows)


def score_gse183841(signatures: list[Signature]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(PROJECT_ROOT / "results" / "public_data_gse183841" / "sample_metadata.csv")
    expr = read_gse183841_tpm(GSE183841_DIR / "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz")
    all_symbols = dedupe_keep_order([gene for signature in signatures for gene in signature.genes])
    rat_mapping = load_or_build_rat_mapping(all_symbols, OUT_DIR / "rat_symbol_to_ensembl_mapping.csv")
    score_rows = []
    presence_rows = []
    for signature in signatures:
        scores, present = signature_score(expr, signature.genes, species="rat", rat_mapping=rat_mapping)
        presence_rows.append(
            {
                "dataset": "GSE183841",
                "species": "rat",
                "model_arm": "DOCA",
                "signature": signature.name,
                "n_signature_genes": len(signature.genes),
                "n_present_genes": len(present),
                "present_genes": ";".join(present),
                "construction": signature.construction,
            }
        )
        for _, row in metadata.iterrows():
            sample_id = row["sample_id"]
            if sample_id in scores.index:
                score_rows.append(
                    {
                        "dataset": "GSE183841",
                        "species": "rat",
                        "model_arm": "DOCA",
                        "sample_id": sample_id,
                        "group": row["normalized_group"],
                        "signature": signature.name,
                        "score": float(scores[sample_id]),
                    }
                )
    return pd.DataFrame(score_rows), pd.DataFrame(presence_rows)


def contrast_summary(score_df: pd.DataFrame) -> pd.DataFrame:
    contrasts = [
        ("GSE216376", "adenine", "control", "adenine", "disease_induction"),
        ("GSE216376", "uuo", "sham", "uuo", "disease_induction"),
        ("GSE183841", "DOCA", "control", "doca", "disease_induction"),
        ("GSE183841", "DOCA", "doca", "finerenone", "treatment_reversal"),
        ("GSE183841", "DOCA", "doca", "spironolactone", "treatment_reversal"),
    ]
    rows = []
    for dataset, model_arm, reference, contrast, role in contrasts:
        sub = score_df[(score_df["dataset"] == dataset) & (score_df["model_arm"] == model_arm)]
        for signature, sig_df in sub.groupby("signature", sort=False):
            a = sig_df[sig_df["group"] == reference]["score"]
            b = sig_df[sig_df["group"] == contrast]["score"]
            delta = float(b.mean() - a.mean()) if len(a) and len(b) else 0.0
            rows.append(
                {
                    "dataset": dataset,
                    "model_arm": model_arm,
                    "signature": signature,
                    "contrast": f"{contrast}_vs_{reference}",
                    "contrast_role": role,
                    "reference_group": reference,
                    "contrast_group": contrast,
                    "n_reference": len(a),
                    "n_contrast": len(b),
                    "reference_mean_score": float(a.mean()) if len(a) else 0.0,
                    "contrast_mean_score": float(b.mean()) if len(b) else 0.0,
                    "delta_score": delta,
                    "cohen_d": cohen_d(a, b),
                    "welch_p": welch_p(a, b),
                    "direction": "up" if delta > 0.2 else "down" if delta < -0.2 else "flat",
                }
            )
    out = pd.DataFrame(rows)
    out["bh_fdr"] = multipletests(out["welch_p"].fillna(1.0), method="fdr_bh")[1]
    return out


def module_correlation(score_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, model_arm in score_df[["dataset", "model_arm"]].drop_duplicates().itertuples(index=False):
        wide = score_df[(score_df["dataset"] == dataset) & (score_df["model_arm"] == model_arm)].pivot_table(
            index="sample_id", columns="signature", values="score", aggfunc="first"
        )
        for a in wide.columns:
            for b in wide.columns:
                if a >= b:
                    continue
                x = wide[a].dropna()
                y = wide[b].dropna()
                common = x.index.intersection(y.index)
                if len(common) < 3:
                    continue
                pearson_r, pearson_p = pearsonr(wide.loc[common, a], wide.loc[common, b])
                spearman_r, spearman_p = spearmanr(wide.loc[common, a], wide.loc[common, b])
                rows.append(
                    {
                        "dataset": dataset,
                        "model_arm": model_arm,
                        "signature_a": a,
                        "signature_b": b,
                        "n_samples": len(common),
                        "pearson_r": float(pearson_r),
                        "pearson_p": float(pearson_p),
                        "spearman_r": float(spearman_r),
                        "spearman_p": float(spearman_p),
                    }
                )
    return pd.DataFrame(rows)


def save_boxplot(score_df: pd.DataFrame, out_path: Path) -> None:
    plot_sigs = [
        "DTL_obese_up_enalapril_down_signature",
        "Mono_Macro_obese_up_enalapril_down_signature",
        "Spp1_Cd44_tubular_immune_program",
        "curated_fibro_inflammatory_context",
    ]
    panels = [
        ("GSE216376", "adenine", ["control", "adenine"]),
        ("GSE216376", "uuo", ["sham", "uuo"]),
        ("GSE183841", "DOCA", ["control", "doca", "finerenone", "spironolactone"]),
    ]
    fig, axes = plt.subplots(len(plot_sigs), len(panels), figsize=(13.5, 10), sharey="row")
    for row_idx, signature in enumerate(plot_sigs):
        for col_idx, (dataset, model_arm, groups) in enumerate(panels):
            ax = axes[row_idx, col_idx]
            sub = score_df[
                (score_df["dataset"] == dataset)
                & (score_df["model_arm"] == model_arm)
                & (score_df["signature"] == signature)
            ]
            data = [sub[sub["group"] == group]["score"].values for group in groups]
            ax.boxplot(data, tick_labels=groups, patch_artist=True, medianprops={"color": "black"})
            for i, values in enumerate(data, start=1):
                ax.scatter(np.full(len(values), i), values, s=22, color="#2f4858", alpha=0.8, zorder=3)
            if row_idx == 0:
                ax.set_title(f"{dataset}\n{model_arm}", fontsize=10)
            if col_idx == 0:
                ax.set_ylabel(signature.replace("_", "\n"), fontsize=8)
            ax.tick_params(axis="x", labelrotation=30, labelsize=8)
            ax.axhline(0, lw=0.8, color="#999999", alpha=0.7)
    fig.suptitle("Bulk projection of GSE233078-derived tubular-immune signatures", fontsize=14, y=0.995)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_delta_heatmap(summary: pd.DataFrame, out_path: Path) -> None:
    plot_df = summary.pivot_table(index="signature", columns="contrast", values="delta_score", aggfunc="first").fillna(0.0)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    vmax = max(0.5, float(np.nanmax(np.abs(plot_df.values))))
    im = ax.imshow(plot_df.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(plot_df.shape[1]))
    ax.set_xticklabels(plot_df.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(plot_df.shape[0]))
    ax.set_yticklabels([idx.replace("_", " ") for idx in plot_df.index], fontsize=8)
    for i in range(plot_df.shape[0]):
        for j in range(plot_df.shape[1]):
            ax.text(j, i, f"{plot_df.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Signature score shifts across public bulk cohorts")
    fig.colorbar(im, ax=ax, shrink=0.82, label="Mean score delta")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signatures = build_signatures(DEG_PATH)
    signature_rows = [
        {
            "signature": sig.name,
            "focus_group": sig.focus_group,
            "n_genes": len(sig.genes),
            "genes": ";".join(sig.genes),
            "construction": sig.construction,
        }
        for sig in signatures
    ]
    pd.DataFrame(signature_rows).to_csv(OUT_DIR / "signature_definitions.csv", index=False, encoding="utf-8-sig")

    gse216_scores, gse216_presence = score_gse216376(signatures)
    gse183_scores, gse183_presence = score_gse183841(signatures)
    score_df = pd.concat([gse216_scores, gse183_scores], ignore_index=True)
    presence_df = pd.concat([gse216_presence, gse183_presence], ignore_index=True)
    summary = contrast_summary(score_df)
    corrs = module_correlation(score_df)

    score_df.to_csv(OUT_DIR / "bulk_signature_scores_long.csv", index=False, encoding="utf-8-sig")
    presence_df.to_csv(OUT_DIR / "bulk_signature_gene_presence.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT_DIR / "bulk_signature_contrast_summary.csv", index=False, encoding="utf-8-sig")
    corrs.to_csv(OUT_DIR / "bulk_signature_score_correlations.csv", index=False, encoding="utf-8-sig")
    save_boxplot(score_df, OUT_DIR / "FigureBL_bulk_signature_projection_boxplots.png")
    save_delta_heatmap(summary, OUT_DIR / "FigureBM_bulk_signature_projection_delta_heatmap.png")

    key = summary[summary["signature"].isin(["Spp1_Cd44_tubular_immune_program", "DTL_obese_up_enalapril_down_signature", "Mono_Macro_obese_up_enalapril_down_signature"])]
    positive_disease = key[key["contrast_role"] == "disease_induction"].sort_values(["signature", "contrast"])
    treatment = key[key["contrast_role"] == "treatment_reversal"].sort_values(["signature", "contrast"])
    lines = [
        "Bulk signature projection of GSE233078-derived Spp1/Cd44-related programs",
        "",
        "Purpose:",
        "This analysis projects mouse single-cell-derived DTL and Mono/Macro disease signatures into independent public bulk cohorts.",
        "It is intended as cross-platform support, not as deconvolution or wet-lab validation.",
        "",
        "Signature definitions:",
    ]
    for sig in signatures:
        lines.append(f"- {sig.name}: n={len(sig.genes)}; {sig.construction}")
    lines.extend(["", "Key disease-induction readouts:"])
    for _, row in positive_disease.iterrows():
        lines.append(
            f"- {row['dataset']} {row['contrast']} | {row['signature']}: delta={row['delta_score']:.3f}, "
            f"d={row['cohen_d']:.3f}, p={row['welch_p']:.4g}, FDR={row['bh_fdr']:.4g}, direction={row['direction']}"
        )
    lines.extend(["", "Key treatment/reversal readouts:"])
    for _, row in treatment.iterrows():
        lines.append(
            f"- {row['dataset']} {row['contrast']} | {row['signature']}: delta={row['delta_score']:.3f}, "
            f"d={row['cohen_d']:.3f}, p={row['welch_p']:.4g}, FDR={row['bh_fdr']:.4g}, direction={row['direction']}"
        )
    lines.extend(
        [
            "",
            "Interpretation boundary:",
            "1. These scores are rank/z-score-style bulk projections of single-cell-derived signatures.",
            "2. They do not prove cell proportions or ligand-receptor communication in bulk tissue.",
            "3. A disease-up and treatment-down pattern can be written as cross-platform directional support for the pathological program.",
            "4. If a signature fails in a cohort, it should be reported as model/context sensitivity rather than hidden.",
        ]
    )
    write_text(OUT_DIR / "status.txt", "completed\n")


if __name__ == "__main__":
    main()


