# -*- coding: utf-8 -*-
"""
Step 13: public longitudinal/progression outcome screen and Acoba 2025 progression-DE alignment.
"""
from pathlib import Path
import os
import re, gzip, json, math
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[3]
BASE = Path(os.environ.get("AJP_DATA_ROOT", REPO_ROOT))
OUT = Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
OUT.mkdir(parents=True, exist_ok=True)
WORK = OUT / "step13_public_outcome_screen"
SUPP = WORK / "Acoba_2025_BMCNephrol_supp" / "12882_2025_4364_MOESM2_ESM.xlsx"

CANDIDATE = ["HAVCR1","LCN2","VCAM1","CLU","COL1A1","FN1","ACTA2","TGFB1","TLR4","NLRP3"]
INJURY = ["HAVCR1","LCN2","VCAM1","CLU"]
OLD_ANCHOR = ["SPP1","CD44"]

SHEET_MAP = {
    "steep eGFR slope (Glom)": ("steep_eGFR_slope", "glomeruli"),
    "steep eGFR slope (Tub)": ("steep_eGFR_slope", "tubulointerstitium"),
    "CKD stage advancement (Glom)": ("CKD_stage_advancement", "glomeruli"),
    "CKD stage advancement (Tub)": ("CKD_stage_advancement", "tubulointerstitium"),
    "UACR increase (Glom)": ("UACR_increase", "glomeruli"),
    "UACR increase (Tub)": ("UACR_increase", "tubulointerstitium"),
    "NRA (Glom)": ("NRA", "glomeruli"),
    "NRA (Tub)": ("NRA", "tubulointerstitium"),
    "Composite outcome (Glom)": ("composite_outcome", "glomeruli"),
    "Composite outcome (Tub)": ("composite_outcome", "tubulointerstitium"),
}

# Parse all DE sheets
xl = pd.ExcelFile(SUPP)
all_rows = []
for sheet, (endpoint, compartment) in SHEET_MAP.items():
    df = xl.parse(sheet)
    df.columns = [str(c).strip() for c in df.columns]
    if "symbol" not in df.columns:
        continue
    df = df.dropna(subset=["symbol"]).copy()
    df["symbol"] = df["symbol"].astype(str).str.upper().str.strip()
    df["endpoint"] = endpoint
    df["compartment"] = compartment
    df["sheet"] = sheet
    all_rows.append(df[["endpoint","compartment","sheet","ensembl_id","symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]])
all_de = pd.concat(all_rows, ignore_index=True)
all_de.to_csv(OUT / "13_acoba2025_all_progression_de_tables.csv.gz", index=False, compression="gzip")

# Signature gene extraction
gene_rows=[]
for sig_name, genes in [("candidate_10_gene", CANDIDATE), ("injury_core", INJURY), ("SPP1_CD44_anchor", OLD_ANCHOR)]:
    genes_upper=[g.upper() for g in genes]
    sub=all_de[all_de["symbol"].isin(genes_upper)].copy()
    present=set(sub["symbol"])
    # add rows for missing per endpoint/compartment for audit
    for endpoint in sorted(all_de.endpoint.unique()):
        for comp in sorted(all_de.compartment.unique()):
            chunk=sub[(sub.endpoint==endpoint)&(sub.compartment==comp)]
            for _,r in chunk.iterrows():
                gene_rows.append({
                    "signature": sig_name,
                    "gene": r.symbol,
                    "endpoint": endpoint,
                    "compartment": comp,
                    "log2FoldChange_rapid_vs_nonrapid": r.log2FoldChange,
                    "pvalue": r.pvalue,
                    "padj": r.padj,
                    "direction": "up_in_rapid_progression" if r.log2FoldChange>0 else "down_in_rapid_progression",
                    "nominal_p_lt_0_05": bool(r.pvalue < 0.05),
                    "fdr_lt_0_10": bool(r.padj < 0.10) if pd.notna(r.padj) else False,
                    "fdr_lt_0_05": bool(r.padj < 0.05) if pd.notna(r.padj) else False,
                })
            for g in genes_upper:
                if g not in set(chunk.symbol):
                    # if absent from a sheet, record absence; should be rare if all genes retained
                    gene_rows.append({"signature": sig_name, "gene": g, "endpoint": endpoint, "compartment": comp, "log2FoldChange_rapid_vs_nonrapid": np.nan, "pvalue": np.nan, "padj": np.nan, "direction": "not_found", "nominal_p_lt_0_05": False, "fdr_lt_0_10": False, "fdr_lt_0_05": False})

gene_df=pd.DataFrame(gene_rows).drop_duplicates()
gene_df.to_csv(OUT / "13_acoba2025_signature_gene_progression_alignment.csv", index=False, encoding="utf-8-sig")

# Endpoint summary by signature and compartment
summary=[]
for (sig, endpoint, comp), g in gene_df.groupby(["signature","endpoint","compartment"]):
    gg=g[g.direction!="not_found"].copy()
    if gg.empty:
        continue
    n=len(gg)
    summary.append({
        "signature": sig,
        "endpoint": endpoint,
        "compartment": comp,
        "n_genes_found": n,
        "n_up_in_rapid": int((gg.log2FoldChange_rapid_vs_nonrapid>0).sum()),
        "n_down_in_rapid": int((gg.log2FoldChange_rapid_vs_nonrapid<0).sum()),
        "median_log2FC": float(np.nanmedian(gg.log2FoldChange_rapid_vs_nonrapid)),
        "mean_log2FC": float(np.nanmean(gg.log2FoldChange_rapid_vs_nonrapid)),
        "n_nominal_p_lt_0_05": int(gg.nominal_p_lt_0_05.sum()),
        "n_fdr_lt_0_10": int(gg.fdr_lt_0_10.sum()),
        "n_fdr_lt_0_05": int(gg.fdr_lt_0_05.sum()),
        "nominal_genes": ";".join(gg.loc[gg.nominal_p_lt_0_05,"gene"].tolist()),
        "fdr05_genes": ";".join(gg.loc[gg.fdr_lt_0_05,"gene"].tolist()),
    })
summary_df=pd.DataFrame(summary)
summary_df.to_csv(OUT / "13_acoba2025_signature_progression_summary.csv", index=False, encoding="utf-8-sig")

# ORA: enrichment among nominal/up rapid DE genes for candidate genes.
# Universe = symbols in each endpoint/compartment table. Test overlap with nominal p<0.05 and logFC>0.
ora=[]
for (endpoint, comp), d in all_de.groupby(["endpoint","compartment"]):
    universe=set(d.symbol.dropna().astype(str))
    for sig, genes in [("candidate_10_gene", CANDIDATE), ("injury_core", INJURY), ("SPP1_CD44_anchor", OLD_ANCHOR)]:
        gs=set([g.upper() for g in genes]) & universe
        for label, mask in [
            ("nominal_up_in_rapid", (d.pvalue<0.05) & (d.log2FoldChange>0)),
            ("fdr05_up_in_rapid", (d.padj<0.05) & (d.log2FoldChange>0)),
            ("nominal_any_direction", (d.pvalue<0.05)),
            ("fdr05_any_direction", (d.padj<0.05)),
        ]:
            hits=set(d.loc[mask,"symbol"])
            a=len(gs & hits)
            b=len(gs - hits)
            c=len(hits - gs)
            dd=len(universe - gs - hits)
            odds,p=fisher_exact([[a,b],[c,dd]], alternative="greater") if len(gs)>0 and len(hits)>0 else (np.nan,np.nan)
            ora.append({
                "signature":sig,"endpoint":endpoint,"compartment":comp,"set_tested":label,
                "n_signature_present":len(gs),"n_hit_genes":len(hits),"overlap_n":a,
                "overlap_genes":";".join(sorted(gs&hits)),"odds_ratio":odds,"fisher_p":p
            })
ora_df=pd.DataFrame(ora)
ora_df.to_csv(OUT / "13_acoba2025_signature_progression_enrichment.csv", index=False, encoding="utf-8-sig")

# Parse GSE180 matrices metadata compactly for inventory
def parse_geo_meta(matrix_path):
    records={}
    sample_ids=[]
    with gzip.open(matrix_path,'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if not line.startswith('!'):
                continue
            parts=line.rstrip('\n').split('\t')
            key=parts[0]
            vals=[x.strip('"') for x in parts[1:]]
            if key == '!Sample_geo_accession':
                sample_ids=vals
                for sid in sample_ids:
                    records.setdefault(sid,{})['geo_accession']=sid
            elif key == '!Sample_title' and sample_ids:
                for sid,v in zip(sample_ids,vals): records.setdefault(sid,{})['title']=v
            elif key == '!Sample_source_name_ch1' and sample_ids:
                for sid,v in zip(sample_ids,vals): records.setdefault(sid,{})['source_name']=v
            elif key == '!Sample_characteristics_ch1' and sample_ids:
                for sid,v in zip(sample_ids,vals):
                    if ':' in v:
                        k,val=v.split(':',1); k=k.strip().lower().replace(' ','_'); val=val.strip()
                    else:
                        k='characteristics'; val=v
                    # avoid overwrite
                    old=records.setdefault(sid,{}).get(k)
                    records[sid][k]=val if old is None else str(old)+'; '+val
    return pd.DataFrame(records.values())

inv=[]
for acc in ["GSE180393","GSE180394","GSE180395"]:
    p=WORK / f"{acc}_series_matrix.txt.gz"
    if p.exists():
        meta=parse_geo_meta(p)
        meta.to_csv(WORK / f"{acc}_sample_metadata_parsed.csv", index=False, encoding="utf-8-sig")
        fields=';'.join(meta.columns)
        has_outcome=bool(re.search(r'(egfr|gfr|esrd|dialysis|follow|outcome|event|slope|survival|progress)', fields, re.I))
        sample_group_counts=''
        if 'sample_group' in meta.columns:
            sample_group_counts='; '.join([f"{k}:{v}" for k,v in meta['sample_group'].value_counts().head(12).items()])
        tissue_counts=''
        if 'tissue' in meta.columns:
            tissue_counts='; '.join([f"{k}:{v}" for k,v in meta['tissue'].value_counts().items()])
        inv.append({
            'resource':acc,
            'n_samples':len(meta),
            'fields':fields,
            'sample_group_counts_head':sample_group_counts,
            'tissue_counts':tissue_counts,
            'outcome_fields_in_geo_matrix':has_outcome,
            'current_use':'expression possible, but no individual longitudinal outcomes in GEO matrix' if not has_outcome else 'inspect outcome fields',
        })
inv_df=pd.DataFrame(inv)
inv_df.to_csv(OUT / "13_public_longitudinal_dataset_inventory.csv", index=False, encoding="utf-8-sig")

# Dataset recommendation table
recommendations = pd.DataFrame([
    {
        "resource":"Acoba et al. 2025 / KaroKidney RNA-seq-DN",
        "fit":"highest",
        "why":"Human DKD kidney biopsy RNA-seq with 5-year follow-up/progression endpoints reported; supplementary tables contain endpoint-specific glomerular and tubulointerstitial DE results.",
        "available_now":"Supplementary gene-level DE tables downloaded; individual expression matrix not yet retrieved because KaroKidney site timed out in this session.",
        "analysis_done_now":"Signature-gene and enrichment alignment against rapid-progression DE tables.",
        "can_do_sample_level_survival_now":"No",
        "next_action":"Try KaroKidney data download again or contact/download mirror; if expression matrix has same sample IDs and progression labels, compute sample-level signature score vs rapid progression/eGFR slope.",
    },
    {
        "resource":"GSE180393/GSE180394/GSE180395 C-PROBE",
        "fit":"medium",
        "why":"Large human microdissected kidney biopsy transcriptome; GEO summaries mention association with disease outcomes.",
        "available_now":"Series matrices downloaded and parsed from NCBI GEO.",
        "analysis_done_now":"Metadata audit; sample group/tissue fields parsed.",
        "can_do_sample_level_survival_now":"No",
        "next_action":"Need external C-PROBE clinical outcome table; GEO matrix does not expose eGFR slope/time-to-event fields.",
    },
    {
        "resource":"GSE175759",
        "fit":"baseline clinical anchor only",
        "why":"Human tubulointerstitial RNA-seq with baseline eGFR; already analyzed in Step 12.",
        "available_now":"Expression and baseline eGFR available locally.",
        "analysis_done_now":"Cross-sectional association and eGFR-tertile analysis.",
        "can_do_sample_level_survival_now":"No",
        "next_action":"Keep as cross-sectional clinical relevance layer; do not write progression prediction.",
    },
])
recommendations.to_csv(OUT / "13_public_outcome_dataset_recommendations.csv", index=False, encoding="utf-8-sig")

# Markdown report
cand_summ = summary_df[(summary_df.signature=='candidate_10_gene')].copy()
# choose key lines: tubulointerstitium eGFR slope, CKD stage, composite
key = cand_summ[(cand_summ.compartment=='tubulointerstitium') & (cand_summ.endpoint.isin(['steep_eGFR_slope','CKD_stage_advancement','composite_outcome']))]
lines=[]
lines.append("# Step 13. Public longitudinal/progression outcome screen\n\n")
lines.append("Date: 2026-10-05\n\n")
lines.append("## Question\n\nCan we add a true longitudinal/prognostic bioinformatics layer after the NDT desk rejection, without overclaiming from cross-sectional eGFR data?\n\n")
lines.append("## Main conclusion\n\n")
lines.append("A true sample-level survival/eGFR-slope analysis is not yet possible from the files currently available locally. The best immediately usable progression resource is Acoba et al. 2025 / KaroKidney, but the accessible material in this session is gene-level progression differential-expression tables rather than individual expression plus outcome metadata. Therefore the safe next evidence layer is progression-DE alignment, not Cox/Kaplan-Meier analysis.\n\n")
lines.append("## Acoba 2025 progression-DE alignment\n\n")
lines.append("Downloaded and parsed the BMC Nephrology supplementary differential-expression workbook for rapid DKD progression endpoints. The workbook includes glomerular and tubulointerstitial results for steep eGFR slope, CKD stage advancement, UACR increase, NRA, and composite outcome.\n\n")
lines.append("For the current 10-gene candidate, key tubulointerstitial summaries are:\n\n")
for r in key.itertuples(index=False):
    lines.append(f"- {r.endpoint}: {r.n_genes_found} genes found; {r.n_up_in_rapid} up and {r.n_down_in_rapid} down in rapid progression; median log2FC = {r.median_log2FC:.3f}; nominal p<0.05 genes = {r.nominal_genes or 'none'}; FDR<0.05 genes = {r.fdr05_genes or 'none'}.\n")
lines.append("\nThis is useful as an external progression-context check, but it is not equivalent to validating a sample-level risk score.\n\n")
lines.append("## Candidate datasets screened\n\n")
for r in recommendations.itertuples(index=False):
    lines.append(f"### {r.resource}\n\n")
    lines.append(f"Fit: {r.fit}.\n\n")
    lines.append(f"Why: {r.why}\n\n")
    lines.append(f"Available now: {r.available_now}\n\n")
    lines.append(f"Current conclusion: sample-level survival/eGFR-slope analysis now? {r.can_do_sample_level_survival_now}.\n\n")
    lines.append(f"Next action: {r.next_action}\n\n")
lines.append("## Manuscript-safe wording\n\n")
lines.append("Use: `External DKD progression-associated transcriptomic data provided gene-level support for parts of the injury-remodeling candidate, particularly as a progression-context check.`\n\n")
lines.append("Avoid: `the signature predicts DKD progression`, `Cox analysis confirmed risk stratification`, or `validated as an ESRD survival marker`, unless individual-level expression and longitudinal outcomes are later retrieved and analyzed.\n\n")
lines.append("## Outputs\n\n")
for fn in [
    "13_acoba2025_signature_gene_progression_alignment.csv",
    "13_acoba2025_signature_progression_summary.csv",
    "13_acoba2025_signature_progression_enrichment.csv",
    "13_public_longitudinal_dataset_inventory.csv",
    "13_public_outcome_dataset_recommendations.csv",
    "13_public_outcome_screen.md",
]:
    lines.append(f"- `{fn}`\n")
(OUT / "13_public_outcome_screen.md").write_text("".join(lines), encoding="utf-8")

print('DONE')
print(OUT / '13_public_outcome_screen.md')
print('\nKEY SUMMARY')
print(key.to_string(index=False))
print('\nCANDIDATE INVENTORY')
print(recommendations[['resource','fit','can_do_sample_level_survival_now']].to_string(index=False))
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
