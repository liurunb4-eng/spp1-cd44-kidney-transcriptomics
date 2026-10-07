# AJP 2026 analysis and figure run order

The filenames retain their original project step numbers. Before running upstream analyses, set `AJP_DATA_ROOT` to the directory containing the expected `data/` and `results/` trees. Outputs default to `ajp_2026/work/`; set `AJP_WORKDIR` to override that location. Set `GSE137570_XLSX` when the normalized GSE137570 workbook is stored outside the default work directory.

## Upstream analyses

1. Run `scripts/analysis/06_frozen_signature_selection_audit.py` to assemble candidate-signature scores in GSE175759 and GSE30122.
2. Run `scripts/analysis/07_candidate_sensitivity_random_control.py` for component sensitivity and 5,000 expression-matched random controls.
3. Run `scripts/analysis/12_clinical_egfr_stratification.py` for eGFR tertile and continuous-association summaries.
4. Run `scripts/analysis/13_public_outcome_screen.py` to parse external progression-associated differential-expression tables and calculate gene-level enrichment.
5. Run `scripts/analysis/16_pan_ckd_specificity_audit.py` for C-PROBE tubulointerstitial and glomerular disease-context analyses.
6. Run `scripts/analysis/35_analyze_gse137570_ckd_progression.py` for GSE137570 progression-group, GFR, and tubulointerstitial-fibrosis analyses.

## Final figures

1. Figure 1: `scripts/figures/31_main_figure1_signature_audit_reframed.py`, using outputs from steps 06 and 07.
2. Figure 2: `scripts/figures/20_main_figure2_clinical_pan_ckd_audit.py`, using outputs from steps 12 and 16.
3. Figure 3: `scripts/figures/37_main_figure3_gse137570_patient_level_support.py`, using outputs from step 35.
4. Figure 4: `scripts/figures/45_main_figure4_acoba_progression_context_bhq.py`, using outputs from step 13.
5. Figure 5: `scripts/figures/43_rebuild_figure5_non_circular.py`, using derived GSE131882/GSE195460 PT/PT-like cell tables and the immune-context summary from GSE211785.

The `source_data/figure1` through `source_data/figure5` folders contain the final tables used to audit the plotted values. They are not substitutes for the original public expression matrices.

## Path configuration

PowerShell example:

```powershell
$env:AJP_DATA_ROOT = "D:\kidney_project"
$env:AJP_WORKDIR = "D:\kidney_project\ajp_work"
$env:GSE137570_XLSX = "D:\kidney_project\inputs\GSE137570_Doyle_et_al.Normalized_read_counts.xlsx"
python ajp_2026/scripts/analysis/06_frozen_signature_selection_audit.py
```

Do not commit downloaded raw or large processed third-party data to this repository.
