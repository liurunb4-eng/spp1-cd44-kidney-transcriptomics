# Input-data manifest

The full pipeline uses public, de-identified transcriptomic data. Large matrices
and third-party supplementary files are not redistributed. The table below
records the exact filenames expected under `data/input/`.

| Resource | Expected relative path | Role | Source |
| --- | --- | --- | --- |
| GSE175759 | `GSE175759/gse175759_merged_counts_by_ensembl.csv.gz` | Human tubulointerstitial expression | GEO GSE175759; derived merged count matrix |
| GSE175759 | `data/derived/GSE175759/gse175759_ensembl_to_symbol_map.csv` | Ensembl-symbol mapping | Included; derived from the GSE175759 feature annotation |
| GSE175759 | `data/derived/GSE175759/gse175759_metadata_clean.csv` | eGFR and technical-QC metadata | Included; derived from public GEO metadata |
| GSE30122 | `GSE30122/GSE30122_series_matrix.txt.gz` | Human DKD tubule/control expression | GEO GSE30122 |
| GPL571 | `GSE30122/GPL571.soft.gz` | GSE30122 probe annotation | GEO platform GPL571 |
| GSE137570 | `data/derived/GSE137570/GSE137570_cohort2_component_gene_scores_long.csv` | Public-cohort progression component scores | Included; derived from the public GSE137570 cohort-2 workbook using the prior analysis release |
| GSE216376 | `GSE216376/GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz` | Rat adenine nephropathy | GEO GSE216376 processed file |
| GSE216376 | `GSE216376/GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz` | Rat UUO | GEO GSE216376 processed file |
| GSE183841 | `GSE183841/GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz` | Rat DOCA-salt and antagonist arms | GEO GSE183841 processed file |
| GSE183841 | `data/derived/GSE183841/sample_metadata.csv` | Treatment-group metadata | Included; derived from public GSE183841 metadata |
| Rat annotation | `data/derived/mappings/rat_symbol_to_ensembl_mapping.csv` | Rat gene identifier mapping | Included; derived from the public feature annotation |
| GSE303653 / Li et al. 2025 | `GSE303653/42003_2025_9164_MOESM3_ESM.xlsx` | Naturally occurring feline CKD cortex and medulla | Article supplementary workbook; DOI `10.1038/s42003-025-09164-8` |

## Configuration options

- `KIDNEY_COMPARATIVE_INPUT_DIR`: replaces the default `data/input/` root.
- `KIDNEY_PROJECT_DIR`: points to the legacy analysis project root and enables
  the earlier directory layout for human and rat inputs.
- `GSE137570_COMPONENTS_FILE`: explicit path to the derived GSE137570 table.
- `FELINE_CKD_XLSX`: explicit path to the feline supplementary workbook.

Example for PowerShell:

```powershell
$env:KIDNEY_COMPARATIVE_INPUT_DIR = "D:\public_kidney_inputs"
Rscript run_all.R
```

Example for a legacy project tree:

```powershell
$env:KIDNEY_PROJECT_DIR = "D:\kidney_project"
$env:FELINE_CKD_XLSX = "D:\public_kidney_inputs\GSE303653\42003_2025_9164_MOESM3_ESM.xlsx"
Rscript run_all.R
```
