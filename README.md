# SPP1/CD44 kidney transcriptomics analysis

This repository contains analysis scripts and selected supplementary tables for the manuscript:

`SPP1/CD44 Tubular-Immune Program in Kidney Injury and Renal Functional Decline`

The study re-analyzes public rat and human kidney transcriptomic datasets to nominate an SPP1/CD44-related tubular-immune program associated with kidney injury and kidney function decline.

## Scope

The repository supports the reported secondary analyses of public datasets, including:

- rat bulk kidney injury module analyses;
- rat single-cell kidney atlas mapping and cell-type-level program localization;
- focused ligand-receptor analysis of the Spp1-Cd44 candidate axis;
- permutation and multi-method robustness checks;
- state-shift and graph-based transcriptional distance analyses;
- human cohort endpoint analyses using eGFR and DKD tubule comparisons;
- supplementary NicheNet receiver-target program projection, matched random gene-set control and compartment specificity checks;
- bounded human sc/snRNA compartment-level corroboration summaries.

The analyses are designed for transcriptomic nomination and contextualization. They do not establish protein-level ligand secretion, receptor activation, spatial contact, downstream signaling, or causal function.

## Public Datasets

Raw or processed input files should be downloaded from the original public repositories before running the scripts.

| Dataset | Use in this repository |
| --- | --- |
| GSE216376 | Rat adenine nephropathy and UUO bulk RNA-seq module analysis |
| GSE183841/GSE183842 | Rat DOCA-salt and mineralocorticoid receptor antagonist bulk RNA-seq analysis |
| GSE233078 | Rat ZSF1 kidney single-cell RNA-seq discovery dataset |
| GSE175759 | Human tubulointerstitial RNA-seq cohort with eGFR metadata |
| GSE30122 | Human diabetic kidney disease tubule and glomerular microarray dataset |
| GSE131882/GSE195460/GSE211785 | Human kidney sc/snRNA or atlas-derived compartment summaries used for bounded corroboration |

Large expression matrices and raw data objects are not included in this repository.

## Repository Layout

```text
scripts/                 Analysis scripts used for the manuscript figures and tables
data/                    Placeholder for downloaded public input data
results/                 Placeholder for regenerated analysis outputs
supplementary_tables/    Supplementary table files used in the manuscript package
docs/                    Reproducibility notes and suggested run order
requirements_python.txt  Python package versions recorded from the analysis environment
requirements_r.txt       R package notes recorded from the analysis environment
```

## Suggested Run Order

The scripts retain the original project numbering. They should be run from the repository root after placing public input files under `data/`.

1. Bulk disease-context analyses: scripts `14`, `15`, `17`, `18`.
2. GSE233078 single-cell intake, atlas mapping and pseudobulk analyses: scripts `21`, `22`, `24`, `26`.
3. Focused CellChat and ligand-receptor analyses: scripts `32`, `34`, `36`.
4. DTL annotation audit and state-context analyses: scripts `37`, `38`, `40`, `42`, `60`, `69`.
5. Bulk signature projection and Spp1-Cd44 robustness checks: scripts `57`, `58`, `61`.
6. Focused LIANA validation: script `62`.
7. Human cohort analyses: scripts `63`, `67`, `68`.
8. NicheNet supplementary receiver-target analyses: scripts `68`, `70`, `71`, `72` with the longer NicheNet filenames.
9. Human sc/snRNA compartment-level corroboration summary: script `73`.

See `docs/run_order.md` for a more detailed script map.

## Environment Notes

The original analyses used a mixed Python/R workflow. Package versions recorded from the analysis environment are listed in `requirements_python.txt` and `requirements_r.txt`.

The R package list is a record of tested packages rather than a strict lockfile. For fully locked reproduction, create an `renv.lock` file from the execution environment used for final reruns.

## Data Availability

All analyzed datasets are publicly available from the NCBI Gene Expression Omnibus. This repository does not redistribute raw sequencing data, large matrices, or third-party atlas files.

## Code Availability Statement

Analysis scripts, module definitions, signature scoring code and supplementary target-set projection scripts are provided in this repository for review and reuse.


