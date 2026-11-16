# Comparative veterinary kidney transcriptomics (2026)

This release contains the analysis code and audited derived outputs for the
cross-species extension of a prespecified kidney injury-remodeling signature in
human, rat, and naturally occurring feline chronic kidney disease (CKD).

## What this analysis tests

1. Cross-species direction and effect-size consistency of a harmonized 9-gene
   score across human, rat, and feline kidney transcriptomic datasets.
2. Direct comparison of the composite score with individual component genes
   and prespecified injury and remodeling submodules.
3. Stability of feline results to score scaling, published quality-control
   exclusions, and leave-one-sample-out analysis.
4. Placement of disease and intervention contrasts in a two-dimensional
   injury-remodeling state map.

The 10-gene set is `HAVCR1`, `LCN2`, `VCAM1`, `CLU`, `COL1A1`, `FN1`,
`ACTA2`, `TGFB1`, `TLR4`, and `NLRP3`. The harmonized cross-species score omits
`LCN2`, which is absent from the feline expression matrix.

## Immediate release verification

The repository includes the final result tables, figure source data, and figure
files. To verify their internal contracts without downloading raw data:

```sh
Rscript verify_release.R
```

## Full regeneration

Install the packages listed in `requirements_r.txt`, acquire the public inputs
listed in `docs/input_data_manifest.md`, and then run:

```sh
Rscript run_all.R
```

By default, inputs are read from `data/input/`. A different directory can be
selected with `KIDNEY_COMPARATIVE_INPUT_DIR`. Users of the earlier project tree
may instead set `KIDNEY_PROJECT_DIR`; the feline workbook and the derived
GSE137570 component table can be supplied separately through `FELINE_CKD_XLSX`
and `GSE137570_COMPONENTS_FILE`.

The pipeline stops at the first failed data, numerical, or figure contract. It
writes only within this release directory.

## Main outputs

- `results/harmonized_sample_scores.csv`: analysis-ready scores and public
  sample identifiers for 10 endpoints across three species.
- `results/benchmark_endpoint_effects.csv`: bounded effects, uncertainty, and
  statistical tests for composites, modules, and individual genes.
- `results/cross_species_effect_synthesis*.csv`: descriptive cross-species
  synthesis and leave-one-dataset summaries.
- `results/feline_*.csv`: scaling, quality-control, and influence analyses.
- `results/injury_remodeling_state_map.csv`: two-axis state coordinates.
- `figures/source_data/`: one source-data table per figure panel.
- `figures/Figure_competitive_cross_species_validation.{svg,pdf,png}`.

## Interpretation boundary

These secondary analyses evaluate robustness, portability, and tissue-state
stratification. Disease effects are synthesized descriptively with bounded
statistics because datasets differ in species, platform, design, and sample
size. Intervention contrasts are contextual tests rather than evidence of
clinical treatment efficacy.

## Data and provenance

All sample-level inputs are from public, de-identified datasets. Large public
matrices and third-party supplementary workbooks are not redistributed here.
Small derived identifier maps and public metadata tables needed by the workflow
are versioned under `data/derived/`. Accession numbers, filenames, provenance,
and expected placement are specified in `docs/input_data_manifest.md`.
