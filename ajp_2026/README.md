# AJP 2026 human CKD injury-remodeling analysis

This directory contains the analysis-code supplement for the manuscript:

`A conserved injury-dominant tubulointerstitial remodeling signature marks kidney dysfunction in human chronic kidney disease`

It is an additive update to the repository. The historical SPP1/CD44 analysis and the frozen `v1.0.1` release remain unchanged.

## Frozen candidate

The AJP analysis evaluates a 10-gene injury-remodeling candidate comprising `HAVCR1`, `LCN2`, `VCAM1`, `CLU`, `COL1A1`, `FN1`, `ACTA2`, `TGFB1`, `TLR4`, and `NLRP3`. `SPP1` and `CD44` are used only for cellular-context analyses and are not part of the 10-gene score.

## Contents

- `scripts/analysis/`: upstream selection, matched-random-control, clinical, pan-CKD, progression, and clinicopathologic analyses.
- `scripts/figures/`: the final scripts used for Figures 1-5.
- `source_data/`: final tabular source data supporting the displayed panels.
- `docs/run_order.md`: dependency and execution map.
- `docs/input_data_manifest.md`: public input datasets and redistribution boundaries.
- `source_data/SHA256SUMS.txt`: integrity hashes for all deposited source-data files.
- `validate_deposition.py`: structural, syntax, and checksum validation.

## Reproducibility boundary

Raw sequencing matrices, large processed expression objects, and third-party atlas files are not redistributed. The deposited analysis scripts preserve the computations used in the manuscript. Set `AJP_DATA_ROOT` to the local project or public-data root and, if desired, set `AJP_WORKDIR` to the directory used for intermediate outputs and regenerated figures. `GSE137570_XLSX` can point directly to the downloaded normalized workbook. The final source-data tables are provided independently for result auditing.

## Validation

From the repository root, run:

```bash
python ajp_2026/validate_deposition.py
```

The expected terminal line is:

```text
AJP 2026 deposition validation: PASS
```

## Versioning

The intended GitHub release for this additive deposition is `v2.0.0-ajp`. The release and its archival DOI should be cited in the final manuscript Code availability statement once created.
