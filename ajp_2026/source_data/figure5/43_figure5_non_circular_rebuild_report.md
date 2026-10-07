# Figure 5 non-circular rebuild, 2026-10-07

## Core conclusion

Selected injury-related components align with PT/PT-like cellular injury-state context after removing the mathematical circularity in the original marker-defined coordinate.

## What changed

- Panel B now uses leave-one-marker-out coordinates: the marker being tested is excluded from the PCA coordinate used for its correlation.
- Panel C uses an independent PT-state coordinate built from Injury score, PT score, LRP2 and SLC34A1, excluding SPP1, HAVCR1, LCN2 and VCAM1.
- CD44 is summarized as separate immune-context evidence and is not used to define the tubular coordinate.

## Leave-one-marker-out correlations

| dataset | marker_label | cell_level_spearman_rho | cell_level_spearman_p | sample_n | sample_mean_high_minus_low | sample_wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- |
| GSE131882 | SPP1 | 0.5298 | 0 | 6 | 0.9652 | 0.03125 |
| GSE131882 | HAVCR1 | 0.1509 | 9.49e-55 | 6 | 0.08023 | 0.03125 |
| GSE131882 | LCN2 | 0.03774 | 0.0001056 | 6 | 0.003474 | 0.125 |
| GSE131882 | VCAM1 | 0.3571 | 0 | 6 | 0.652 | 0.03125 |
| GSE195460 | SPP1 | 0.5421 | 0 | 8 | 0.9553 | 0.007812 |
| GSE195460 | HAVCR1 | 0.1897 | 1.005e-83 | 8 | 0.174 | 0.007812 |
| GSE195460 | LCN2 | 0.02896 | 0.003348 | 8 | 0.001019 | 0.25 |
| GSE195460 | VCAM1 | 0.2965 | 3.469e-207 | 8 | 0.4656 | 0.007812 |

## Old versus non-circular comparison

| dataset | marker_label | old_circular_rho | cell_level_spearman_rho |
| --- | --- | --- | --- |
| GSE131882 | SPP1 | 0.562 | 0.5298 |
| GSE131882 | HAVCR1 | 0.1752 | 0.1509 |
| GSE131882 | LCN2 | 0.04329 | 0.03774 |
| GSE131882 | VCAM1 | 0.4099 | 0.3571 |
| GSE195460 | SPP1 | 0.6131 | 0.5421 |
| GSE195460 | HAVCR1 | 0.2387 | 0.1897 |
| GSE195460 | LCN2 | 0.03445 | 0.02896 |
| GSE195460 | VCAM1 | 0.3698 | 0.2965 |

## Output files

- `43_main_figure5_singlecell_context_non_circular.svg`
- `43_main_figure5_singlecell_context_non_circular.pdf`
- `43_main_figure5_singlecell_context_non_circular.png`
- `43_main_figure5_singlecell_context_non_circular.tiff`
- `43_figure5_panelB_leave_one_marker_out_correlations.csv`
- `43_figure5_panelC_independent_coordinate_trends.csv`
- `43_non_circular_independent_coordinate_loadings.csv`
- `43_leave_one_marker_out_coordinate_loadings.csv`
- `43_old_vs_non_circular_marker_correlations.csv`
