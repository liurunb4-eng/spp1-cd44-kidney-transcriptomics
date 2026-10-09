# Figure 5 marker-excluded rebuild, 2026-10-09

## Core conclusion

Selected injury-related components align with PT/PT-like cellular injury-state context after the tested marker is removed both directly and from the precomputed injury score.

## What changed

- Panels B-C use marker-excluded coordinates: the tested marker is removed from both the direct PCA features and its contribution to the injury score.
- A target-free PT-state coordinate is built from the residual injury genes, PT score, LRP2 and SLC34A1, excluding SPP1, HAVCR1, LCN2 and VCAM1.
- CD44 is summarized as separate immune-context evidence and is not used to define the tubular coordinate.

## Fully marker-excluded correlations

| dataset | marker_label | cell_level_spearman_rho | cell_level_spearman_p | sample_n | sample_mean_high_minus_low | sample_wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- |
| GSE131882 | SPP1 | 0.05365 | 3.527e-08 | 6 | 0.04325 | 0.6875 |
| GSE131882 | HAVCR1 | 0.09302 | 1.033e-21 | 6 | 0.03619 | 0.03125 |
| GSE131882 | LCN2 | 0.03233 | 0.0008964 | 6 | 0.002633 | 0.125 |
| GSE131882 | VCAM1 | 0.3571 | 0 | 6 | 0.652 | 0.03125 |
| GSE195460 | SPP1 | 0.0554 | 1.965e-08 | 8 | 0.01482 | 0.6406 |
| GSE195460 | HAVCR1 | 0.08173 | 1.125e-16 | 8 | 0.04931 | 0.007812 |
| GSE195460 | LCN2 | 0.02739 | 0.005532 | 8 | 0.001019 | 0.25 |
| GSE195460 | VCAM1 | 0.06732 | 8.787e-12 | 8 | 0.07299 | 0.02344 |

## Old versus fully marker-excluded comparison

| dataset | marker_label | old_circular_rho | cell_level_spearman_rho |
| --- | --- | --- | --- |
| GSE131882 | SPP1 | 0.562 | 0.05365 |
| GSE131882 | HAVCR1 | 0.1752 | 0.09302 |
| GSE131882 | LCN2 | 0.04329 | 0.03233 |
| GSE131882 | VCAM1 | 0.4099 | 0.3571 |
| GSE195460 | SPP1 | 0.6131 | 0.0554 |
| GSE195460 | HAVCR1 | 0.2387 | 0.08173 |
| GSE195460 | LCN2 | 0.03445 | 0.02739 |
| GSE195460 | VCAM1 | 0.3698 | 0.06732 |

## Output files

- `71_main_figure5_marker_excluded_context.svg`
- `71_main_figure5_marker_excluded_context.pdf`
- `71_main_figure5_marker_excluded_context.png`
- `71_main_figure5_marker_excluded_context.tiff`
- `71_figure5_panelB_marker_excluded_correlations.csv`
- `71_figure5_panelC_target_free_coordinate_trends.csv`
- `71_marker_excluded_sample_high_low.csv`
- `71_target_free_pt_state_bin_summary.csv`
- `71_target_free_coordinate_loadings.csv`
- `71_marker_excluded_coordinate_loadings.csv`
- `71_old_vs_marker_excluded_correlations.csv`
