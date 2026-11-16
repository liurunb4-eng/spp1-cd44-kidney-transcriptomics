# Input-data directory

Large raw or processed public input files are intentionally not versioned in
this repository. Place them under `data/input/` using the layout documented in
`../docs/input_data_manifest.md`, or set the environment variables described
there before running `run_all.R`. Small public-data-derived identifier maps,
metadata, and component-score tables are included under `data/derived/`.

The versioned `results/` and `figures/source_data/` directories are sufficient
to audit the numerical claims and published figure without redistributing large
third-party files.
