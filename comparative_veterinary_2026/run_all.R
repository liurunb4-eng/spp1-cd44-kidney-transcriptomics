try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
root <- normalizePath(dirname(.this), winslash = "/", mustWork = TRUE)
rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

steps <- c(
  "tests/test_data_contracts.R",
  "R/01_build_scores.R",
  "tests/test_score_build.R",
  "R/02_benchmark_robustness.R",
  "tests/test_benchmark_robustness.R",
  "R/03_effect_synthesis_state_map.R",
  "tests/test_effect_synthesis_state_map.R",
  "R/04_make_competitive_figure.R",
  "tests/test_figure_outputs.R"
)

for (step in steps) {
  target <- normalizePath(file.path(root, step), winslash = "/", mustWork = TRUE)
  message("\n=== RUNNING: ", step, " ===")
  status <- system2(rscript, c("--vanilla", shQuote(target)))
  if (!identical(status, 0L)) stop("Pipeline failed at ", step, " with exit status ", status)
}

message("\nPASS: complete analysis, figure, and output-contract pipeline finished successfully.")
