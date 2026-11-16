try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
root <- normalizePath(dirname(.this), winslash = "/", mustWork = TRUE)
rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

# These checks use only the versioned result tables and figure files. They let
# reviewers verify the release immediately without downloading the raw inputs.
steps <- c(
  "tests/test_score_build.R",
  "tests/test_benchmark_robustness.R",
  "tests/test_effect_synthesis_state_map.R",
  "tests/test_figure_outputs.R"
)

for (step in steps) {
  target <- normalizePath(file.path(root, step), winslash = "/", mustWork = TRUE)
  message("\n=== VERIFYING: ", step, " ===")
  status <- system2(rscript, c("--vanilla", shQuote(target)))
  if (!identical(status, 0L)) stop("Release verification failed at ", step)
}

message("\nPASS: versioned numerical outputs and figure deliverables satisfy all release contracts.")
