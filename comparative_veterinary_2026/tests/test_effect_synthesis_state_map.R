try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
root <- normalizePath(file.path(dirname(.this), ".."), winslash = "/")
suppressPackageStartupMessages({library(readr); library(dplyr)})

s <- read_csv(file.path(root, "results", "cross_species_effect_synthesis.csv"), show_col_types = FALSE)
ss <- read_csv(file.path(root, "results", "cross_species_effect_synthesis_summary.csv"), show_col_types = FALSE)
sm <- read_csv(file.path(root, "results", "injury_remodeling_state_map.csv"), show_col_types = FALSE)

stopifnot(nrow(s) == 7L, n_distinct(s$species) == 3L)
stopifnot(all(s$cliffs_delta > 0), all(abs(s$cliffs_delta) <= 1))
stopifnot(ss$positive_datasets == 7L, ss$sign_test_one_sided_p < 0.01)
stopifnot(nrow(sm) == 9L, all(is.finite(sm$injury_delta)), all(is.finite(sm$remodeling_delta)))
stopifnot(all(c("joint induction", "joint attenuation") %in% sm$quadrant))

cat("PASS: bounded cross-species synthesis and injury-remodeling state map satisfy expected contracts.\n")
