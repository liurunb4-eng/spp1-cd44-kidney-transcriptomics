try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
source(file.path(dirname(normalizePath(.this, winslash = "/")), "00_config.R"))

effects <- readr::read_csv(file.path(OUT_DIR, "benchmark_endpoint_effects.csv"), show_col_types = FALSE)

labels <- tribble(
  ~dataset, ~display_label, ~context,
  "GSE30122_DKD_tubules", "Human DKD tubules", "disease",
  "GSE137570_progression", "Human CKD progression", "progression",
  "GSE216376_adenine", "Rat adenine injury", "disease",
  "GSE216376_UUO", "Rat UUO", "disease",
  "GSE183841_DOCA", "Rat DOCA-salt", "disease",
  "GSE183841_finerenone", "Rat finerenone", "intervention",
  "GSE183841_spironolactone", "Rat spironolactone", "intervention",
  "GSE303653_feline_cortex", "Cat natural CKD cortex", "disease",
  "GSE303653_feline_medulla", "Cat natural CKD medulla", "disease"
)

# A bounded, model-light synthesis is used because several rat contrasts have
# n=3 per group and near-zero within-group variance. A conventional pooled SMD
# would be dominated by unstable extreme values. Cliff's delta and sign
# concordance preserve direction and remain interpretable across platforms.
disease_datasets <- labels |> filter(context != "intervention") |> pull(dataset)
synthesis <- effects |>
  filter(dataset %in% disease_datasets, score_name == "common9", transform == "z", endpoint_type == "binary") |>
  left_join(labels, by = "dataset") |>
  transmute(dataset, display_label, species, context, n_reference, n_contrast,
            cliffs_delta = raw_effect, ci_low, ci_high, p_value)

if (nrow(synthesis) != length(disease_datasets)) stop("Missing disease endpoints from synthesis")

loo_medians <- vapply(synthesis$dataset, function(ds) {
  median(synthesis$cliffs_delta[synthesis$dataset != ds], na.rm = TRUE)
}, numeric(1))

synthesis_summary <- tibble(
  n_datasets = nrow(synthesis),
  n_species = n_distinct(synthesis$species),
  positive_datasets = sum(synthesis$cliffs_delta > 0),
  median_cliffs_delta = median(synthesis$cliffs_delta),
  q1 = quantile(synthesis$cliffs_delta, 0.25),
  q3 = quantile(synthesis$cliffs_delta, 0.75),
  min_effect = min(synthesis$cliffs_delta),
  max_effect = max(synthesis$cliffs_delta),
  sign_test_one_sided_p = binom.test(sum(synthesis$cliffs_delta > 0), nrow(synthesis), p = 0.5,
                                     alternative = "greater")$p.value,
  leave_one_dataset_out_median_min = min(loo_medians),
  leave_one_dataset_out_median_max = max(loo_medians)
)

species_summary <- synthesis |>
  group_by(species) |>
  summarise(n_datasets = n(), median_cliffs_delta = median(cliffs_delta),
            min_effect = min(cliffs_delta), max_effect = max(cliffs_delta), .groups = "drop")

write_csv_stable(synthesis, file.path(OUT_DIR, "cross_species_effect_synthesis.csv"))
write_csv_stable(synthesis_summary, file.path(OUT_DIR, "cross_species_effect_synthesis_summary.csv"))
write_csv_stable(species_summary, file.path(OUT_DIR, "cross_species_effect_by_species.csv"))

state_map <- effects |>
  filter(dataset %in% labels$dataset, transform == "z", endpoint_type == "binary",
         score_name %in% c("injury", "remodeling")) |>
  select(dataset, species, score_name, raw_effect, ci_low, ci_high, p_value) |>
  pivot_wider(names_from = score_name, values_from = c(raw_effect, ci_low, ci_high, p_value)) |>
  left_join(labels, by = "dataset") |>
  transmute(
    dataset, display_label, species, context,
    injury_delta = raw_effect_injury,
    remodeling_delta = raw_effect_remodeling,
    injury_ci_low = ci_low_injury,
    injury_ci_high = ci_high_injury,
    remodeling_ci_low = ci_low_remodeling,
    remodeling_ci_high = ci_high_remodeling,
    injury_p = p_value_injury,
    remodeling_p = p_value_remodeling
  ) |>
  mutate(
    quadrant = case_when(
      injury_delta >= 0 & remodeling_delta >= 0 ~ "joint induction",
      injury_delta < 0 & remodeling_delta < 0 ~ "joint attenuation",
      injury_delta >= 0 & remodeling_delta < 0 ~ "injury-predominant",
      TRUE ~ "remodeling-predominant"
    )
  )

write_csv_stable(state_map, file.path(OUT_DIR, "injury_remodeling_state_map.csv"))
message("Wrote bounded cross-species synthesis for ", nrow(synthesis),
        " disease/progression endpoints and a ", nrow(state_map), "-contrast state map.")
