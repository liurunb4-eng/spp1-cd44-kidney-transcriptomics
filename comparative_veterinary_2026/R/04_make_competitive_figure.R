try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
source(file.path(dirname(normalizePath(.this, winslash = "/")), "00_config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

source_dir <- file.path(FIG_DIR, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

species_cols <- c(Human = "#0072B2", Rat = "#D55E00", Cat = "#009E73")
class_cols <- c(composite = "#6A3D9A", module = "#E69F00", `single gene` = "#7A7A7A")
scale_cols <- c(global_z = "#0072B2", control_referenced_z = "#009E73", within_sample_rank = "#CC79A7")

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, colour = "black"),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    plot.margin = margin(4, 5, 4, 4),
    plot.tag = element_text(face = "bold", size = 11)
  )

# Panel A: bounded cross-species effects.
a_dat <- readr::read_csv(file.path(OUT_DIR, "cross_species_effect_synthesis.csv"), show_col_types = FALSE) |>
  mutate(display_label = factor(display_label, levels = rev(display_label)))
write_csv_stable(a_dat, file.path(source_dir, "panel_A_cross_species_effects.csv"))

p_a <- ggplot(a_dat, aes(x = cliffs_delta, y = display_label, colour = species)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.6, orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = species_cols) +
  scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  labs(x = "Cliff's delta (disease/progression vs reference)", y = NULL, colour = NULL) +
  base_theme + theme(legend.position = "none")

# Panel B: benchmark against individual genes and prespecified modules.
b_dat <- readr::read_csv(file.path(OUT_DIR, "feature_consistency_summary.csv"), show_col_types = FALSE) |>
  mutate(
    display = sub("^gene_", "", score_name),
    display = recode(display, full_available = "10-gene available", common9 = "Common 9-gene"),
    display = factor(display, levels = display[order(median_aligned_effect)])
  )
write_csv_stable(b_dat, file.path(source_dir, "panel_B_feature_benchmark.csv"))

p_b <- ggplot(b_dat, aes(x = median_aligned_effect, y = display, colour = feature_class)) +
  geom_segment(aes(x = 0, xend = median_aligned_effect, yend = display), colour = "#D0D0D0", linewidth = 0.45) +
  geom_point(size = 2.7, alpha = 0.95) +
  scale_colour_manual(values = class_cols) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  labs(x = "Median direction-aligned bounded effect", y = NULL, colour = NULL) +
  base_theme +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))

# Panel C: feline scaling, published-QC exclusion and leave-one-out stability.
c_full <- readr::read_csv(file.path(OUT_DIR, "feline_scaling_qc_sensitivity.csv"), show_col_types = FALSE)
c_loo <- readr::read_csv(file.path(OUT_DIR, "feline_leave_one_out_summary.csv"), show_col_types = FALSE)
c_dat <- c_full |>
  filter(analysis == "all_processed_samples") |>
  left_join(c_loo, by = c("dataset", "scale")) |>
  mutate(
    tissue = if_else(grepl("cortex", dataset), "Cortex", "Medulla"),
    scale_label = recode(scale,
                         global_z = "Global z",
                         control_referenced_z = "Control-ref. z",
                         within_sample_rank = "Within-sample rank"),
    scale_label = factor(scale_label, levels = c("Global z", "Control-ref. z", "Within-sample rank"))
  )
c_qc <- c_full |>
  filter(analysis == "exclude_published_QC_flags") |>
  mutate(
    tissue = if_else(grepl("cortex", dataset), "Cortex", "Medulla"),
    scale_label = recode(scale,
                         global_z = "Global z",
                         control_referenced_z = "Control-ref. z",
                         within_sample_rank = "Within-sample rank"),
    scale_label = factor(scale_label, levels = c("Global z", "Control-ref. z", "Within-sample rank"))
  )
write_csv_stable(c_dat, file.path(source_dir, "panel_C_feline_scaling_LOO.csv"))
write_csv_stable(c_qc, file.path(source_dir, "panel_C_feline_QC_exclusion.csv"))

p_c <- ggplot(c_dat, aes(x = scale_label, y = ordinal_rho, colour = scale)) +
  geom_errorbar(aes(ymin = min_rho, ymax = max_rho), width = 0.12, linewidth = 0.65) +
  geom_point(size = 2.4) +
  geom_point(data = c_qc, aes(x = scale_label, y = ordinal_rho), shape = 1, size = 3, stroke = 0.75) +
  facet_wrap(~tissue, nrow = 1) +
  scale_colour_manual(values = scale_cols, guide = "none") +
  scale_y_continuous(limits = c(0.45, 0.85), breaks = seq(0.5, 0.8, 0.1)) +
  labs(x = NULL, y = "Spearman rho with CKD stage",
       caption = "Filled: all processed samples; open: after excluding published QC flags") +
  base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        plot.caption = element_text(size = 6.5, colour = "#555555", hjust = 0.5))

# Panel D: injury-remodeling state space.
d_dat <- readr::read_csv(file.path(OUT_DIR, "injury_remodeling_state_map.csv"), show_col_types = FALSE)
write_csv_stable(d_dat, file.path(source_dir, "panel_D_injury_remodeling_state_map.csv"))

d_plot <- d_dat |>
  group_by(species, context, injury_delta, remodeling_delta) |>
  summarise(
    plot_label = if_else(n() > 1, "Rat injury models", first(display_label)),
    .groups = "drop"
  ) |>
  mutate(
    plot_label = recode(
      plot_label,
      `Cat natural CKD medulla` = "Cat CKD medulla",
      `Cat natural CKD cortex` = "Cat CKD cortex"
    )
  )

p_d <- ggplot(d_plot, aes(x = remodeling_delta, y = injury_delta, colour = species, shape = context)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#999999") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#999999") +
  geom_point(size = 2.6) +
  geom_text_repel(
    aes(label = plot_label),
    seed = 20261113,
    size = 2.05,
    colour = "#222222",
    box.padding = 0.28,
    point.padding = 0.20,
    force = 1.3,
    force_pull = 0.7,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.size = 0.3,
    segment.colour = "#9A9A9A"
  ) +
  scale_colour_manual(values = species_cols) +
  scale_shape_manual(values = c(disease = 16, progression = 17, intervention = 15)) +
  scale_x_continuous(limits = c(-1.12, 1.12), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1.12, 1.12), breaks = c(-1, 0, 1)) +
  coord_equal(clip = "off") +
  labs(x = "Remodeling-module Cliff's delta", y = "Injury-module Cliff's delta",
       colour = NULL, shape = NULL) +
  base_theme + theme(legend.position = "none")

combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.08, 1))

stem <- file.path(FIG_DIR, "Figure_competitive_cross_species_validation")
ggsave(paste0(stem, ".svg"), combined, device = svglite::svglite, width = 7.205, height = 7.35, units = "in")
ggsave(paste0(stem, ".pdf"), combined, device = cairo_pdf, width = 7.205, height = 7.35, units = "in")
ragg::agg_png(paste0(stem, ".png"), width = 7.205, height = 7.35, units = "in", res = 300, background = "white")
print(combined)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = 7.205, height = 7.35, units = "in", res = 600,
               compression = "lzw", background = "white")
print(combined)
dev.off()

message("Exported publication figure to ", FIG_DIR)
