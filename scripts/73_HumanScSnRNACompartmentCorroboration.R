#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

paper_font <- "Arial"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}
theme_set(theme_bw(base_size = 10.5, base_family = paper_font))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  project_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  project_root <- getwd()
}
if (basename(project_root) == "scripts") {
  project_root <- dirname(project_root)
}

out_dir <- file.path(project_root, "results", "gpb_human_scsnrna_corroboration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gse131882_path <- file.path(
  project_root, "results", "human_ckd_gse131882_atlas_lite",
  "human_ckd_gse131882_atlas_lite_contrasts.csv"
)
gse195460_path <- file.path(
  project_root, "results", "human_ckd_gse195460_first_pass",
  "human_ckd_gse195460_validation_calls.csv"
)
gse211785_calls_path <- file.path(
  project_root, "results", "human_ckd_gse211785_myeloid_refinement",
  "human_ckd_gse211785_myeloid_calls.csv"
)
gse211785_scope_path <- file.path(
  project_root, "results", "human_ckd_gse211785_myeloid_refinement",
  "human_ckd_gse211785_myeloid_scope_summary.csv"
)

required <- c(gse131882_path, gse195460_path, gse211785_calls_path, gse211785_scope_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"))
}

gse131882 <- readr::read_csv(gse131882_path, show_col_types = FALSE)
gse195460 <- readr::read_csv(gse195460_path, show_col_types = FALSE)
gse211785_calls <- readr::read_csv(gse211785_calls_path, show_col_types = FALSE)
gse211785_scope <- readr::read_csv(gse211785_scope_path, show_col_types = FALSE)

pick <- function(df, key) {
  row <- df %>% filter(comparison == key)
  if (nrow(row) != 1) {
    stop("Expected one row for ", key, ", found ", nrow(row))
  }
  row
}

rows <- list()

r <- pick(gse131882, "PT_mean_SPP1_logcp10k")
rows[[length(rows) + 1]] <- tibble(
  dataset = "GSE131882",
  technology_context = "human DKD snRNA, atlas-lite",
  compartment = "PT-like tubules",
  endpoint = "SPP1 mean expression",
  gene = "SPP1",
  disease = "diabetic",
  control = "control",
  n_disease = r$n_diabetic_samples,
  n_control = r$n_control_samples,
  delta = r$delta_disease_minus_control,
  support_call = ifelse(r$delta_disease_minus_control > 0, "modest_positive", "not_positive"),
  interpretation = "Tubular SPP1 shows a modest positive disease-control shift."
)

r <- pick(gse131882, "PT_pct_tubular_injury_like")
rows[[length(rows) + 1]] <- tibble(
  dataset = "GSE131882",
  technology_context = "human DKD snRNA, atlas-lite",
  compartment = "PT-like tubules",
  endpoint = "Injury-like tubular fraction",
  gene = "injury markers",
  disease = "diabetic",
  control = "control",
  n_disease = r$n_diabetic_samples,
  n_control = r$n_control_samples,
  delta = r$delta_disease_minus_control,
  support_call = ifelse(r$delta_disease_minus_control > 0, "positive", "not_positive"),
  interpretation = "PT-like injury-context fraction is higher in disease."
)

r <- pick(gse131882, "Myeloid_mean_CD44_logcp10k")
rows[[length(rows) + 1]] <- tibble(
  dataset = "GSE131882",
  technology_context = "human DKD snRNA, atlas-lite",
  compartment = "Myeloid-like cells",
  endpoint = "CD44 mean expression",
  gene = "CD44",
  disease = "diabetic",
  control = "control",
  n_disease = r$n_diabetic_samples,
  n_control = r$n_control_samples,
  delta = r$delta_disease_minus_control,
  support_call = ifelse(r$delta_disease_minus_control > 0, "positive", "not_positive"),
  interpretation = "Myeloid CD44 shows a strong positive disease-control shift, but the marker-defined myeloid compartment is sparse."
)

r <- pick(gse195460, "PT_like_mean_SPP1_logcp10k")
rows[[length(rows) + 1]] <- tibble(
  dataset = "GSE195460",
  technology_context = "human DKD snRNA, first-pass marker-defined",
  compartment = "PT-like tubules",
  endpoint = "SPP1 mean expression",
  gene = "SPP1",
  disease = "diabetic/DN",
  control = "control",
  n_disease = r$n_diabetic_samples,
  n_control = r$n_control_samples,
  delta = r$delta_disease_minus_control,
  support_call = ifelse(r$delta_disease_minus_control > 0, "modest_positive", "not_positive"),
  interpretation = "PT-like SPP1 shows a small positive disease-control shift."
)

r <- pick(gse195460, "PT_like_pct_PT_injury_like")
rows[[length(rows) + 1]] <- tibble(
  dataset = "GSE195460",
  technology_context = "human DKD snRNA, first-pass marker-defined",
  compartment = "PT-like tubules",
  endpoint = "Injury-like tubular fraction",
  gene = "injury markers",
  disease = "diabetic/DN",
  control = "control",
  n_disease = r$n_diabetic_samples,
  n_control = r$n_control_samples,
  delta = r$delta_disease_minus_control,
  support_call = ifelse(r$delta_disease_minus_control > 0, "modest_positive", "not_positive"),
  interpretation = "PT-like injury-context fraction is mildly higher in disease."
)

for (cluster_name in c("Mac", "CD14_Mono", "CD16_Mono", "Neutrophil", "cDC")) {
  ccall <- gse211785_calls %>% filter(cluster == cluster_name)
  sc <- gse211785_scope %>% filter(cluster == cluster_name, tech_scope == "sc_rna_only")
  all <- gse211785_scope %>% filter(cluster == cluster_name, tech_scope == "all_tech")
  if (nrow(ccall) == 1 && nrow(sc) == 1 && nrow(all) == 1) {
    rows[[length(rows) + 1]] <- tibble(
      dataset = "GSE211785",
      technology_context = "human CKD/DKD kidney multi-omics atlas, myeloid refinement",
      compartment = cluster_name,
      endpoint = "CD44 mean expression",
      gene = "CD44",
      disease = "DKD",
      control = "control",
      n_disease = sc$dkd_n,
      n_control = sc$control_n,
      delta = sc$delta_cd44,
      support_call = ccall$support_call,
      interpretation = paste0(
        "SC_RNA-only CD44 delta=", sprintf("%.3f", sc$delta_cd44),
        "; all-tech delta=", sprintf("%.3f", all$delta_cd44),
        "."
      )
    )
  }
}

combined <- bind_rows(rows) %>%
  mutate(
    evidence_axis = case_when(
      gene == "SPP1" ~ "Tubular-side SPP1",
      gene == "CD44" ~ "Immune/myeloid-side CD44",
      TRUE ~ "Tubular injury context"
    ),
    support_tier = case_when(
      support_call %in% c("positive", "robust_positive") ~ "supportive",
      support_call %in% c("modest_positive", "positive_low_count") ~ "modest/bounded",
      support_call %in% c("near_neutral") ~ "boundary",
      TRUE ~ "not supportive"
    )
  )

readr::write_csv(combined, file.path(out_dir, "human_scsnrna_compartment_corroboration_summary.csv"))

plot_df <- combined %>%
  filter(endpoint %in% c("SPP1 mean expression", "CD44 mean expression", "Injury-like tubular fraction")) %>%
  mutate(
    panel = factor(evidence_axis, levels = c("Tubular-side SPP1", "Tubular injury context", "Immune/myeloid-side CD44")),
    dataset_compartment = paste(dataset, compartment, sep = "\n"),
    support_tier = factor(support_tier, levels = c("supportive", "modest/bounded", "boundary", "not supportive"))
  )

p <- ggplot(plot_df, aes(x = delta, y = reorder(dataset_compartment, delta), color = support_tier)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) +
  geom_point(size = 2.7) +
  facet_wrap(~panel, scales = "free_x", ncol = 3) +
  scale_color_manual(
    values = c(
      "supportive" = "#1B9E77",
      "modest/bounded" = "#D95F02",
      "boundary" = "#7570B3",
      "not supportive" = "#666666"
    )
  ) +
  labs(
    title = "Human sc/snRNA compartment-level corroboration of the SPP1/CD44 program",
    subtitle = "Positive delta indicates higher disease/DKD signal than control; analyses are lightweight compartment checks, not mechanistic validation",
    x = "Disease - Control delta",
    y = NULL,
    color = "Support tier"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.5)
  )
ggsave(file.path(out_dir, "Supplementary_Figure_S5_human_scsnrna_compartment_corroboration.png"), p, width = 11.2, height = 5.8, dpi = 300)
ggsave(file.path(out_dir, "Supplementary_Figure_S5_human_scsnrna_compartment_corroboration.svg"), p, width = 11.2, height = 5.8)

summary_lines <- c(
  "Human sc/snRNA compartment-level corroboration for GPB submission",
  paste0("Date: ", Sys.Date()),
  "",
  "Datasets summarized:",
  "- GSE131882: human DKD/control snRNA atlas-lite summary.",
  "- GSE195460: human DKD/DN-control snRNA first-pass marker-defined summary.",
  "- GSE211785: human kidney multi-omics atlas myeloid refinement summary.",
  "",
  "Core findings:",
  paste0(
    "- GSE131882 PT-like tubules showed a modest positive SPP1 shift (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE131882" & combined$endpoint == "SPP1 mean expression"]),
    ") and higher injury-like tubular fraction (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE131882" & combined$endpoint == "Injury-like tubular fraction"]),
    ")."
  ),
  paste0(
    "- GSE195460 PT-like tubules showed a small positive SPP1 shift (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE195460" & combined$endpoint == "SPP1 mean expression"]),
    ") and mildly higher injury-like tubular fraction (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE195460" & combined$endpoint == "Injury-like tubular fraction"]),
    ")."
  ),
  paste0(
    "- GSE131882 Myeloid-like cells showed a positive CD44 shift (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE131882" & combined$compartment == "Myeloid-like cells"]),
    "), with sparse marker-defined myeloid sample support."
  ),
  paste0(
    "- GSE211785 macrophages showed robust CD44 elevation in SC_RNA-only cells (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE211785" & combined$compartment == "Mac"]),
    "), and CD14 monocytes also showed a positive shift (delta=",
    sprintf("%.3f", combined$delta[combined$dataset == "GSE211785" & combined$compartment == "CD14_Mono"]),
    ")."
  ),
  "",
  "Recommended manuscript interpretation:",
  "- These human single-cell/single-nucleus results provide compartment-level corroboration, not full replication of the rat Spp1-Cd44 program.",
  "- The human data support a bounded pattern: modest injured-tubular/SPP1 directionality and more consistent immune/myeloid CD44 remodeling.",
  "- This layer should be placed in Supplementary Results or a short Results paragraph and should not be used as evidence of spatial contact, protein-level signaling or causal mechanism.",
  "- If included for GPB, use conservative wording such as 'human sc/snRNA compartment-level corroboration' and 'bounded cross-dataset support'."
)
writeLines(enc2utf8(summary_lines), con = file.path(out_dir, "human_scsnrna_compartment_corroboration_summary.txt"), useBytes = TRUE)

writeLines(
  enc2utf8(c(
    "gpb_human_scsnrna_corroboration completed",
    paste0("timestamp: ", Sys.time())
  )),
  con = file.path(out_dir, "status.txt"),
  useBytes = TRUE
)

message("Done. Outputs written to: ", out_dir)
