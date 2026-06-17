suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  project_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  project_root <- getwd()
}

out_dir <- file.path(project_root, "results", "gse233078_nichenet_ligand_target")
deg_dir <- file.path(project_root, "results", "public_data_gse233078_spp1_cd44_deg_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lt_path <- file.path(out_dir, "ligand_target_matrix.rds")
lr_path <- file.path(out_dir, "lr_network.rds")
counts_path <- file.path(deg_dir, "pseudobulk_counts_focus_groups.csv")
meta_path <- file.path(deg_dir, "pseudobulk_group_metadata.csv")
deg_path <- file.path(deg_dir, "spp1_cd44_axis_deg_all.csv")

required <- c(lt_path, lr_path, counts_path, meta_path, deg_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
}

to_upper_gene <- function(x) toupper(trimws(as.character(x)))

make_logcpm <- function(counts) {
  gene <- to_upper_gene(counts$gene)
  numeric <- as.matrix(counts[, setdiff(colnames(counts), "gene"), drop = FALSE])
  storage.mode(numeric) <- "numeric"
  summed <- rowsum(numeric, group = gene, reorder = FALSE)
  lib <- colSums(summed)
  log2(t(t(summed) / lib * 1e6) + 0.5)
}

average_precision <- function(scores, positives) {
  ok <- is.finite(scores) & !is.na(positives)
  scores <- scores[ok]
  positives <- positives[ok]
  n_pos <- sum(positives)
  if (n_pos == 0 || n_pos == length(positives)) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  pos_ord <- positives[ord]
  precision <- cumsum(pos_ord) / seq_along(pos_ord)
  mean(precision[pos_ord])
}

roc_auc <- function(scores, positives) {
  ok <- is.finite(scores) & !is.na(positives)
  scores <- scores[ok]
  positives <- positives[ok]
  n_pos <- sum(positives)
  n_neg <- length(positives) - n_pos
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[positives]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

write_utf8 <- function(path, lines) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

ligand_target_matrix <- readRDS(lt_path)
lr_network <- readRDS(lr_path)
counts <- read.csv(counts_path, check.names = FALSE)
metadata <- read.csv(meta_path, check.names = FALSE)
deg <- read.csv(deg_path, check.names = FALSE)

rownames(ligand_target_matrix) <- to_upper_gene(rownames(ligand_target_matrix))
colnames(ligand_target_matrix) <- to_upper_gene(colnames(ligand_target_matrix))
lr_network$from <- to_upper_gene(lr_network$from)
lr_network$to <- to_upper_gene(lr_network$to)
deg$gene_upper <- to_upper_gene(deg$gene)

logcpm <- make_logcpm(counts)
common_targets <- intersect(rownames(ligand_target_matrix), rownames(logcpm))

mean_by <- function(focus_group, condition) {
  ids <- metadata$group_id[metadata$focus_group == focus_group & metadata$condition == condition]
  ids <- intersect(ids, colnames(logcpm))
  if (length(ids) == 0) stop("No pseudobulk columns for ", focus_group, " / ", condition)
  rowMeans(logcpm[, ids, drop = FALSE])
}

sender_obese <- mean_by("DTL", "obese")
receiver_obese <- mean_by("Mono_Macro", "obese")
receiver_lean <- mean_by("Mono_Macro", "lean")

expressed_sender <- names(sender_obese)[sender_obese > 1]
expressed_receiver <- names(receiver_obese)[receiver_obese > 1]

deg_receiver <- deg[deg$focus_group == "Mono_Macro" & deg$contrast == "obese_vs_lean", ]
genes_interest <- deg_receiver$gene_upper[
  deg_receiver$logFC > 0 &
    deg_receiver$welch_p < 0.1 &
    deg_receiver$gene_upper %in% common_targets
]
genes_interest <- unique(genes_interest)
if (length(genes_interest) < 30) {
  fallback <- deg_receiver[deg_receiver$logFC > 0 & deg_receiver$gene_upper %in% common_targets, ]
  fallback <- fallback[order(-fallback$signed_score), ]
  genes_interest <- unique(head(fallback$gene_upper, 250))
}

background <- intersect(expressed_receiver, common_targets)
genes_interest <- intersect(genes_interest, background)

candidate_ligands <- unique(lr_network$from[
  lr_network$from %in% expressed_sender &
    lr_network$from %in% colnames(ligand_target_matrix) &
    lr_network$to %in% expressed_receiver
])
candidate_ligands <- unique(c(candidate_ligands, intersect("SPP1", colnames(ligand_target_matrix))))
random_baseline_aupr <- length(genes_interest) / length(background)

activity_rows <- lapply(candidate_ligands, function(ligand) {
  scores <- ligand_target_matrix[background, ligand]
  positives <- background %in% genes_interest
  data.frame(
    ligand = ligand,
    aupr = average_precision(scores, positives),
    auroc = roc_auc(scores, positives),
    mean_positive_score = mean(scores[positives], na.rm = TRUE),
    mean_background_score = mean(scores, na.rm = TRUE),
    n_targets_positive = sum(positives),
    n_background = length(background),
    ligand_sender_obese_logcpm = unname(sender_obese[ligand] %||% NA_real_),
    stringsAsFactors = FALSE
  )
})
activity <- do.call(rbind, activity_rows)
activity$score_delta <- activity$mean_positive_score - activity$mean_background_score
activity <- activity[order(-activity$aupr, -activity$auroc, -activity$score_delta), ]
activity$rank_aupr <- seq_len(nrow(activity))
activity$is_focus_ligand <- activity$ligand == "SPP1"

receiver_reversal <- deg[deg$focus_group == "Mono_Macro" & deg$contrast == "enalapril_vs_obese", ]
receiver_reversal <- receiver_reversal[, c("gene_upper", "logFC", "cohen_d", "welch_p", "fdr")]
colnames(receiver_reversal) <- c("gene_upper", "enalapril_vs_obese_logFC", "enalapril_vs_obese_cohen_d", "enalapril_vs_obese_p", "enalapril_vs_obese_fdr")

receiver_disease <- deg_receiver[, c("gene_upper", "gene", "logFC", "cohen_d", "welch_p", "fdr", "mean_logcpm_a", "mean_logcpm_b")]
colnames(receiver_disease) <- c("gene_upper", "rat_gene_symbol", "obese_vs_lean_logFC", "obese_vs_lean_cohen_d", "obese_vs_lean_p", "obese_vs_lean_fdr", "obese_mean_logcpm", "lean_mean_logcpm")

spp1_scores <- ligand_target_matrix[background, "SPP1"]
targets <- data.frame(
  gene_upper = background,
  spp1_ligand_target_score = as.numeric(spp1_scores),
  in_receiver_obese_up_response = background %in% genes_interest,
  stringsAsFactors = FALSE
)
targets <- merge(targets, receiver_disease, by = "gene_upper", all.x = TRUE)
targets <- merge(targets, receiver_reversal, by = "gene_upper", all.x = TRUE)
targets$reversal_compatible <- targets$in_receiver_obese_up_response & targets$enalapril_vs_obese_logFC < 0
targets <- targets[order(-targets$in_receiver_obese_up_response, -targets$spp1_ligand_target_score, targets$obese_vs_lean_p), ]
top_targets <- head(targets[targets$in_receiver_obese_up_response, ], 100)

lr_spp1_cd44 <- lr_network[lr_network$from == "SPP1" & lr_network$to == "CD44", ]
spp1_rank <- activity[activity$ligand == "SPP1", ]
if (nrow(spp1_rank) == 0) {
  spp1_rank_line <- "SPP1 was absent from the candidate-ligand activity table."
} else {
  spp1_rank_line <- sprintf(
    "SPP1 ranked %d/%d candidate DTL-expressed ligands by AUPR (AUPR=%.4f; AUROC=%.4f; positive-target/background score delta=%.5f).",
    spp1_rank$rank_aupr[1], nrow(activity), spp1_rank$aupr[1], spp1_rank$auroc[1], spp1_rank$score_delta[1]
  )
}

reversal_count <- sum(top_targets$reversal_compatible, na.rm = TRUE)
top_n <- nrow(top_targets)

write.csv(activity, file.path(out_dir, "nichenet_ligand_activity_DTL_to_MonoMacro_obese_vs_lean.csv"), row.names = FALSE)
write.csv(top_targets, file.path(out_dir, "nichenet_SPP1_predicted_MonoMacro_targets_top100.csv"), row.names = FALSE)
write.csv(targets, file.path(out_dir, "nichenet_SPP1_all_receiver_target_scores.csv"), row.names = FALSE)
write.csv(lr_spp1_cd44, file.path(out_dir, "nichenet_lr_network_SPP1_CD44_evidence.csv"), row.names = FALSE)

plot_activity <- head(activity, 20)
plot_activity$ligand <- factor(plot_activity$ligand, levels = rev(plot_activity$ligand))
p1 <- ggplot(plot_activity, aes(x = ligand, y = aupr, fill = is_focus_ligand)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = random_baseline_aupr, linetype = "dashed", color = "#111827", linewidth = 0.35) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#b91c1c", `FALSE` = "#6b7280"), guide = "none") +
  labs(
    title = "NicheNet ligand activity: DTL sender to Mono/Macro receiver",
    subtitle = "Receiver response: Mono/Macro genes up in Obese vs Lean",
    x = NULL,
    y = "Average precision over receiver response genes"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "Figure_NicheNet_ligand_activity_DTL_to_MonoMacro.png"), p1, width = 7.2, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "Figure_NicheNet_ligand_activity_DTL_to_MonoMacro.pdf"), p1, width = 7.2, height = 5.2)

plot_targets <- head(top_targets, 30)
plot_targets$gene_upper <- factor(plot_targets$gene_upper, levels = rev(plot_targets$gene_upper))
p2 <- ggplot(plot_targets, aes(x = gene_upper, y = spp1_ligand_target_score)) +
  geom_col(aes(fill = reversal_compatible), width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#047857", `FALSE` = "#9ca3af"), na.value = "#9ca3af", name = "Enalapril-down") +
  labs(
    title = "Top SPP1-prioritized Mono/Macro receiver targets",
    subtitle = "Targets are selected from Mono/Macro Obese-vs-Lean up-response genes",
    x = NULL,
    y = "NicheNet ligand-target potential"
  ) +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "Figure_NicheNet_SPP1_predicted_targets.png"), p2, width = 7.2, height = 6.0, dpi = 300)
ggsave(file.path(out_dir, "Figure_NicheNet_SPP1_predicted_targets.pdf"), p2, width = 7.2, height = 6.0)

method_text <- c(
  "NicheNet ligand-to-target analysis (GSE233078)",
  "",
  "Purpose:",
  "This analysis asks whether ligands expressed by DTL sender cells, with emphasis on SPP1, can prioritize the Mono/Macro transcriptional response observed in obese metabolic CKD.",
  "",
  "Inputs:",
  "- Sender: DTL pseudobulk expression in Obese samples.",
  "- Receiver: Mono_Macro pseudobulk expression in Obese samples.",
  "- Receiver response genes: Mono_Macro Obese-vs-Lean up genes with nominal Welch p < 0.1, intersected with the NicheNet ligand-target matrix background.",
  "- Prior model: precomputed NicheNet ligand-target matrix and ligand-receptor network from the NicheNet Zenodo resource (https://zenodo.org/records/3260758).",
  paste0("- Random-baseline AUPR: ", sprintf("%.4f", random_baseline_aupr), " (receiver response genes / receiver background genes)."),
  "",
  "Boundary:",
  "This is a ligand-to-target prior-based computational prioritization analysis. It does not experimentally validate SPP1/CD44 signaling or causal perturbation."
)
write_utf8(file.path(out_dir, "nichenet_methods_note.txt"), method_text)

summary_lines <- c(
  "GSE233078 NicheNet ligand-to-target supplement summary",
  "",
  "Primary question:",
  "Can DTL-expressed ligands, especially SPP1, prioritize the Mono/Macro Obese-vs-Lean transcriptional response?",
  "",
  paste0("Candidate ligands tested: ", nrow(activity)),
  paste0("Receiver background genes: ", length(background)),
  paste0("Receiver response genes: ", length(genes_interest)),
  paste0("Random-baseline AUPR: ", sprintf("%.4f", random_baseline_aupr)),
  paste0("NicheNet LR evidence for SPP1->CD44 rows: ", nrow(lr_spp1_cd44)),
  spp1_rank_line,
  paste0("Top SPP1-positive receiver targets listed: ", top_n),
  paste0("Top SPP1-positive receiver targets with enalapril-down compatible direction: ", reversal_count, "/", top_n),
  "",
  "Top 10 ligand activity rows:",
  paste(capture.output(print(head(activity[, c("rank_aupr", "ligand", "aupr", "auroc", "score_delta", "ligand_sender_obese_logcpm")], 10), row.names = FALSE)), collapse = "\n"),
  "",
  "Top 20 SPP1-prioritized Mono/Macro targets:",
  paste(capture.output(print(head(top_targets[, c("gene_upper", "rat_gene_symbol", "spp1_ligand_target_score", "obese_vs_lean_logFC", "obese_vs_lean_p", "enalapril_vs_obese_logFC", "reversal_compatible")], 20), row.names = FALSE)), collapse = "\n"),
  "",
  "Writing boundary:",
  "Use wording such as 'ligand-to-target prediction prioritized SPP1-compatible Mono/Macro response genes'. Do not write 'NicheNet validated SPP1/CD44 causality' or 'virtual perturbation confirmed SPP1 function'."
)
write_utf8(file.path(out_dir, "nichenet_ligand_target_summary.txt"), summary_lines)
write_utf8(file.path(out_dir, "status.txt"), c("completed", Sys.time()))

