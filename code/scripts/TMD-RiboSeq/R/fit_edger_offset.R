#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10) {
  stop(paste(
    "Usage: fit_edger_offset.R <profiles.tsv.gz> <design.tsv> <results.tsv>",
    "<summary.tsv> <factor_a> <factor_b> <factor_a_reference>",
    "<factor_b_reference> <primary_center_nt> <fdr_threshold>"
  ))
}

profiles_path <- args[[1]]
design_path <- args[[2]]
results_path <- args[[3]]
summary_path <- args[[4]]
factor_a_name <- args[[5]]
factor_b_name <- args[[6]]
factor_a_reference <- args[[7]]
factor_b_reference <- args[[8]]
primary_center <- as.integer(args[[9]])
fdr_threshold <- as.numeric(args[[10]])

profiles <- read.delim(profiles_path, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)

required_profile <- c(
  "sample_id", "tmd_id", "gene_id", "center_nt", "target_count",
  "target_length_nt", "background_count", "background_length_nt",
  "coverage_pass"
)
missing_profile <- setdiff(required_profile, names(profiles))
if (length(missing_profile) > 0) {
  stop(paste("Profile table is missing:", paste(missing_profile, collapse = ", ")))
}
required_design <- c("sample_id", factor_a_name, factor_b_name)
missing_design <- setdiff(required_design, names(samples))
if (length(missing_design) > 0) {
  stop(paste("Design table is missing:", paste(missing_design, collapse = ", ")))
}

profiles <- profiles[
  profiles$center_nt == primary_center & profiles$sample_id %in% samples$sample_id,
]
if (nrow(profiles) == 0) stop("No primary-center model rows were found.")
if (anyDuplicated(profiles[, c("sample_id", "tmd_id")])) {
  stop("Profile table contains duplicate sample_id × tmd_id rows at the primary center.")
}

samples <- samples[match(unique(samples$sample_id), samples$sample_id),]
samples[[factor_a_name]] <- relevel(factor(samples[[factor_a_name]]), factor_a_reference)
samples[[factor_b_name]] <- relevel(factor(samples[[factor_b_name]]), factor_b_reference)
if (nlevels(samples[[factor_a_name]]) != 2 || nlevels(samples[[factor_b_name]]) != 2) {
  stop("The validated v0.1 backend requires exactly two levels for each factor.")
}

sample_order <- samples$sample_id
tmd_ids <- sort(unique(profiles$tmd_id))
complete_tmds <- vapply(tmd_ids, function(id) {
  rows <- profiles[profiles$tmd_id == id,]
  all(sample_order %in% rows$sample_id) &&
    nrow(rows) == length(sample_order) &&
    all(as.logical(rows$coverage_pass))
}, logical(1))
tmd_ids <- tmd_ids[complete_tmds]
if (length(tmd_ids) < 2) {
  stop("Fewer than two TMDs pass the strict all-model-sample coverage rule.")
}
coverage_eligible_tmds <- length(tmd_ids)

profiles <- profiles[
  profiles$tmd_id %in% tmd_ids,
]
profiles$sample_id <- factor(profiles$sample_id, levels = sample_order)
profiles$tmd_id <- factor(profiles$tmd_id, levels = tmd_ids)

make_matrix <- function(column) {
  formula <- as.formula(paste(column, "~ tmd_id + sample_id"))
  matrix <- xtabs(formula, data = profiles)
  matrix[tmd_ids, sample_order, drop = FALSE]
}

counts <- make_matrix("target_count")
storage.mode(counts) <- "integer"
background <- make_matrix("background_count")
target_length <- make_matrix("target_length_nt")
background_length <- make_matrix("background_length_nt")

# A TMD can satisfy the background-coverage rule while having zero target-window
# footprints in every model sample. Such a row has no information for estimating
# a target-window differential effect and is excluded consistently across the
# primary, sensitivity, and positional models.
nonzero_target <- rowSums(counts) > 0
if (sum(nonzero_target) < 2) {
  stop("Fewer than two TMDs have nonzero target-window counts across model samples.")
}
tmd_ids <- tmd_ids[nonzero_target]
counts <- counts[nonzero_target, , drop = FALSE]
background <- background[nonzero_target, , drop = FALSE]
target_length <- target_length[nonzero_target, , drop = FALSE]
background_length <- background_length[nonzero_target, , drop = FALSE]

if (any(background_length <= 0)) stop("Non-positive background length in eligible data.")
exposure <- (background + 0.5) * target_length / background_length
if (any(!is.finite(exposure)) || any(exposure <= 0)) stop("Invalid offset exposure.")

formula <- as.formula(paste("~", factor_a_name, "*", factor_b_name))
design <- model.matrix(formula, data = samples)
if (qr(design)$rank < ncol(design)) stop("The design matrix is not full rank.")

factor_a_levels <- levels(samples[[factor_a_name]])
factor_b_levels <- levels(samples[[factor_b_name]])
factor_a_alternate <- factor_a_levels[[2]]
factor_b_alternate <- factor_b_levels[[2]]
new_data <- data.frame(
  a = factor(rep(factor_a_levels, each = 2), levels = factor_a_levels),
  b = factor(rep(factor_b_levels, times = 2), levels = factor_b_levels)
)
names(new_data) <- c(factor_a_name, factor_b_name)
cell_design <- model.matrix(formula, data = new_data)
rownames(cell_design) <- c("ref_base", "ref_alt_b", "alt_a_base", "alt_a_alt_b")
reference_b_contrast <- cell_design["ref_alt_b",] - cell_design["ref_base",]
alternate_b_contrast <- cell_design["alt_a_alt_b",] - cell_design["alt_a_base",]
interaction_contrast <- alternate_b_contrast - reference_b_contrast

y <- DGEList(counts = counts)
y$offset <- log(exposure)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

reference_test <- topTags(
  glmQLFTest(fit, contrast = reference_b_contrast), n = Inf, sort.by = "none"
)$table
alternate_test <- topTags(
  glmQLFTest(fit, contrast = alternate_b_contrast), n = Inf, sort.by = "none"
)$table
interaction_test <- topTags(
  glmQLFTest(fit, contrast = interaction_contrast), n = Inf, sort.by = "none"
)$table

annotation <- unique(profiles[, c("tmd_id", "gene_id", "transcript_id")])
annotation <- annotation[match(tmd_ids, annotation$tmd_id),]
results <- data.frame(
  tmd_id = tmd_ids,
  gene_id = annotation$gene_id,
  transcript_id = annotation$transcript_id,
  factor_a = factor_a_name,
  factor_a_reference = factor_a_reference,
  factor_a_alternate = factor_a_alternate,
  factor_b = factor_b_name,
  factor_b_reference = factor_b_reference,
  factor_b_alternate = factor_b_alternate,
  reference_b_log2_effect = reference_test$logFC,
  reference_b_p = reference_test$PValue,
  reference_b_fdr = reference_test$FDR,
  alternate_b_log2_effect = alternate_test$logFC,
  alternate_b_p = alternate_test$PValue,
  alternate_b_fdr = alternate_test$FDR,
  interaction_log2_effect = interaction_test$logFC,
  interaction_p = interaction_test$PValue,
  interaction_fdr = interaction_test$FDR
)

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, results_path, sep = "\t", row.names = FALSE, quote = FALSE)
scalar_or_na <- function(value) if (length(value) == 0) NA_real_ else as.numeric(value[[1]])
summary <- data.frame(
  metric = c(
    "fdr_threshold", "coverage_eligible_tmds", "tested_tmds",
    "reference_b_significant", "alternate_b_significant",
    "interaction_significant", "common_dispersion", "median_tagwise_dispersion",
    "median_ql_dispersion"
  ),
  value = c(
    fdr_threshold, coverage_eligible_tmds, nrow(results),
    sum(results$reference_b_fdr < fdr_threshold),
    sum(results$alternate_b_fdr < fdr_threshold),
    sum(results$interaction_fdr < fdr_threshold),
    scalar_or_na(y$common.dispersion), scalar_or_na(median(y$tagwise.dispersion)),
    scalar_or_na(median(fit$s2.post))
  )
)
write.table(summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
