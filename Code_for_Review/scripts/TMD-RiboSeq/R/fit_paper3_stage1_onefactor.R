#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: fit_paper3_stage1_onefactor.R <input.tsv> <results.tsv> <summary.tsv>")
}

input_path <- args[[1]]
results_path <- args[[2]]
summary_path <- args[[3]]

d <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "analysis_family", "feature_id", "transcript", "gene_id", "sample",
  "condition", "target_count", "background_count", "target_length",
  "background_length", "report_inference"
)
missing <- setdiff(required, names(d))
if (length(missing) > 0) {
  stop(paste("Stage 1 model input is missing:", paste(missing, collapse = ", ")))
}
if (nrow(d) == 0) stop("Stage 1 model input is empty.")

fit_family <- function(family) {
  x <- d[d$analysis_family == family, , drop = FALSE]
  samples <- unique(x[, c("sample", "condition")])
  samples <- samples[order(match(samples$sample, unique(d$sample))), , drop = FALSE]
  if (anyDuplicated(samples$sample)) stop(paste("Duplicate design sample in", family))
  if (length(unique(samples$condition)) != 2) {
    stop(paste("Family does not contain exactly two conditions:", family))
  }
  samples$condition <- relevel(factor(samples$condition), "Control")
  if (!identical(levels(samples$condition), c("Control", "Knockdown"))) {
    stop(paste("Unexpected condition levels in", family))
  }
  design <- model.matrix(~ condition, data = samples)
  if (!identical(colnames(design), c("(Intercept)", "conditionKnockdown"))) {
    stop(paste("Unexpected design columns in", family))
  }
  if (qr(design)$rank < ncol(design)) stop(paste("Rank-deficient design in", family))

  feature_ids <- sort(unique(x$feature_id))
  sample_order <- samples$sample
  if (anyDuplicated(x[, c("feature_id", "sample")])) {
    stop(paste("Duplicate feature x sample row in", family))
  }
  complete <- vapply(feature_ids, function(id) {
    y <- x[x$feature_id == id, , drop = FALSE]
    nrow(y) == length(sample_order) && all(sample_order %in% y$sample)
  }, logical(1))
  if (!all(complete)) stop(paste("Incomplete feature rows in", family))
  x$feature_id <- factor(x$feature_id, levels = feature_ids)
  x$sample <- factor(x$sample, levels = sample_order)
  make_matrix <- function(column) {
    z <- xtabs(as.formula(paste(column, "~ feature_id + sample")), data = x)
    z[feature_ids, sample_order, drop = FALSE]
  }
  counts <- make_matrix("target_count")
  storage.mode(counts) <- "integer"
  background <- make_matrix("background_count")
  target_length <- make_matrix("target_length")
  background_length <- make_matrix("background_length")
  if (any(background < 20)) stop(paste("Ineligible background count in", family))
  if (any(background_length <= 0)) stop(paste("Non-positive background length in", family))
  if (any(target_length != 31)) stop(paste("Non-frozen target length in", family))
  if (any(rowSums(counts) == 0)) stop(paste("All-zero target entered inference in", family))
  exposure <- (background + 0.5) * target_length / background_length
  if (any(!is.finite(exposure)) || any(exposure <= 0)) {
    stop(paste("Invalid exposure in", family))
  }

  y <- DGEList(counts = counts)
  y$offset <- log(exposure)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- topTags(
    glmQLFTest(fit, coef = "conditionKnockdown"),
    n = Inf,
    sort.by = "none"
  )$table
  report_inference <- unique(x$report_inference)
  if (length(report_inference) != 1 || !report_inference %in% c("YES", "NO")) {
    stop(paste("Invalid report_inference policy in", family))
  }
  reported_p <- if (report_inference == "YES") test$PValue else rep(NA_real_, nrow(test))
  reported_fdr <- if (report_inference == "YES") {
    p.adjust(test$PValue, method = "BH")
  } else {
    rep(NA_real_, nrow(test))
  }
  annotation <- unique(x[, c("feature_id", "transcript", "gene_id")])
  annotation <- annotation[match(feature_ids, annotation$feature_id), , drop = FALSE]
  data.frame(
    analysis_family = family,
    feature_id = feature_ids,
    transcript = annotation$transcript,
    gene_id = annotation$gene_id,
    coefficient = "conditionKnockdown",
    contrast = "Knockdown-Control",
    logFC = test$logFC,
    logCPM = test$logCPM,
    F = test$F,
    PValue = reported_p,
    FDR = reported_fdr,
    common_dispersion = if (length(y$common.dispersion)) y$common.dispersion else NA_real_,
    median_tagwise_dispersion = median(y$tagwise.dispersion),
    median_ql_dispersion = median(fit$s2.post)
  )
}

families <- unique(d$analysis_family)
results <- do.call(rbind, lapply(families, fit_family))
write.table(results, results_path, sep = "\t", row.names = FALSE, quote = FALSE)
summary <- do.call(rbind, lapply(split(results, results$analysis_family), function(x) {
  data.frame(
    analysis_family = x$analysis_family[[1]],
    tested_features = nrow(x),
    bh_fdr_lt_0_05 = if (all(is.na(x$FDR))) NA_integer_ else sum(x$FDR < 0.05),
    median_logFC = median(x$logFC),
    coefficient = "conditionKnockdown",
    contrast = "Knockdown-Control"
  )
}))
write.table(summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
