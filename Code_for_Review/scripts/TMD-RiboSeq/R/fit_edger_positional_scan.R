#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste(
    "Usage: fit_edger_positional_scan.R <profiles.tsv.gz> <design.tsv>",
    "<results.tsv.gz> <summary.tsv> <factor_a> <factor_b>",
    "<factor_a_reference> <factor_b_reference> <fdr_threshold>"
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
fdr_threshold <- as.numeric(args[[9]])

profiles <- read.delim(profiles_path, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
samples[[factor_a_name]] <- relevel(factor(samples[[factor_a_name]]), factor_a_reference)
samples[[factor_b_name]] <- relevel(factor(samples[[factor_b_name]]), factor_b_reference)
if (nlevels(samples[[factor_a_name]]) != 2 || nlevels(samples[[factor_b_name]]) != 2) {
  stop("The positional backend requires exactly two levels for each factor.")
}
sample_order <- samples$sample_id

formula <- as.formula(paste("~", factor_a_name, "*", factor_b_name))
design <- model.matrix(formula, data = samples)
if (qr(design)$rank < ncol(design)) stop("The design matrix is not full rank.")
factor_a_levels <- levels(samples[[factor_a_name]])
factor_b_levels <- levels(samples[[factor_b_name]])
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

fit_center <- function(center) {
  d <- profiles[profiles$center_nt == center & profiles$sample_id %in% sample_order,]
  if (anyDuplicated(d[, c("sample_id", "tmd_id")])) {
    stop(paste("Duplicate sample × TMD rows at center", center))
  }
  tmd_ids <- sort(unique(d$tmd_id))
  keep <- vapply(tmd_ids, function(id) {
    rows <- d[d$tmd_id == id,]
    all(sample_order %in% rows$sample_id) && nrow(rows) == length(sample_order) &&
      all(as.logical(rows$coverage_pass))
  }, logical(1))
  tmd_ids <- tmd_ids[keep]
  if (length(tmd_ids) < 2) return(NULL)
  d <- d[d$tmd_id %in% tmd_ids,]
  d$sample_id <- factor(d$sample_id, levels = sample_order)
  d$tmd_id <- factor(d$tmd_id, levels = tmd_ids)
  make_matrix <- function(column) {
    matrix <- xtabs(as.formula(paste(column, "~ tmd_id + sample_id")), data = d)
    matrix[tmd_ids, sample_order, drop = FALSE]
  }
  counts <- make_matrix("target_count")
  storage.mode(counts) <- "integer"
  background <- make_matrix("background_count")
  target_length <- make_matrix("target_length_nt")
  background_length <- make_matrix("background_length_nt")
  nonzero <- rowSums(counts) > 0
  if (sum(nonzero) < 2) return(NULL)
  tmd_ids <- tmd_ids[nonzero]
  counts <- counts[nonzero,, drop = FALSE]
  background <- background[nonzero,, drop = FALSE]
  target_length <- target_length[nonzero,, drop = FALSE]
  background_length <- background_length[nonzero,, drop = FALSE]
  exposure <- (background + 0.5) * target_length / background_length
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
  annotation <- unique(d[, c("tmd_id", "gene_id", "transcript_id")])
  annotation <- annotation[match(tmd_ids, annotation$tmd_id),]
  data.frame(
    center_nt = center,
    tmd_id = tmd_ids,
    gene_id = annotation$gene_id,
    transcript_id = annotation$transcript_id,
    reference_b_log2_effect = reference_test$logFC,
    reference_b_p = reference_test$PValue,
    alternate_b_log2_effect = alternate_test$logFC,
    alternate_b_p = alternate_test$PValue,
    interaction_log2_effect = interaction_test$logFC,
    interaction_p = interaction_test$PValue
  )
}

results <- lapply(sort(unique(profiles$center_nt)), fit_center)
results <- results[!vapply(results, is.null, logical(1))]
if (length(results) == 0) stop("No positional center had sufficient strict coverage.")
results <- do.call(rbind, results)
results$reference_b_center_fdr <- ave(
  results$reference_b_p, results$center_nt, FUN = function(x) p.adjust(x, "BH")
)
results$alternate_b_center_fdr <- ave(
  results$alternate_b_p, results$center_nt, FUN = function(x) p.adjust(x, "BH")
)
results$interaction_center_fdr <- ave(
  results$interaction_p, results$center_nt, FUN = function(x) p.adjust(x, "BH")
)
results$reference_b_global_fdr <- p.adjust(results$reference_b_p, "BH")
results$alternate_b_global_fdr <- p.adjust(results$alternate_b_p, "BH")
results$interaction_global_fdr <- p.adjust(results$interaction_p, "BH")

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(results_path, "wt")
write.table(results, connection, sep = "\t", row.names = FALSE, quote = FALSE)
close(connection)
summary <- do.call(rbind, lapply(split(results, results$center_nt), function(d) {
  data.frame(
    center_nt = d$center_nt[[1]],
    tested_tmds = nrow(d),
    fdr_threshold = fdr_threshold,
    reference_center_significant = sum(d$reference_b_center_fdr < fdr_threshold),
    alternate_center_significant = sum(d$alternate_b_center_fdr < fdr_threshold),
    interaction_center_significant = sum(d$interaction_center_fdr < fdr_threshold),
    interaction_global_significant = sum(d$interaction_global_fdr < fdr_threshold)
  )
}))
write.table(summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
