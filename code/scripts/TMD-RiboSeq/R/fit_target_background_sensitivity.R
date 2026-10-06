#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste(
    "Usage: fit_target_background_sensitivity.R <profiles.tsv.gz> <design.tsv>",
    "<results.tsv> <factor_a> <factor_b> <factor_a_reference>",
    "<factor_b_reference> <primary_center_nt>"
  ))
}

profiles_path <- args[[1]]
design_path <- args[[2]]
results_path <- args[[3]]
factor_a_name <- args[[4]]
factor_b_name <- args[[5]]
factor_a_reference <- args[[6]]
factor_b_reference <- args[[7]]
primary_center <- as.integer(args[[8]])

profiles <- read.delim(profiles_path, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
profiles <- profiles[
  profiles$center_nt == primary_center & profiles$sample_id %in% samples$sample_id,
]
if (nrow(profiles) == 0) stop("No primary-center model rows were found.")
if (anyDuplicated(profiles[, c("sample_id", "tmd_id")])) {
  stop("Duplicate sample_id × tmd_id rows at the primary center.")
}

samples[[factor_a_name]] <- relevel(factor(samples[[factor_a_name]]), factor_a_reference)
samples[[factor_b_name]] <- relevel(factor(samples[[factor_b_name]]), factor_b_reference)
if (nlevels(samples[[factor_a_name]]) != 2 || nlevels(samples[[factor_b_name]]) != 2) {
  stop("The sensitivity backend requires exactly two levels for each factor.")
}
sample_order <- samples$sample_id
tmd_ids <- sort(unique(profiles$tmd_id))
complete_tmds <- vapply(tmd_ids, function(id) {
  rows <- profiles[profiles$tmd_id == id,]
  all(sample_order %in% rows$sample_id) && nrow(rows) == length(sample_order) &&
    all(as.logical(rows$coverage_pass))
}, logical(1))
tmd_ids <- tmd_ids[complete_tmds]
if (length(tmd_ids) < 2) stop("Fewer than two strict-coverage TMDs are available.")
profiles <- profiles[profiles$tmd_id %in% tmd_ids,]

# Match the primary and positional models by removing TMDs with no
# target-window footprints in any model sample. Background-only information
# cannot identify a target-versus-background differential-usage effect.
target_totals <- tapply(profiles$target_count, profiles$tmd_id, sum)
nonzero_target <- as.numeric(target_totals[tmd_ids]) > 0
tmd_ids <- tmd_ids[nonzero_target]
if (length(tmd_ids) < 2) {
  stop("Fewer than two TMDs have nonzero target-window counts across model samples.")
}
profiles <- profiles[profiles$tmd_id %in% tmd_ids,]

target <- data.frame(
  tmd_id = profiles$tmd_id,
  sample_id = profiles$sample_id,
  region = "target",
  count = profiles$target_count
)
background <- data.frame(
  tmd_id = profiles$tmd_id,
  sample_id = profiles$sample_id,
  region = "background",
  count = profiles$background_count
)
features <- rbind(target, background)
features$feature_id <- paste(features$tmd_id, features$region, sep = "__")
features$sample_id <- factor(features$sample_id, levels = sample_order)
counts <- xtabs(count ~ feature_id + sample_id, data = features)
counts <- counts[, sample_order, drop = FALSE]
storage.mode(counts) <- "integer"

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

y <- DGEList(counts = counts)
y <- normLibSizes(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
gene_ids <- sub("__(target|background)$", "", rownames(counts))
region_ids <- sub("^.*__", "", rownames(counts))

run_usage_test <- function(contrast, prefix) {
  test <- diffSpliceDGE(
    fit,
    contrast = contrast,
    geneid = gene_ids,
    exonid = region_ids,
    robust = TRUE,
    verbose = FALSE
  )
  genes <- topSpliceDGE(test, test = "Simes", number = Inf)
  names(genes)[names(genes) == "GeneID"] <- "tmd_id"
  names(genes)[names(genes) == "P.Value"] <- paste0(prefix, "_p")
  names(genes)[names(genes) == "FDR"] <- paste0(prefix, "_fdr")
  genes <- genes[, c("tmd_id", paste0(prefix, "_p"), paste0(prefix, "_fdr"))]

  coefficient <- as.numeric(test$coefficients)
  names(coefficient) <- paste(gene_ids, region_ids, sep = "__")
  effect <- vapply(unique(gene_ids), function(id) {
    coefficient[paste0(id, "__target")] - coefficient[paste0(id, "__background")]
  }, numeric(1))
  effect_table <- data.frame(
    tmd_id = names(effect),
    value = as.numeric(effect),
    stringsAsFactors = FALSE
  )
  names(effect_table)[2] <- paste0(prefix, "_target_minus_background_log2_effect")
  merge(genes, effect_table, by = "tmd_id", all.x = TRUE, sort = FALSE)
}

reference_result <- run_usage_test(reference_b_contrast, "reference_b")
alternate_result <- run_usage_test(alternate_b_contrast, "alternate_b")
interaction_result <- run_usage_test(interaction_contrast, "interaction")
result <- Reduce(
  function(x, y) merge(x, y, by = "tmd_id", all = TRUE, sort = FALSE),
  list(reference_result, alternate_result, interaction_result)
)
annotation <- unique(profiles[, c("tmd_id", "gene_id", "transcript_id")])
result <- merge(annotation, result, by = "tmd_id", all.y = TRUE, sort = FALSE)
result <- result[match(tmd_ids, result$tmd_id),]

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
write.table(result, results_path, sep = "\t", row.names = FALSE, quote = FALSE)
