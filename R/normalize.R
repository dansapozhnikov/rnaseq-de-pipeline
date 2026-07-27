# =============================================================================
# R/normalize.R -- build the DESeq2 object, filter, normalize, transform
# -----------------------------------------------------------------------------
# One concern: everything from a raw integer count matrix + metadata to a fitted
# DESeqDataSet plus a variance-stabilized matrix for visualization.
#
# Key correctness point: DESeq2 models RAW counts and estimates its own size
# factors internally; we therefore feed it raw counts and NEVER pre-normalize.
# VST is used ONLY for distance-based visualization (PCA/heatmap), never for DE.
# =============================================================================

suppressPackageStartupMessages(library(DESeq2))

#' Build a DESeqDataSet from a raw count matrix and sample sheet.
#'
#' @param counts Integer matrix (genes x samples), columns aligned to metadata.
#' @param metadata data.frame of sample annotations (row-aligned to counts cols).
#' @param design_str Model formula string, tested variable last.
#' @param contrast [factor, numerator, denominator]; denominator becomes the
#'   reference level so log2FoldChange has the documented direction.
#' @return A DESeqDataSet (not yet fitted).
build_dds <- function(counts, metadata, design_str, contrast) {
  design <- stats::as.formula(design_str)
  md <- metadata
  # Coerce every variable used in the design to a factor; DESeq2 needs factors to
  # build the model matrix, and character columns would be silently dropped.
  for (v in all.vars(design)) {
    if (v %in% colnames(md)) md[[v]] <- factor(md[[v]])
  }
  tested_factor <- contrast[[1]]
  denom <- contrast[[3]]
  # Set the reference level to the contrast's denominator. WHY: DESeq2 reports
  # log2(numerator/denominator); making the denominator the reference keeps the
  # sign convention explicit and matches the results() contrast we request later.
  if (tested_factor %in% colnames(md)) {
    if (!denom %in% levels(md[[tested_factor]])) {
      stop_pipeline(sprintf("Contrast denominator '%s' is not a level of factor '%s'.",
                            denom, tested_factor))
    }
    md[[tested_factor]] <- stats::relevel(md[[tested_factor]], ref = denom)
  }
  rownames(md) <- colnames(counts)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = md, design = design)
  dds
}

#' Filter out low-count genes.
#'
#' @param dds A DESeqDataSet.
#' @param cfg Config (uses qc$min_count and qc$min_replicates_per_group).
#' @return list(dds = filtered object, n_before, n_after).
#' Keep a gene only if it has at least `min_count` reads in at least
#' `min_replicates_per_group` samples. WHY: genes that are essentially zero
#' everywhere carry no information, inflate multiple-testing correction, and
#' destabilize dispersion estimation. Using the smallest group size as the
#' sample threshold means a gene expressed in even one full group is retained.
filter_low_counts <- function(dds, cfg) {
  n_before <- nrow(dds)
  min_count <- cfg$qc$min_count
  min_samples <- cfg$qc$min_replicates_per_group
  keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
  dds <- dds[keep, ]
  list(dds = dds, n_before = n_before, n_after = nrow(dds))
}

#' Fit the DESeq2 model (size factors, dispersions, negative-binomial GLM).
#'
#' @param dds A DESeqDataSet.
#' @param test "Wald" (default) or "LRT".
#' @param reduced Reduced-model formula string; REQUIRED when test="LRT".
#' @return A fitted DESeqDataSet.
#' Wald tests a single coefficient/contrast (good for two-group comparisons and
#' specific effects). LRT (likelihood-ratio test) compares the full design to a
#' `reduced` one and tests the SET of terms dropped in a single p-value -- the
#' right tool for multi-level factors, time courses, and ANOVA-style questions.
run_deseq <- function(dds, test = "Wald", reduced = NULL) {
  test <- match.arg(test, c("Wald", "LRT"))
  if (test == "LRT") {
    if (is.null(reduced)) {
      stop_pipeline("de$test is 'LRT' but de$reduced (the reduced-model formula) is not set.")
    }
    log_info(sprintf("Fitting DESeq2 (LRT: full design vs reduced '%s')...", reduced))
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = stats::as.formula(reduced))
  } else {
    log_info("Fitting DESeq2 model (Wald: size factors -> dispersions -> NB GLM)...")
    dds <- DESeq2::DESeq(dds)
  }
  log_success(sprintf("Model fitted (%s). Size factors: %s", test,
                     paste(sprintf("%s=%.2f", names(DESeq2::sizeFactors(dds)),
                                   DESeq2::sizeFactors(dds)), collapse = ", ")))
  dds
}

#' Variance-stabilizing transform for visualization/clustering.
#'
#' @param dds A (fitted or unfitted) DESeqDataSet.
#' @param blind TRUE = ignore the design (unsupervised QC view, the default).
#' @return A DESeqTransform whose assay() is the VST matrix.
#' WHY VST and not log-CPM: VST removes the strong mean-variance dependence of
#' counts so that PCA/heatmap distances are not dominated by a handful of very
#' highly expressed genes. blind=TRUE keeps the QC view honest (design-agnostic).
get_vst <- function(dds, blind = TRUE) {
  # vst() needs enough genes to fit the dispersion trend; fall back to the exact
  # varianceStabilizingTransformation for very small inputs (e.g. tiny fixtures).
  vsd <- tryCatch(
    DESeq2::vst(dds, blind = blind),
    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = blind))
  vsd
}
