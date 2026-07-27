# =============================================================================
# R/batch.R -- batch/confounding handling
# -----------------------------------------------------------------------------
# One concern: (1) identify the nuisance (batch) variables in the design, and
# (2) produce a batch-corrected matrix FOR VISUALIZATION ONLY.
#
# CRITICAL distinction: batch is handled for DIFFERENTIAL EXPRESSION by INCLUDING
# it in the DESeq2 design (e.g. ~ batch + condition) -- DESeq2 regresses it out
# inside the GLM. We do NOT feed batch-corrected values to DESeq2. limma's
# removeBatchEffect is used purely to make PCA/heatmaps show the biological
# structure after nuisance variation is visually subtracted. Using corrected
# values as DE input would invalidate the variance estimates.
# =============================================================================

suppressPackageStartupMessages(library(limma))

#' Identify batch (nuisance) variables: every design term except the tested one.
#'
#' @param design_str Model formula string.
#' @param tested_factor The variable of interest (last term).
#' @return Character vector of nuisance variable names (possibly empty).
identify_batch_vars <- function(design_str, tested_factor) {
  vars <- all.vars(stats::as.formula(design_str))
  setdiff(vars, tested_factor)
}

#' Remove batch effects from a VST matrix for visualization.
#'
#' @param vsd A DESeqTransform (from get_vst()).
#' @param metadata Sample sheet, row-aligned to the VST columns.
#' @param batch_vars Nuisance variable name(s) to remove.
#' @param tested_factor Biological variable to PRESERVE.
#' @return A DESeqTransform copy whose assay() has the batch effect removed.
#' We pass the biological design to removeBatchEffect so it does not accidentally
#' regress out the treatment effect while removing batch. Supports one or two
#' batch variables (limma's documented capability).
remove_batch_for_viz <- function(vsd, metadata, batch_vars, tested_factor) {
  mat <- SummarizedExperiment::assay(vsd)
  md <- metadata[match(colnames(mat), rownames(metadata)), , drop = FALSE]

  # Preserve the biological contrast: model.matrix of the tested factor.
  design <- if (tested_factor %in% colnames(md)) {
    stats::model.matrix(stats::reformulate(tested_factor), data = md)
  } else {
    matrix(1, nrow = ncol(mat))
  }

  b1 <- factor(md[[batch_vars[1]]])
  b2 <- if (length(batch_vars) >= 2) factor(md[[batch_vars[2]]]) else NULL

  corrected <- limma::removeBatchEffect(mat, batch = b1, batch2 = b2, design = design)

  vsd2 <- vsd
  SummarizedExperiment::assay(vsd2) <- corrected
  log_success(sprintf("Batch effect removed for visualization (batch var(s): %s).",
                     paste(batch_vars, collapse = ", ")))
  vsd2
}

#' Convenience: run the identifiability QC check and, if a batch variable exists,
#' return a batch-corrected VST for a side-by-side "before/after" PCA.
#'
#' @return list(batch_vars, vsd_corrected or NULL).
handle_batch <- function(qc, vsd, metadata, design_str, tested_factor) {
  # (1) Identifiability is a QC FAIL check (confounded batch cannot be separated).
  qc_check_design_identifiable(qc, as.data.frame(metadata), design_str, tested_factor)

  # (2) If there is a nuisance variable, build the corrected view.
  batch_vars <- identify_batch_vars(design_str, tested_factor)
  batch_vars <- batch_vars[batch_vars %in% colnames(metadata)]
  if (length(batch_vars) == 0) {
    log_info("No batch variable in the design; skipping batch-correction view.")
    return(list(batch_vars = character(0), vsd_corrected = NULL))
  }
  vsd_c <- remove_batch_for_viz(vsd, as.data.frame(metadata), batch_vars, tested_factor)
  list(batch_vars = batch_vars, vsd_corrected = vsd_c)
}
