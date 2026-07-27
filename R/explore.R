# =============================================================================
# R/explore.R — exploratory analysis: PCA + variance-vs-metadata association
# -----------------------------------------------------------------------------
# One concern: unsupervised views that answer "does the data cluster the way the
# biology says it should, and is any technical variable (batch) driving the
# major axes of variation?" Both plots are saved as PNGs for the HTML report.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(matrixStats)
})

#' Compute a PCA on the most variable genes of a VST matrix.
#'
#' @param vsd A DESeqTransform (from get_vst()).
#' @param ntop Number of top-variance genes to use (default 500).
#' @return list(pca = prcomp object, percentVar = variance % per PC, scores = df).
#' WHY top-variance genes: PCA on all genes is dominated by noise from the long
#' tail of low-variance genes; restricting to the most variable genes surfaces
#' the structure (treatment, batch) that actually separates samples.
compute_pca <- function(vsd, ntop = 500) {
  mat <- SummarizedExperiment::assay(vsd)
  vars <- matrixStats::rowVars(mat)
  ntop <- min(ntop, nrow(mat))
  top <- order(vars, decreasing = TRUE)[seq_len(ntop)]
  pca <- stats::prcomp(t(mat[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  scores <- as.data.frame(pca$x)
  scores$sample <- rownames(scores)
  list(pca = pca, percentVar = percent, scores = scores)
}

#' Draw and save a PC1-vs-PC2 scatter colored by the tested factor.
#'
#' @param pca_res Output of compute_pca().
#' @param metadata Sample sheet (row-aligned to samples).
#' @param color_by Factor mapped to color (usually the tested variable).
#' @param shape_by Optional factor mapped to shape (usually batch).
#' @param outdir Output root; the plot is written under <outdir>/plots/.
#' @param filename Output PNG name under <outdir>/plots/ (default "pca.png").
#' @param save_scores If TRUE, also persist PCA scores for the interactive report
#'   (only wanted for the primary PCA, not the batch-corrected variant).
#' @return Path to the saved PNG.
plot_pca <- function(pca_res, metadata, color_by, shape_by = NULL, outdir,
                     filename = "pca.png", save_scores = TRUE) {
  df <- pca_res$scores
  # Bind metadata columns (avoid clobbering the existing 'sample' column).
  md <- metadata[match(df$sample, rownames(metadata)), , drop = FALSE]
  md <- md[, setdiff(colnames(md), colnames(df)), drop = FALSE]
  df <- cbind(df, md)
  pv <- pca_res$percentVar

  # Build the aesthetic with the .data pronoun so column names come from strings
  # robustly across ggplot2 versions (2.x tidy-eval and 4.x alike).
  mapping <- aes(x = .data[["PC1"]], y = .data[["PC2"]], color = .data[[color_by]])
  if (!is.null(shape_by) && shape_by %in% colnames(df)) {
    mapping <- aes(x = .data[["PC1"]], y = .data[["PC2"]],
                   color = .data[[color_by]], shape = .data[[shape_by]])
  }

  p <- ggplot(df, mapping) +
    geom_point(size = 4, alpha = 0.9) +
    labs(x = sprintf("PC1: %s%% variance", pv[1]),
         y = sprintf("PC2: %s%% variance", pv[2]),
         title = "PCA of variance-stabilized counts (top 500 variable genes)") +
    theme_bw(base_size = 13)

  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(plots_dir, filename)
  ggplot2::ggsave(path, p, width = 7, height = 5.5, dpi = 150)

  if (save_scores) {
    # Persist the scores + metadata + variance percentages so the HTML report can
    # rebuild this as an INTERACTIVE plotly scatter (hover = sample id + group).
    tables_dir <- file.path(outdir, "tables")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    scores_out <- df[, unique(c("sample", grep("^PC[1-9]", colnames(df), value = TRUE),
                                colnames(metadata)))]
    utils::write.table(scores_out, file.path(tables_dir, "pca_scores.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(paste(pv, collapse = "\t"), file.path(tables_dir, "pca_percentvar.txt"))
  }
  log_success(sprintf("PCA plot saved to %s", path))
  path
}

#' Associate each principal component with each metadata variable.
#'
#' @param pca_res Output of compute_pca().
#' @param metadata Sample sheet.
#' @param npc Number of leading PCs to test (default 5).
#' @param outdir Output root.
#' @return list(assoc = matrix of R^2 [PC x variable], plot = png path, table = tsv path).
#' For every (PC, variable) pair we fit PC ~ variable and record the R^2, i.e.
#' the fraction of that PC's spread explained by the variable. WHY: this is the
#' quantitative batch-effect diagnostic — if a technical variable explains most
#' of PC1, the major axis of variation is technical, not biological.
variance_vs_metadata <- function(pca_res, metadata, npc = 5, outdir) {
  scores <- pca_res$scores
  npc <- min(npc, sum(grepl("^PC", colnames(scores))))
  pcs <- paste0("PC", seq_len(npc))
  md <- metadata[match(scores$sample, rownames(metadata)), , drop = FALSE]
  n <- nrow(md)

  # Keep only INFORMATIVE variables. WHY: a metadata column that is a per-sample
  # unique identifier (SampleName, Run, BioSample, ...) has one level per sample,
  # so regressing any PC on it fits perfectly (R^2 = 1) and tells us nothing about
  # structure. We therefore drop constants (1 level) AND identifiers (n levels),
  # keeping genuine experimental/technical factors and continuous covariates.
  informative <- function(x) {
    u <- length(unique(x[!is.na(x)]))
    is.numeric(x) && u > 1 || (u > 1 && u < n)
  }
  vars <- colnames(md)[vapply(md, informative, logical(1))]
  vars <- setdiff(vars, "sample")
  if (length(vars) == 0) {
    log_warn("No informative metadata variables for variance-vs-metadata; skipping.")
    return(list(assoc = NULL, plot = NULL, table = NULL))
  }

  assoc <- matrix(NA_real_, nrow = length(pcs), ncol = length(vars),
                  dimnames = list(pcs, vars))
  for (pc in pcs) for (v in vars) {
    # Continuous covariates enter the model as-is; categorical ones as factors.
    # Coercing a continuous variable to a factor would give it n-1 dummy columns
    # and again trivially inflate R^2, so we branch on type.
    xv <- md[[v]]
    predictor <- if (is.numeric(xv)) xv else factor(xv)
    fit <- stats::lm(scores[[pc]] ~ predictor)
    assoc[pc, v] <- summary(fit)$r.squared
  }

  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  tbl_path <- file.path(outdir, "qc", "variance_vs_metadata.tsv")
  dir.create(dirname(tbl_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(round(assoc, 3), tbl_path, sep = "\t", quote = FALSE, col.names = NA)

  png_path <- file.path(plots_dir, "variance_vs_metadata.png")
  # A heatmap of R^2; pheatmap handles the single-row/col edge cases of tiny data.
  if (length(vars) >= 1) {
    pheatmap::pheatmap(
      assoc, cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = TRUE, number_format = "%.2f",
      color = grDevices::colorRampPalette(c("white", "firebrick"))(100),
      main = "Variance explained (R^2) of each PC by each metadata variable",
      filename = png_path, width = 6, height = 4)
  }
  log_success(sprintf("Variance-vs-metadata written to %s and %s", tbl_path, png_path))
  list(assoc = assoc, plot = png_path, table = tbl_path)
}
