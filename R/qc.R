# =============================================================================
# R/qc.R — the loud PASS / WARN / FAIL quality-control engine
# -----------------------------------------------------------------------------
# One concern: run configurable QC checks, record each with a severity, print
# them to the console in COLOR as they happen, draw a boxed summary banner that
# a human cannot miss, persist the results to TSV for the report, and gate the
# run (any FAIL => stop + non-zero exit, unless --force downgrades it).
#
# Severity semantics:
#   WARN = proceed, but surface prominently.
#   FAIL = stop the run (exit 1) unless --force, which downgrades FAIL -> WARN.
# =============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(crayon)
})

# ---------------------------------------------------------------------------
# QC accumulator
# ---------------------------------------------------------------------------

#' Create a fresh QC accumulator.
#' @return An environment holding a growing results data.frame. WHY an env: QC
#'   checks are added incrementally across the run (some before the model, some
#'   after DESeq()), so we need a mutable object that threads through by reference.
qc_init <- function() {
  qc <- new.env(parent = emptyenv())
  qc$results <- data.frame(
    check = character(), status = character(), severity = character(),
    value = character(), message = character(), stringsAsFactors = FALSE)
  qc$forced <- FALSE
  qc
}

#' Record one QC check and print it immediately in color.
#'
#' @param qc A QC accumulator from qc_init().
#' @param name Short check name (e.g. "library_size").
#' @param severity Severity if the check does NOT pass: "WARN" or "FAIL".
#' @param condition Logical; TRUE => PASS, FALSE => status becomes `severity`.
#' @param message Human explanation (shown on non-pass, and stored always).
#' @param value Optional observed value string for the record/report.
#' @return `qc` (invisibly), with one row appended.
#' The immediate colored print is what makes failures impossible to miss while
#' the run is happening; the stored row feeds the banner, the TSV, and the report.
qc_check <- function(qc, name, severity, condition, message = "", value = "") {
  severity <- match.arg(severity, c("WARN", "FAIL"))
  status <- if (isTRUE(condition)) "PASS" else severity

  qc$results <- rbind(qc$results, data.frame(
    check = name, status = status, severity = severity,
    value = as.character(value), message = message, stringsAsFactors = FALSE))

  # Loud, colored, per-check console line.
  label <- sprintf("%-26s %s", name, if (nzchar(value)) paste0("(", value, ")") else "")
  if (status == "PASS") {
    cli::cli_alert_success(crayon::green(paste0("PASS  ", label)))
  } else if (status == "WARN") {
    cli::cli_alert_warning(crayon::yellow(paste0("WARN  ", label, " - ", message)))
  } else {
    cli::cli_alert_danger(crayon::red(crayon::bold(paste0("FAIL  ", label, " - ", message))))
  }
  .log_to_file(status, sprintf("QC %s [%s] %s %s", name, status, value, message))
  invisible(qc)
}

# ---------------------------------------------------------------------------
# Individual checks (thresholds come from config; no magic numbers here)
# ---------------------------------------------------------------------------

#' Check 1 — count columns line up with metadata rows (count only, not identity;
#' identity is enforced separately by reconcile_samples()).
qc_check_sample_alignment <- function(qc, counts, metadata) {
  n_c <- ncol(counts); n_m <- nrow(metadata)
  qc_check(qc, "sample_alignment", "FAIL", n_c == n_m,
           sprintf("counts have %d samples but metadata has %d rows.", n_c, n_m),
           value = sprintf("%d vs %d", n_c, n_m))
}

#' Check 2 — raw-count sanity: integer and non-negative. The loaders already
#' enforce this, so a PASS here is a belt-and-braces confirmation for the report.
qc_check_raw_counts <- function(qc, counts) {
  ok <- is.integer(counts) && !any(counts < 0)
  qc_check(qc, "raw_count_sanity", "FAIL", ok,
           "counts must be non-negative integers (raw counts, not TPM/CPM).",
           value = if (ok) "integer, non-negative" else "NON-INTEGER/NEGATIVE")
}

#' Check 3 — per-sample library size. WARN below min, FAIL below the hard floor.
qc_check_library_size <- function(qc, counts, cfg) {
  libsize <- colSums(counts)
  worst <- min(libsize)
  worst_id <- names(libsize)[which.min(libsize)]
  # Hard floor first: a sample below this is unusable.
  qc_check(qc, "library_size_floor", "FAIL",
           worst >= cfg$qc$min_library_size_fail,
           sprintf("sample '%s' has only %s reads (hard floor %s).",
                   worst_id, format(worst, big.mark=","),
                   format(cfg$qc$min_library_size_fail, big.mark=",", scientific=FALSE)),
           value = sprintf("min=%s", format(worst, big.mark=",")))
  # Soft threshold: WARN.
  qc_check(qc, "library_size_min", "WARN",
           worst >= cfg$qc$min_library_size,
           sprintf("sample '%s' below recommended %s reads.",
                   worst_id, format(cfg$qc$min_library_size, big.mark=",", scientific=FALSE)),
           value = sprintf("min=%s", format(worst, big.mark=",")))
}

#' Check 4 — genes detected (count > 0) per sample. WARN if any sample is low.
qc_check_genes_detected <- function(qc, counts, cfg) {
  detected <- colSums(counts > 0)
  worst <- min(detected)
  worst_id <- names(detected)[which.min(detected)]
  qc_check(qc, "genes_detected", "WARN",
           worst >= cfg$qc$min_genes_detected,
           sprintf("sample '%s' detects only %d genes (min %d).",
                   worst_id, worst, cfg$qc$min_genes_detected),
           value = sprintf("min=%d", worst))
}

#' Check 5 — replicates per group in the TESTED factor (last term of the design).
#' FAIL if any group has fewer than the configured minimum: DE with <2 replicates
#' per group cannot estimate within-group dispersion.
qc_check_replicates <- function(qc, metadata, tested_factor, cfg) {
  if (!tested_factor %in% colnames(metadata)) {
    return(qc_check(qc, "replicates_per_group", "FAIL", FALSE,
                    sprintf("tested factor '%s' not found in metadata.", tested_factor),
                    value = "missing factor"))
  }
  tab <- table(metadata[[tested_factor]])
  worst <- min(tab)
  qc_check(qc, "replicates_per_group", "FAIL",
           worst >= cfg$qc$min_replicates_per_group,
           sprintf("group '%s' has %d replicate(s); need >= %d. Counts: %s",
                   names(tab)[which.min(tab)], worst, cfg$qc$min_replicates_per_group,
                   paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse=", ")),
           value = sprintf("min group n=%d", worst))
}

#' Check 6 — design identifiability. Build the model matrix and check it is full
#' rank. A rank-deficient matrix means a covariate is confounded with the tested
#' variable (e.g. every treated sample is batch B) and DESeq2 cannot separate
#' their effects. On FAIL we print the offending cross-tabulation.
qc_check_design_identifiable <- function(qc, metadata, design_str, tested_factor) {
  design <- stats::as.formula(design_str)
  # Coerce all design variables to factors so model.matrix builds contrasts.
  vars <- all.vars(design)
  md <- metadata
  for (v in vars) if (v %in% colnames(md)) md[[v]] <- factor(md[[v]])
  mm <- tryCatch(stats::model.matrix(design, data = md),
                 error = function(e) NULL)
  if (is.null(mm)) {
    return(qc_check(qc, "design_identifiable", "FAIL", FALSE,
                    sprintf("could not build model.matrix for design %s.", design_str),
                    value = "model.matrix error"))
  }
  full_rank <- qr(mm)$rank == ncol(mm)
  msg <- "design is full rank."
  if (!full_rank) {
    # Show the tested factor crossed with each other design variable to reveal
    # the confound. Printed to the log so the user can see the culprit.
    others <- setdiff(vars, tested_factor)
    xtabs_txt <- vapply(others, function(o) {
      paste0(o, " x ", tested_factor, ":\n",
             paste(utils::capture.output(print(table(md[[o]], md[[tested_factor]]))),
                   collapse = "\n"))
    }, character(1))
    msg <- sprintf("design '%s' is rank-deficient (confounded covariate). %s",
                   design_str, paste(xtabs_txt, collapse = "\n"))
  }
  qc_check(qc, "design_identifiable", "FAIL", full_rank, msg,
           value = sprintf("rank %d / %d cols", qr(mm)$rank, ncol(mm)))
}

#' Check 7 — low-count gene filter report. WARN if the filter removes more than
#' the configured percentage of genes (a huge fraction removed can signal a
#' shallow library or a mismatched annotation).
qc_check_low_count_filter <- function(qc, n_before, n_after, cfg) {
  pct <- 100 * (n_before - n_after) / n_before
  qc_check(qc, "low_count_filter", "WARN",
           pct <= cfg$qc$max_pct_genes_filtered,
           sprintf("low-count filter removed %.1f%% of genes (> %d%% threshold).",
                   pct, cfg$qc$max_pct_genes_filtered),
           value = sprintf("%d -> %d genes (%.1f%% removed)", n_before, n_after, pct))
}

#' Check 8 — dispersion-fit sanity (post-DESeq). WARN if the fitted gene-wise
#' dispersions look abnormal (all near-zero, or an implausibly large median),
#' which usually indicates too few replicates or a mis-specified design.
qc_check_dispersion <- function(qc, dds) {
  disp <- tryCatch(DESeq2::dispersions(dds), error = function(e) NA_real_)
  med <- suppressWarnings(stats::median(disp, na.rm = TRUE))
  # Heuristic sane band for bulk RNA-seq gene-wise dispersion medians.
  ok <- is.finite(med) && med > 1e-4 && med < 10
  qc_check(qc, "dispersion_fit", "WARN", ok,
           sprintf("median dispersion %.4g is outside the expected 1e-4..10 band; check replicates/design.", med),
           value = sprintf("median disp=%.4g", med))
}

#' Check 9 — outlier flag from Cook's distance (post-DESeq). WARN and LIST the
#' candidate samples. We do NOT auto-remove; dropping samples requires explicit
#' opt-in (documented in USAGE), because silent removal changes the model.
qc_check_cooks_outliers <- function(qc, dds, cfg) {
  if (!isTRUE(cfg$qc$cooks_outlier)) return(invisible(qc))
  cooks <- tryCatch(SummarizedExperiment::assays(dds)[["cooks"]],
                    error = function(e) NULL)
  if (is.null(cooks)) {
    return(qc_check(qc, "cooks_outliers", "WARN", TRUE,
                    "Cook's distances unavailable (skipped).", value = "n/a"))
  }
  # A sample with many genes above the F-based Cook's cutoff is a candidate.
  m <- ncol(dds); p <- length(DESeq2::resultsNames(dds))
  cutoff <- stats::qf(0.99, p, m - p)
  frac_high <- colMeans(cooks > cutoff, na.rm = TRUE)
  flagged <- names(frac_high)[frac_high > 0.01]  # >1% of genes flagged
  qc_check(qc, "cooks_outliers", "WARN", length(flagged) == 0,
           sprintf("candidate outlier sample(s): %s (many high-Cook's genes). Not removed automatically.",
                   paste(flagged, collapse = ", ")),
           value = if (length(flagged)) paste(flagged, collapse=",") else "none")
}

# ---------------------------------------------------------------------------
# Summary banner, persistence, and the gate
# ---------------------------------------------------------------------------

#' Draw a boxed, colored summary banner listing every WARN and FAIL.
#' WHY a hand-drawn box: guarantees a visually unmissable block at the end of QC
#' regardless of terminal width quirks; green if all clear, else lists problems.
qc_summary <- function(qc) {
  res <- qc$results
  n_fail <- sum(res$status == "FAIL")
  n_warn <- sum(res$status == "WARN")
  n_pass <- sum(res$status == "PASS")

  lines <- c(sprintf("QC SUMMARY:  %d PASS   %d WARN   %d FAIL", n_pass, n_warn, n_fail))
  probs <- res[res$status != "PASS", , drop = FALSE]
  if (nrow(probs)) {
    lines <- c(lines, "")
    for (i in seq_len(nrow(probs))) {
      lines <- c(lines, sprintf("  [%s] %s: %s", probs$status[i], probs$check[i], probs$message[i]))
    }
  }
  # Box drawing.
  width <- max(nchar(unlist(strsplit(lines, "\n")))) + 2
  top <- paste0("+", strrep("-", width), "+")
  colorize <- if (n_fail > 0) crayon::red else if (n_warn > 0) crayon::yellow else crayon::green

  cat("\n")
  cat(colorize(crayon::bold(top)), "\n")
  for (ln in lines) {
    for (sub in strsplit(ln, "\n")[[1]]) {
      cat(colorize(crayon::bold(sprintf("| %-*s |", width - 2, sub))), "\n")
    }
  }
  cat(colorize(crayon::bold(top)), "\n\n")
  .log_to_file("SUMMARY", paste(lines, collapse = " | "))
  invisible(qc)
}

#' Write the QC results table to <outdir>/qc/qc_results.tsv for the report.
qc_write <- function(qc, outdir) {
  qc_dir <- file.path(outdir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(qc_dir, "qc_results.tsv")
  utils::write.table(qc$results, path, sep = "\t", quote = FALSE, row.names = FALSE)
  log_info(sprintf("QC results written to %s", path))
  invisible(path)
}

#' The gate: decide whether the run may continue.
#'
#' @param qc QC accumulator.
#' @param force If TRUE, downgrade FAIL -> WARN and continue (recorded).
#' Any FAIL stops the run with exit status 1 (via stop_pipeline). With --force we
#' proceed but stamp the record so the report shows the run was forced past a FAIL.
qc_gate <- function(qc, force = FALSE) {
  fails <- qc$results[qc$results$status == "FAIL", , drop = FALSE]
  if (nrow(fails) == 0) return(invisible(qc))

  if (force) {
    qc$forced <- TRUE
    qc$results$status[qc$results$status == "FAIL"] <- "WARN"
    log_warn(sprintf("--force: downgraded %d FAIL(s) to WARN and continuing. This is recorded.",
                     nrow(fails)))
    return(invisible(qc))
  }
  stop_pipeline(sprintf(
    "QC gate: %d check(s) FAILED (%s). Fix the input/design, or re-run with --force to override.",
    nrow(fails), paste(fails$check, collapse = ", ")))
}
