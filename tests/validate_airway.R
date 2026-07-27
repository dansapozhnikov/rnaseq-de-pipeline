#!/usr/bin/env Rscript
# =============================================================================
# tests/validate_airway.R -- positive-control self-test on the Bioconductor
# `airway` dataset (dexamethasone-treated airway smooth-muscle cells).
# -----------------------------------------------------------------------------
# WHY this exists: `airway` has a well-characterised glucocorticoid response, so
# it is a KNOWN-ANSWER end-to-end test. It runs the WHOLE pipeline through the
# same run_pipeline_core() users run, then asserts the pipeline recovers the
# expected biology: CRISPLD2, DUSP1 and KLF15 are significant and UP in the
# treated group. This is both the CI smoke test and the first green E2E run.
# =============================================================================

#' Run the airway end-to-end validation.
#'
#' @return TRUE if all biological assertions pass, FALSE otherwise.
validate_airway <- function(cfg_path, outdir, run_timestamp = NULL,
                            force = FALSE, project_root = ".",
                            pipeline_version = NULL, resume = TRUE) {
  # If invoked standalone (not via run_pipeline.R), source the modules first.
  if (!exists("run_pipeline_core", mode = "function")) {
    for (f in c("logging.R", "cache.R", "config.R", "io_detect.R", "qc.R", "normalize.R",
                "explore.R", "batch.R", "de.R", "enrich.R", "report.R", "pipeline.R")) {
      source(file.path(project_root, "R", f))
    }
  }
  if (is.null(run_timestamp)) run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (is.null(pipeline_version)) pipeline_version <- get_pipeline_version(project_root)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  init_logging(outdir, run_timestamp)
  log_section("AIRWAY VALIDATION (positive control)")

  if (!requireNamespace("airway", quietly = TRUE)) {
    log_error("Package 'airway' is not installed; cannot run validation."); return(FALSE)
  }

  # --- Load airway and coerce to the pipeline's inputs -----------------------
  env <- new.env(); utils::data("airway", package = "airway", envir = env)
  se <- env$airway
  counts <- SummarizedExperiment::assay(se)
  storage.mode(counts) <- "integer"
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  meta$sample <- rownames(meta)

  # --- Config: start from the example, override to the airway experiment -----
  cfg <- load_config(cfg_path)
  cfg$organism <- "human"
  cfg$design   <- "~ cell + dex"                 # cell line is the nuisance/batch term
  cfg$contrast <- list("dex", "trt", "untrt")    # treated vs untreated (untrt = reference)
  echo_config(cfg)

  # --- Run the WHOLE pipeline (same path as a real run) ----------------------
  result <- run_pipeline_core(counts, meta, cfg, outdir, run_timestamp,
                              sample_col = "sample", force = force,
                              pipeline_version = pipeline_version, resume = resume)

  # ===========================================================================
  # Assertions on the recovered biology
  # ===========================================================================
  log_section("Validation assertions")
  tab <- result$de$table
  all_ok <- TRUE

  # (1) Known glucocorticoid-responsive genes: significant AND up in treated.
  targets <- c("CRISPLD2", "DUSP1", "KLF15")
  for (sym in targets) {
    row <- tab[!is.na(tab$symbol) & tab$symbol == sym, , drop = FALSE]
    if (nrow(row) == 0) {
      log_error(sprintf("ASSERT FAIL: %s not found in results.", sym)); all_ok <- FALSE; next
    }
    r <- row[1, ]
    sig <- !is.na(r$padj) && r$padj < cfg$de$padj_cutoff
    up  <- !is.na(r$log2FoldChange) && r$log2FoldChange > 0
    if (sig && up) {
      log_success(sprintf("ASSERT PASS: %s significant & UP in treated (log2FC=%.2f, padj=%.2g).",
                         sym, r$log2FoldChange, r$padj))
    } else {
      log_error(sprintf("ASSERT FAIL: %s expected significant & up; got log2FC=%.2f, padj=%.2g.",
                        sym, r$log2FoldChange, r$padj)); all_ok <- FALSE
    }
  }

  # (2) Overall DEG count in a sane range (guards against a broken model that
  #     calls everything or nothing significant).
  n_sig <- result$de$n_sig
  sane <- n_sig >= 50 && n_sig <= 15000
  if (sane) {
    log_success(sprintf("ASSERT PASS: DEG count %d is in the sane range [50, 15000].", n_sig))
  } else {
    log_error(sprintf("ASSERT FAIL: DEG count %d is outside the sane range [50, 15000].", n_sig))
    all_ok <- FALSE
  }

  log_session_info()
  if (all_ok) {
    cat(crayon::bold(crayon::green("\n==== AIRWAY VALIDATION: PASS ====\n\n")))
  } else {
    cat(crayon::bold(crayon::red("\n==== AIRWAY VALIDATION: FAIL ====\n\n")))
  }
  all_ok
}

# If executed directly (Rscript tests/validate_airway.R [config] [outdir]), run it.
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (!interactive() && length(grep("--file=validate_airway|validate_airway.R", commandArgs(FALSE)))) {
    root <- tryCatch(dirname(dirname(normalizePath(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), error = function(e) ".")
    cfg  <- if (length(.args) >= 1) .args[1] else file.path(root, "config", "config.example.yaml")
    od   <- if (length(.args) >= 2) .args[2] else file.path("results", "airway_validation")
    ok <- validate_airway(cfg_path = cfg, outdir = od, project_root = root)
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }
}
