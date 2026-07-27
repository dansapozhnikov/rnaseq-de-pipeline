# =============================================================================
# R/report.R -- assemble the self-contained HTML report
# -----------------------------------------------------------------------------
# One concern: render report/report_template.Rmd against a finished run's output
# directory, passing the parameters the template needs. The template reads the
# artifacts (qc_results.tsv, plots, tables) from <outdir>, so this function's job
# is to locate the template, pass config + run metadata, and render to HTML.
# =============================================================================

#' Render the HTML report for a completed run.
#'
#' @param cfg The resolved config list.
#' @param outdir Output root that already contains qc/, plots/, tables/.
#' @param run_timestamp The run's timestamp string (for provenance in the report).
#' @param batch_vars Character vector of batch variables used (for the report).
#' @param forced Logical; TRUE if --force downgraded a QC FAIL.
#' @param template Path to the Rmd template (defaults to report/report_template.Rmd).
#' @return Path to the rendered HTML (invisibly), or NULL on failure (WARN).
#' Rendering is wrapped so a report hiccup never masks an otherwise-successful
#' analysis -- the tables and plots are already on disk regardless.
render_report <- function(cfg, outdir, run_timestamp, batch_vars = character(0),
                          forced = FALSE, template = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    log_warn("rmarkdown not installed; skipping HTML report."); return(invisible(NULL))
  }
  if (is.null(template)) {
    # Resolve relative to this file's project root, robust to the working dir.
    template <- file.path("report", "report_template.Rmd")
  }
  if (!file.exists(template)) {
    log_warn(sprintf("Report template not found at %s; skipping.", template)); return(invisible(NULL))
  }

  out_dir_abs <- normalizePath(outdir, mustWork = FALSE)
  dir.create(out_dir_abs, recursive = TRUE, showWarnings = FALSE)
  contrast_str <- sprintf("%s: %s vs %s", cfg$contrast[[1]], cfg$contrast[[2]], cfg$contrast[[3]])
  # Render with explicit output_dir + a throwaway intermediates_dir. WHY: the Rmd
  # lives under report/, and without these rmarkdown resolves output/intermediate
  # paths relative to the template dir (or the shifting knit working dir), which
  # fails when outdir is a relative path. Absolute dirs make it location-proof.
  intermediates <- file.path(tempdir(), paste0("rmd_", run_timestamp))

  res <- tryCatch({
    out <- rmarkdown::render(
      input = normalizePath(template),
      output_file = "report.html",
      output_dir = out_dir_abs,
      intermediates_dir = intermediates,
      knit_root_dir = out_dir_abs,
      params = list(
        outdir = out_dir_abs,
        cfg = cfg,
        forced = forced,
        run_timestamp = run_timestamp,
        batch_vars = paste(batch_vars, collapse = ", "),
        contrast = contrast_str),
      envir = new.env(parent = globalenv()),  # isolate the render environment
      quiet = TRUE)
    file.path(out_dir_abs, "report.html")
  }, error = function(e) {
    log_warn(sprintf("Report rendering failed: %s", conditionMessage(e))); NULL
  })

  if (!is.null(res)) log_success(sprintf("HTML report rendered to %s", res))
  invisible(res)
}
