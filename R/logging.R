# =============================================================================
# R/logging.R -- structured, colored, tee'd logging + severity helpers
# -----------------------------------------------------------------------------
# One concern: everything the pipeline says to the human. Messages are
# timestamped, written with a severity level, printed to the console in COLOR
# (via cli/crayon) and simultaneously appended to a per-run log file so a run is
# fully reconstructable after the fact.
# =============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(crayon)
})

# Package-local mutable state: the path of the current run's log file. Kept in an
# environment (not a global) so it is encapsulated but writable from log_*().
# WHY: we want a single append target for the whole run without threading a file
# handle through every function signature.
.log_state <- new.env(parent = emptyenv())
.log_state$file <- NULL

# NULL-coalescing helper, defined here (sourced first) so every module can use it
# regardless of source order. `a %||% b` returns `a` unless it is NULL, else `b`.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Initialise the run log
#'
#' @param outdir Root output directory; the log is written under <outdir>/logs/.
#' @param timestamp A filesystem-safe timestamp string identifying this run.
#' @return The full path to the created log file (invisibly).
#' Opens (creates) results/logs/run_<timestamp>.log and records it as the tee
#' target. WHY here: every subsequent log line is mirrored to this file so the
#' terminal session is not the only record of what happened.
init_logging <- function(outdir, timestamp) {
  log_dir <- file.path(outdir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  .log_state$file <- file.path(log_dir, sprintf("run_%s.log", timestamp))
  # Touch the file so downstream appends never fail on a missing path.
  cat("", file = .log_state$file)
  invisible(.log_state$file)
}

# Internal: format one line as "<ISO-timestamp> [LEVEL] message" and append it to
# the run log if one has been initialised. Stripping color codes keeps the file
# plain-text and greppable.
.log_to_file <- function(level, msg) {
  if (is.null(.log_state$file)) return(invisible())
  # Use a fixed, sortable timestamp; avoids locale-dependent formats.
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("%s [%s] %s", stamp, level, crayon::strip_style(msg))
  cat(line, "\n", file = .log_state$file, append = TRUE, sep = "")
}

#' Informational message (green arrow). Routine progress.
log_info <- function(msg) {
  cli::cli_alert_info(msg)
  .log_to_file("INFO", msg)
  invisible()
}

#' Success message (green tick). A step completed as intended.
log_success <- function(msg) {
  cli::cli_alert_success(msg)
  .log_to_file("OK", msg)
  invisible()
}

#' Warning message (amber). Proceed, but the human should see this.
log_warn <- function(msg) {
  cli::cli_alert_warning(crayon::yellow(msg))
  .log_to_file("WARN", msg)
  invisible()
}

#' Error message (red). Something is wrong; usually followed by stop_pipeline().
log_error <- function(msg) {
  cli::cli_alert_danger(crayon::red(msg))
  .log_to_file("ERROR", msg)
  invisible()
}

#' A visually prominent section header, so long console output stays navigable.
log_section <- function(title) {
  cli::cli_h1(title)
  .log_to_file("SECTION", paste0("=== ", title, " ==="))
  invisible()
}

#' Abort the run with an actionable message and a non-zero exit status.
#'
#' @param msg Human-facing explanation of what went wrong and how to fix it.
#' WHY a dedicated function: we want EVERY fatal exit to (a) be logged at ERROR
#' level to the tee file and (b) return a non-zero status so callers / CI can
#' detect failure. `quit(status = 1)` guarantees the shell sees the failure.
stop_pipeline <- function(msg) {
  log_error(msg)
  # In interactive sessions (or test mode) stop() is friendlier and, crucially,
  # is CATCHABLE by testthat's expect_error(); under a real Rscript run we must
  # quit non-zero so `--validate` and CI see the failure. Tests set
  # options(pipeline.test_mode = TRUE).
  if (interactive() || isTRUE(getOption("pipeline.test_mode"))) {
    stop(msg, call. = FALSE)
  } else {
    quit(save = "no", status = 1L)
  }
}

#' Resolve the pipeline version string: contents of the VERSION file, annotated
#' with the short git SHA when the project is a git checkout.
#' WHY: stamping the exact code version on every report/log makes a result
#' traceable to the commit that produced it.
get_pipeline_version <- function(project_root = ".") {
  v <- tryCatch(trimws(readLines(file.path(project_root, "VERSION"), warn = FALSE)[1]),
                error = function(e) "dev")
  sha <- tryCatch(
    suppressWarnings(system2("git", c("-C", project_root, "rev-parse", "--short", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(sha) && nzchar(sha[1])) paste0(v, " (", sha[1], ")") else v
}

#' Record the full session (package versions, R build, platform) to the log.
#' WHY: reproducibility -- sessionInfo() is the ground-truth manifest of what
#' actually ran, complementing renv.lock. Called at the very end of a run.
log_session_info <- function() {
  if (is.null(.log_state$file)) return(invisible())
  cat("\n===== sessionInfo() =====\n", file = .log_state$file, append = TRUE)
  si <- utils::capture.output(utils::sessionInfo())
  cat(si, sep = "\n", file = .log_state$file, append = TRUE)
  cat("\n", file = .log_state$file, append = TRUE)
  log_info(sprintf("sessionInfo() captured to %s", .log_state$file))
  invisible()
}
