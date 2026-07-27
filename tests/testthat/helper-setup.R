# Shared test setup: source the pipeline's R/ modules and enable test mode so
# stop_pipeline() throws a catchable error instead of quitting the process.
options(pipeline.test_mode = TRUE)

# Locate the repo root relative to this helper (tests/testthat/) so tests work
# whether launched from the repo root or via R CMD check.
.repo_root <- normalizePath(file.path(dirname(dirname(getwd()))), mustWork = FALSE)
if (!dir.exists(file.path(.repo_root, "R"))) {
  # Fallback: testthat often sets wd to tests/testthat already.
  .repo_root <- normalizePath("../..", mustWork = FALSE)
}

fixture <- function(...) file.path(.repo_root, "tests", "fixtures", ...)

# Source modules under test (order matters: logging first, it defines stop_pipeline).
for (f in c("logging.R", "config.R", "io_detect.R", "qc.R")) {
  p <- file.path(.repo_root, "R", f)
  if (file.exists(p)) sys.source(p, envir = globalenv())
}
