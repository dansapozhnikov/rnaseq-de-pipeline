# =============================================================================
# R/cache.R -- content-addressed checkpoint cache for resumable runs
# -----------------------------------------------------------------------------
# One concern: make the pipeline resumable. Each expensive, side-effect-free
# stage (DESeq fit, VST, per-contrast shrinkage+annotation, enrichment) is
# memoized to an .rds file under <outdir>/cache/, keyed by a hash of everything
# that determines its result: the input data, the relevant config, the upstream
# stage's key, AND a fingerprint of the compute code itself. A re-run (after a
# crash, an HPC pre-emption, or a manual restart) reloads finished stages
# instantly and only recomputes what actually changed.
#
# CORRECTNESS OVER REUSE: the key includes a deparse-hash of the compute
# functions, so editing the code -- committed or not -- invalidates the cache.
# We would rather recompute than ever serve a stale result.
# =============================================================================

# Hash an arbitrary R object to a short hex string. Uses `digest` when present,
# else falls back to serialize()+md5 (no extra dependency required).
.cache_hash <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(substr(digest::digest(x, algo = "xxhash64"), 1, 16))
  }
  tf <- tempfile(); on.exit(unlink(tf))
  saveRDS(x, tf)
  substr(unname(tools::md5sum(tf)), 1, 16)
}

#' Fingerprint the compute functions so a code change busts the cache.
#' Deparses the bodies of every function whose output we cache; if any changes,
#' the resulting hash changes and all downstream cache entries are invalidated.
cache_code_fingerprint <- function() {
  fns <- c("build_dds", "filter_low_counts", "run_deseq", "get_vst",
           "run_results", "shrink_lfc", "annotate_results",
           "resolve_result_specs", "run_enrichment")
  bodies <- lapply(fns, function(f) {
    if (exists(f, mode = "function")) deparse(get(f, mode = "function")) else f
  })
  .cache_hash(bodies)
}

#' Create a cache handle for a run.
#'
#' @param outdir Output root; cache lives under <outdir>/cache/ (override via dir).
#' @param enable If FALSE, caching is a no-op (every stage computes fresh).
#' @param version Pipeline version string (folded into every key).
#' @param base_deps A list of run-wide inputs (counts, metadata, design, seed,
#'   code fingerprint) that invalidate ALL stages when they change.
#' @param dir Optional explicit cache directory.
#' @return A `pipeline_cache` handle.
new_cache <- function(outdir, enable = TRUE, version = "dev",
                      base_deps = NULL, dir = NULL) {
  cache_dir <- dir %||% file.path(outdir, "cache")
  if (isTRUE(enable)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  structure(list(dir = cache_dir, enable = isTRUE(enable), version = version,
                 base = .cache_hash(list(version, base_deps))),
            class = "pipeline_cache")
}

#' Read the cache key stamped on a value returned by with_cache().
key_of <- function(obj) attr(obj, ".cache_key", exact = TRUE)

#' Memoized execution of one pipeline stage.
#'
#' @param cache A handle from new_cache().
#' @param stage Short stage name (used in the cache filename + logs).
#' @param deps Stage-specific inputs that determine the result.
#' @param compute A zero-arg function that produces the (serializable) result.
#' @param parent Optional upstream cache key, so upstream changes cascade.
#' @param verify Optional function(obj) -> logical; if it returns FALSE on a cache
#'   hit, the entry is treated as invalid and recomputed (used to re-check that a
#'   stage's side-effect output files still exist).
#' @return The stage result, with its cache key stamped as attr ".cache_key".
with_cache <- function(cache, stage, deps, compute, parent = NULL, verify = NULL) {
  if (!isTRUE(cache$enable)) {
    obj <- compute()
    attr(obj, ".cache_key") <- NA_character_
    return(obj)
  }
  key <- .cache_hash(list(cache$base, parent, stage, deps))
  path <- file.path(cache$dir, sprintf("%s.%s.rds", stage, key))

  if (file.exists(path)) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    ok <- !is.null(obj) && (is.null(verify) || isTRUE(tryCatch(verify(obj), error = function(e) FALSE)))
    if (ok) {
      log_success(sprintf("[resume] '%s' loaded from cache (skipped recompute).", stage))
      attr(obj, ".cache_key") <- key
      return(obj)
    }
    log_warn(sprintf("[resume] cached '%s' was invalid/stale; recomputing.", stage))
  }

  obj <- compute()
  # Drop any stale cache files for this stage (different key) to bound growth.
  stale <- list.files(cache$dir, pattern = sprintf("^%s\\.[0-9a-f]+\\.rds$", stage),
                      full.names = TRUE)
  unlink(setdiff(stale, path))
  tryCatch(saveRDS(obj, path),
           error = function(e) log_warn(sprintf("[cache] write failed for '%s': %s",
                                                 stage, conditionMessage(e))))
  log_info(sprintf("[cache] '%s' computed and cached.", stage))
  attr(obj, ".cache_key") <- key
  obj
}

# NULL-coalescing helper (also defined in enrich.R; harmless if sourced twice).
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
