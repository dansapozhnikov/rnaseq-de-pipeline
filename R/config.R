# =============================================================================
# R/config.R — load, validate, and echo the YAML configuration
# -----------------------------------------------------------------------------
# One concern: turn a config.yaml file into a validated R list, fail early and
# clearly on anything missing or malformed, and echo every resolved value so the
# log is a complete record of the parameters a run used.
# =============================================================================

suppressPackageStartupMessages(library(yaml))

# The set of fields the pipeline cannot run without. Expressed as dotted paths
# into the nested list. WHY explicit: a typo'd or absent key should produce a
# named, actionable error ("missing required config field: de.padj_cutoff")
# rather than a cryptic NULL error deep in the DE step.
.required_config_fields <- c(
  "seed", "organism", "input_format", "design", "contrast",
  "qc.min_library_size", "qc.min_library_size_fail", "qc.min_genes_detected",
  "qc.min_replicates_per_group", "qc.min_count", "qc.max_pct_genes_filtered",
  "de.padj_cutoff", "de.lfc_cutoff", "de.shrink",
  "enrichment.enable", "enrichment.ontology",
  "outputs.outdir"
)

# Internal: fetch a dotted path (e.g. "qc.min_count") from a nested list.
# Returns NULL if any level is absent.
.get_path <- function(x, path) {
  for (key in strsplit(path, ".", fixed = TRUE)[[1]]) {
    if (!is.list(x) || is.null(x[[key]])) return(NULL)
    x <- x[[key]]
  }
  x
}

#' Load and validate a pipeline configuration file.
#'
#' @param path Path to a config.yaml file.
#' @return A validated, normalised config list.
#' Reads the YAML, checks every required field is present, and normalises a few
#' values (organism casing, contrast length). WHY validate up front: we would
#' rather stop at second 1 with "you forgot de.shrink" than crash 20 minutes in.
load_config <- function(path) {
  if (!file.exists(path)) {
    stop_pipeline(sprintf(
      "Config file not found: '%s'. Copy config/config.example.yaml and edit it.",
      path))
  }
  cfg <- yaml::read_yaml(path)

  # Check required fields exist (NULL = absent).
  missing <- .required_config_fields[
    vapply(.required_config_fields, function(f) is.null(.get_path(cfg, f)), logical(1))
  ]
  if (length(missing)) {
    stop_pipeline(sprintf(
      "Config is missing required field(s): %s. See config/config.example.yaml.",
      paste(missing, collapse = ", ")))
  }

  # Normalise organism to a known value; anything else is a hard error because it
  # decides which annotation DB we load later.
  cfg$organism <- tolower(trimws(cfg$organism))
  if (!cfg$organism %in% c("human", "mouse")) {
    stop_pipeline(sprintf(
      "config$organism must be 'human' or 'mouse', got '%s'.", cfg$organism))
  }

  # The contrast must be exactly [factor, numerator, denominator].
  if (length(cfg$contrast) != 3L) {
    stop_pipeline(sprintf(
      "config$contrast must have 3 elements [factor, numerator, denominator], got %d.",
      length(cfg$contrast)))
  }

  cfg
}

#' Echo the resolved configuration to the log.
#'
#' @param cfg A config list from load_config().
#' Prints every parameter the run will use. WHY: the log becomes a self-contained
#' record — anyone reading it can see exactly how the run was parameterised
#' without needing the original YAML.
echo_config <- function(cfg) {
  log_section("Resolved configuration")
  log_info(sprintf("seed              : %s", cfg$seed))
  log_info(sprintf("organism          : %s", cfg$organism))
  log_info(sprintf("input_format      : %s", cfg$input_format))
  log_info(sprintf("design            : %s", cfg$design))
  log_info(sprintf("contrast          : %s (numerator) vs %s (denominator) on '%s'",
                   cfg$contrast[[2]], cfg$contrast[[3]], cfg$contrast[[1]]))
  log_info(sprintf("qc.min_library_size      : %s (FAIL floor %s)",
                   cfg$qc$min_library_size, cfg$qc$min_library_size_fail))
  log_info(sprintf("qc.min_genes_detected    : %s", cfg$qc$min_genes_detected))
  log_info(sprintf("qc.min_replicates_group  : %s", cfg$qc$min_replicates_per_group))
  log_info(sprintf("qc.min_count filter      : %s", cfg$qc$min_count))
  log_info(sprintf("qc.max_pct_genes_filtered: %s%%", cfg$qc$max_pct_genes_filtered))
  log_info(sprintf("de.padj_cutoff / lfc     : %s / %s",
                   cfg$de$padj_cutoff, cfg$de$lfc_cutoff))
  log_info(sprintf("de.shrink                : %s", cfg$de$shrink))
  log_info(sprintf("enrichment               : enable=%s ontology=%s",
                   cfg$enrichment$enable, cfg$enrichment$ontology))
  log_info(sprintf("outputs.outdir           : %s", cfg$outputs$outdir))
  invisible()
}
