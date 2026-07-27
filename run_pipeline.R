#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R -- command-line orchestrator
# -----------------------------------------------------------------------------
# Parses CLI arguments, loads config + counts (or the airway self-test), runs the
# shared core, captures sessionInfo(), and exits non-zero on any QC FAIL or on a
# --validate assertion failure.
#
# Usage:
#   Rscript run_pipeline.R --counts <file|dir> --metadata sample_sheet.csv \
#     --config config/config.yaml [--input-format auto] [--tx2gene tx2gene.csv] \
#     [--outdir results/] [--sample-col sample] [--star-strand 2] [--force]
#   Rscript run_pipeline.R --validate --config config/config.example.yaml
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# --- Resolve the directory this script lives in, so sourcing R/ works from any
#     working directory (e.g. invoked from a parent dir or CI).
.this_file <- (function() {
  ca <- commandArgs(FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  normalizePath("run_pipeline.R", mustWork = FALSE)
})()
PROJECT_ROOT <- dirname(.this_file)

# Source modules in dependency order (logging first: defines stop_pipeline).
for (f in c("logging.R", "config.R", "io_detect.R", "qc.R", "normalize.R",
            "explore.R", "batch.R", "de.R", "enrich.R", "report.R", "pipeline.R")) {
  source(file.path(PROJECT_ROOT, "R", f))
}

# ---------------------------------------------------------------------------
# CLI definition
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--counts", type = "character", default = NULL,
              help = "Path to counts (file or directory). Ignored with --validate."),
  make_option("--metadata", type = "character", default = NULL,
              help = "Path to the sample-sheet CSV. Ignored with --validate."),
  make_option("--config", type = "character", default = NULL,
              help = "Path to config.yaml (required)."),
  make_option("--input-format", dest = "input_format", type = "character", default = "auto",
              help = "auto|matrix|featurecounts|star|htseq|salmon|rangedSE [default %default]."),
  make_option("--tx2gene", type = "character", default = NULL,
              help = "transcript->gene map CSV (required for salmon/kallisto)."),
  make_option("--outdir", type = "character", default = NULL,
              help = "Output directory [default: outputs.outdir from config]."),
  make_option("--sample-col", dest = "sample_col", type = "character", default = "sample",
              help = "Metadata column with sample IDs [default %default]."),
  make_option("--star-strand", dest = "star_strand", type = "integer", default = 2L,
              help = "STAR count column: 2=unstranded, 3=fwd, 4=rev [default %default]."),
  make_option("--force", action = "store_true", default = FALSE,
              help = "Downgrade QC FAIL -> WARN and continue (logged)."),
  make_option("--validate", action = "store_true", default = FALSE,
              help = "Run the airway positive-control self-test end-to-end.")
)
opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Reproducible bulk RNA-seq differential-expression pipeline."))

# A single timestamp identifies this run across the log, outputs, and report.
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---------------------------------------------------------------------------
# Validate mode: hand off to the airway self-test (its own asserts + exit code).
# ---------------------------------------------------------------------------
# Resolve the pipeline version once (VERSION file + git SHA) for provenance.
pipeline_version <- get_pipeline_version(PROJECT_ROOT)

if (isTRUE(opt$validate)) {
  source(file.path(PROJECT_ROOT, "tests", "validate_airway.R"))
  cfg_path <- opt$config %||% file.path(PROJECT_ROOT, "config", "config.example.yaml")
  outdir <- opt$outdir %||% file.path("results", "airway_validation")
  ok <- validate_airway(cfg_path = cfg_path, outdir = outdir,
                        run_timestamp = run_timestamp, force = opt$force,
                        project_root = PROJECT_ROOT, pipeline_version = pipeline_version)
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

# ---------------------------------------------------------------------------
# Normal mode
# ---------------------------------------------------------------------------
if (is.null(opt$config)) stop("--config is required (see config/config.example.yaml).")
if (is.null(opt$counts)) stop("--counts is required (or use --validate).")
if (is.null(opt$metadata)) stop("--metadata is required (or use --validate).")

cfg <- load_config(opt$config)
outdir <- opt$outdir %||% cfg$outputs$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Start logging (tee to results/logs/run_<ts>.log) and echo resolved config.
init_logging(outdir, run_timestamp)
log_section("RNA-seq DE pipeline")
log_info(sprintf("Run timestamp: %s", run_timestamp))
echo_config(cfg)

# Load counts (auto-detect or forced) and the sample sheet.
log_section("Input")
loaded <- load_counts(opt$counts, input_format = opt$input_format,
                      tx2gene = opt$tx2gene, star_strand = opt$star_strand)
counts <- loaded$counts

# rangedSE inputs may carry their own colData; otherwise read the sample sheet.
metadata <- if (!is.null(loaded$coldata) && is.null(opt$metadata)) {
  md <- as.data.frame(loaded$coldata); md[[opt$sample_col]] <- rownames(md); md
} else {
  read.csv(opt$metadata, stringsAsFactors = FALSE, check.names = FALSE)
}

# Run everything.
result <- run_pipeline_core(counts, metadata, cfg, outdir, run_timestamp,
                            sample_col = opt$sample_col, tx2gene = opt$tx2gene,
                            force = opt$force, pipeline_version = pipeline_version)

# Capture sessionInfo() to the log for reproducibility, then finish.
log_session_info()
log_success(sprintf("Done. %d significant genes. Report: %s",
                    result$de$n_sig, result$report %||% "(not rendered)"))
quit(save = "no", status = 0L)
