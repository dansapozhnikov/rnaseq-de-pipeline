# =============================================================================
# R/io_detect.R -- input-format detection + loaders for the common bulk formats
# -----------------------------------------------------------------------------
# One concern: turn whatever the user points --counts at (a file or a directory)
# into a clean integer count matrix (genes x samples) plus optional feature
# metadata, and reconcile it against the sample sheet. Every loader returns the
# SAME shape so the rest of the pipeline is format-agnostic:
#     list(counts = <integer matrix, genes x samples>,
#          feature_meta = <data.frame or NULL>)
# We NEVER silently accept normalized data (TPM/CPM/FPKM): DESeq2 models raw
# count overdispersion, so non-integer "counts" are a correctness bug, not a
# rounding nuisance. We FAIL loudly and tell the user to supply raw counts.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
})

# ---------------------------------------------------------------------------
# Coercion guardrail
# ---------------------------------------------------------------------------

#' Coerce a numeric matrix to a validated integer count matrix.
#'
#' @param mat Numeric matrix (genes x samples).
#' @param source_label Human label for messages (e.g. "matrix", "salmon").
#' @param allow_round If TRUE, fractional values are expected (Salmon estimated
#'   counts) and are rounded with a logged note. If FALSE, fractional values are
#'   treated as evidence the input is normalized and trigger a FAIL.
#' @return Integer matrix.
#' WHY: this is the single chokepoint that protects DESeq2 from non-count input.
.coerce_integer_matrix <- function(mat, source_label, allow_round = FALSE) {
  storage.mode(mat) <- "double"

  # Negative values are never valid raw counts.
  if (any(mat < 0, na.rm = TRUE)) {
    stop_pipeline(sprintf(
      "[%s] Count matrix contains negative values; these are not raw counts.",
      source_label))
  }

  # Fraction of entries that are not (near-)integers.
  non_int <- mean(abs(mat - round(mat)) > 1e-8, na.rm = TRUE)

  if (non_int > 0.001) {
    if (allow_round) {
      # Salmon/kallisto estimated counts are legitimately fractional; rounding to
      # integer is the documented way to feed them to DESeq2.
      log_warn(sprintf(
        "[%s] %.1f%% of values are fractional (estimated counts); rounding to integers.",
        source_label, 100 * non_int))
    } else {
      stop_pipeline(sprintf(paste0(
        "[%s] %.1f%% of values are non-integer. This looks like normalized data ",
        "(TPM/CPM/FPKM), not raw counts. DESeq2 requires RAW integer counts. ",
        "Please supply un-normalized counts, or use --input-format salmon for ",
        "transcript-quantifier output."),
        source_label, 100 * non_int))
    }
  }

  mat <- round(mat)
  storage.mode(mat) <- "integer"
  mat
}

# Internal: read a delimited table with data.table's auto-sep detection.
.fread <- function(path, ...) {
  data.table::fread(path, data.table = FALSE, ...)
}

# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

#' Auto-detect the count input format by inspecting the path.
#'
#' @param path A file or directory.
#' @return One of "matrix","featurecounts","star","htseq","salmon","rangedSE".
#' Logs which format was chosen AND why, so an auto-detection decision is always
#' auditable in the run log.
detect_format <- function(path) {
  if (!file.exists(path)) {
    stop_pipeline(sprintf("--counts path does not exist: %s", path))
  }

  # ---- Directory inputs: per-sample files from an aligner/quantifier ---------
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE)
    if (any(grepl("quant\\.sf$", files)) || any(grepl("abundance\\.tsv$", files))) {
      log_info("Detected format: salmon/kallisto (found quant.sf / abundance.tsv).")
      return("salmon")
    }
    if (any(grepl("ReadsPerGene\\.out\\.tab$", files))) {
      log_info("Detected format: star (found *ReadsPerGene.out.tab files).")
      return("star")
    }
    # HTSeq: 2-column per-sample files containing the tell-tale __no_feature row.
    cand <- files[grepl("\\.(txt|tsv|counts|tab)$", files)]
    if (length(cand) &&
        any(vapply(cand, function(f) {
          head_lines <- tryCatch(readLines(f, n = 50), error = function(e) character())
          any(grepl("^__no_feature", head_lines))
        }, logical(1)))) {
      log_info("Detected format: htseq (found per-sample files with __no_feature rows).")
      return("htseq")
    }
    stop_pipeline(sprintf(
      "Could not auto-detect a supported format in directory '%s'. Use --input-format.",
      path))
  }

  # ---- Single-file inputs ----------------------------------------------------
  # Serialized SummarizedExperiment (e.g. airway).
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    log_info("Detected format: rangedSE (.rds serialized SummarizedExperiment).")
    return("rangedSE")
  }

  # Peek at the first lines to distinguish featureCounts from a plain matrix.
  head_lines <- readLines(path, n = 2, warn = FALSE)
  # featureCounts writes a leading '# Program:featureCounts ...' comment line,
  # then a header beginning Geneid<TAB>Chr<TAB>Start<TAB>End<TAB>Strand<TAB>Length.
  if (length(head_lines) >= 2 &&
      grepl("^#", head_lines[1]) &&
      grepl("^Geneid\t.*\tLength", head_lines[2])) {
    log_info("Detected format: featurecounts (leading '#' + Geneid..Length header).")
    return("featurecounts")
  }
  if (grepl("^Geneid\t.*\tLength", head_lines[1])) {
    log_info("Detected format: featurecounts (Geneid..Length header, no comment line).")
    return("featurecounts")
  }

  log_info("Detected format: matrix (generic gene-by-sample table).")
  "matrix"
}

# ---------------------------------------------------------------------------
# Loaders (one per format) -- each returns list(counts, feature_meta)
# ---------------------------------------------------------------------------

#' Load a plain count matrix (CSV/TSV): first column = gene IDs, rest numeric.
load_matrix <- function(path) {
  df <- .fread(path)
  gene_ids <- as.character(df[[1]])            # first column holds gene identifiers
  mat <- as.matrix(df[, -1, drop = FALSE])     # remaining columns are per-sample counts
  rownames(mat) <- gene_ids
  # Raw-count loader: fractional values => normalized data => FAIL.
  mat <- .coerce_integer_matrix(mat, "matrix", allow_round = FALSE)
  list(counts = mat, feature_meta = NULL)
}

#' Load Subread featureCounts output.
#' Layout: optional '#' comment line, then Geneid, Chr, Start, End, Strand,
#' Length, then one column per sample. We drop the 6 annotation columns but keep
#' Length aside (feature_meta) for optional downstream TPM.
load_featurecounts <- function(path) {
  # skip="Geneid" tells fread to start at the header row, ignoring the comment.
  df <- .fread(path, skip = "Geneid")
  annot_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  gene_ids <- as.character(df$Geneid)
  # Preserve gene length: needed only if a user later wants TPM; not used for DE.
  feature_meta <- data.frame(gene_id = gene_ids, length = df$Length,
                             stringsAsFactors = FALSE)
  sample_cols <- setdiff(colnames(df), annot_cols)
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(mat) <- gene_ids
  # featureCounts column names are often full BAM paths; reduce to basenames.
  colnames(mat) <- .clean_sample_names(colnames(mat))
  mat <- .coerce_integer_matrix(mat, "featurecounts", allow_round = FALSE)
  list(counts = mat, feature_meta = feature_meta)
}

# Internal: turn BAM/file paths into tidy sample IDs (strip dir + known suffixes).
.clean_sample_names <- function(x) {
  x <- basename(x)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("\\.(txt|tsv|counts|tab)$", "", x, ignore.case = TRUE)
  x <- sub("ReadsPerGene\\.out$", "", x)
  x <- sub("[._]+$", "", x)   # strip any trailing separator left behind
  x
}

#' Load STAR per-sample *ReadsPerGene.out.tab files from a directory.
#' Each file: 4 columns (gene, unstranded, stranded-fwd, stranded-rev) with the
#' first 4 rows being N_unmapped / N_multimapping / N_noFeature / N_ambiguous.
#' @param strand_col Which count column to use (2=unstranded default, 3=fwd, 4=rev).
load_star <- function(path, strand_col = 2L) {
  files <- list.files(path, pattern = "ReadsPerGene\\.out\\.tab$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop_pipeline("No *ReadsPerGene.out.tab files found for STAR input.")
  # The count column index within the file is strand_col+... : col1 is gene id,
  # cols 2..4 are the three strandedness counts. strand_col=2 => unstranded.
  per_sample <- lapply(files, function(f) {
    d <- .fread(f, header = FALSE)
    # Drop the 4 summary rows STAR prepends before the per-gene counts.
    d <- d[-(1:4), , drop = FALSE]
    setNames(
      data.frame(gene = as.character(d[[1]]), count = d[[strand_col]],
                 stringsAsFactors = FALSE),
      c("gene", "count"))
  })
  genes <- per_sample[[1]]$gene
  # Sanity: all samples must share the same gene order (same STAR reference).
  for (i in seq_along(per_sample)) {
    if (!identical(per_sample[[i]]$gene, genes)) {
      stop_pipeline("STAR files disagree on gene order/identity; were they built against the same reference?")
    }
  }
  mat <- do.call(cbind, lapply(per_sample, `[[`, "count"))
  rownames(mat) <- genes
  colnames(mat) <- .clean_sample_names(files)
  mat <- .coerce_integer_matrix(mat, "star", allow_round = FALSE)
  log_info(sprintf("STAR: used strandedness column %d (2=unstranded, 3=fwd, 4=rev).", strand_col))
  list(counts = mat, feature_meta = NULL)
}

#' Load HTSeq-count per-sample 2-column files from a directory.
#' Each file: gene<TAB>count, with trailing summary rows __no_feature,
#' __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique.
load_htseq <- function(path) {
  files <- list.files(path, pattern = "\\.(txt|tsv|counts|tab)$",
                      recursive = TRUE, full.names = TRUE)
  # Keep only files that actually look like HTSeq output.
  files <- files[vapply(files, function(f) {
    any(grepl("^__", readLines(f, n = 200, warn = FALSE)))
  }, logical(1))]
  if (!length(files)) stop_pipeline("No HTSeq-count files (with __* summary rows) found.")
  per_sample <- lapply(files, function(f) {
    d <- .fread(f, header = FALSE)
    d <- d[!grepl("^__", d[[1]]), , drop = FALSE]   # drop __no_feature etc.
    setNames(data.frame(gene = as.character(d[[1]]), count = d[[2]],
                        stringsAsFactors = FALSE), c("gene", "count"))
  })
  genes <- per_sample[[1]]$gene
  for (i in seq_along(per_sample)) {
    if (!identical(per_sample[[i]]$gene, genes)) {
      stop_pipeline("HTSeq files disagree on gene identity/order.")
    }
  }
  mat <- do.call(cbind, lapply(per_sample, `[[`, "count"))
  rownames(mat) <- genes
  colnames(mat) <- .clean_sample_names(files)
  mat <- .coerce_integer_matrix(mat, "htseq", allow_round = FALSE)
  list(counts = mat, feature_meta = NULL)
}

#' Load Salmon/kallisto quantifications and summarise to gene level via tximport.
#' @param path Directory containing per-sample subdirs with quant.sf/abundance.tsv.
#' @param tx2gene Path to a 2-column transcript->gene map (required).
load_salmon <- function(path, tx2gene) {
  if (is.null(tx2gene) || !file.exists(tx2gene)) {
    stop_pipeline("Salmon/kallisto input requires --tx2gene <transcript,gene map file>.")
  }
  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop_pipeline("Package 'tximport' is required for salmon/kallisto input.")
  }
  quant_files <- list.files(path, pattern = "quant\\.sf$|abundance\\.tsv$",
                            recursive = TRUE, full.names = TRUE)
  if (!length(quant_files)) stop_pipeline("No quant.sf / abundance.tsv files found.")
  # Sample id = the name of the directory that contains each quant file.
  names(quant_files) <- basename(dirname(quant_files))
  type <- if (any(grepl("abundance\\.tsv$", quant_files))) "kallisto" else "salmon"
  t2g <- .fread(tx2gene, header = TRUE)
  # tximport with countsFromAbundance="no" returns library-size-scaled estimated
  # counts; these are fractional and MUST be rounded before DESeq2 (allow_round).
  txi <- tximport::tximport(quant_files, type = type, tx2gene = t2g[, 1:2],
                            countsFromAbundance = "no")
  mat <- .coerce_integer_matrix(txi$counts, "salmon", allow_round = TRUE)
  list(counts = mat, feature_meta = NULL)
}

#' Load counts + metadata from a serialized (Ranged)SummarizedExperiment (.rds).
#' Returns counts AND colData so airway-style objects need no separate metadata.
load_rangedSE <- function(path) {
  se <- readRDS(path)
  if (!methods::is(se, "SummarizedExperiment")) {
    stop_pipeline("--input-format rangedSE expects an .rds holding a (Ranged)SummarizedExperiment.")
  }
  mat <- as.matrix(SummarizedExperiment::assay(se))
  mat <- .coerce_integer_matrix(mat, "rangedSE", allow_round = FALSE)
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  list(counts = mat, feature_meta = NULL, coldata = coldata)
}

# ---------------------------------------------------------------------------
# Dispatch + sample reconciliation
# ---------------------------------------------------------------------------

#' Top-level loader: resolve the format (auto or forced) and call its loader.
#'
#' @param path --counts file or directory.
#' @param input_format One of the config values ("auto" resolves via detect_format).
#' @param tx2gene Optional path (required for salmon/kallisto).
#' @param star_strand STAR strandedness column (default 2 = unstranded).
#' @return list(counts, feature_meta[, coldata]).
load_counts <- function(path, input_format = "auto", tx2gene = NULL, star_strand = 2L) {
  fmt <- if (identical(input_format, "auto")) detect_format(path) else input_format
  res <- switch(fmt,
    matrix        = load_matrix(path),
    featurecounts = load_featurecounts(path),
    star          = load_star(path, strand_col = star_strand),
    htseq         = load_htseq(path),
    salmon        = load_salmon(path, tx2gene = tx2gene),
    rangedSE      = load_rangedSE(path),
    stop_pipeline(sprintf("Unknown --input-format '%s'.", fmt))
  )
  log_success(sprintf("Loaded %d genes x %d samples (format: %s).",
                     nrow(res$counts), ncol(res$counts), fmt))
  res$format <- fmt
  res
}

#' Reconcile count columns against the sample-sheet rows.
#'
#' @param counts Integer matrix (genes x samples).
#' @param metadata data.frame; sample IDs in column `sample_col`.
#' @param sample_col Name of the ID column in metadata.
#' @return Count matrix with columns reordered to match metadata row order.
#' FAILS loudly and prints the offending IDs (setdiff BOTH directions) on any
#' mismatch. WHY: silently intersecting samples is a classic source of wrong
#' results -- a dropped sample changes the model. Better to stop and show the IDs.
reconcile_samples <- function(counts, metadata, sample_col) {
  if (!sample_col %in% colnames(metadata)) {
    stop_pipeline(sprintf("Sample sheet has no '%s' column.", sample_col))
  }
  meta_ids  <- as.character(metadata[[sample_col]])
  count_ids <- colnames(counts)

  only_meta  <- setdiff(meta_ids, count_ids)
  only_count <- setdiff(count_ids, meta_ids)
  if (length(only_meta) || length(only_count)) {
    stop_pipeline(sprintf(paste0(
      "Sample IDs do not match between counts and metadata.\n",
      "  In metadata but NOT in counts (%d): %s\n",
      "  In counts but NOT in metadata (%d): %s\n",
      "Fix the sample sheet or the count column names so they align exactly."),
      length(only_meta),  paste(only_meta,  collapse = ", "),
      length(only_count), paste(only_count, collapse = ", ")))
  }
  # Reorder columns to the metadata order so counts and colData are row-aligned.
  counts[, meta_ids, drop = FALSE]
}
