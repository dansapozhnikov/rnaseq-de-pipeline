# =============================================================================
# R/pipeline.R -- the end-to-end orchestration of a single analysis
# -----------------------------------------------------------------------------
# One concern: given a raw count matrix + metadata + resolved config, run every
# stage in order and produce all artifacts. Shared by real runs (run_pipeline.R)
# and the airway self-test (tests/validate_airway.R) so both exercise the SAME
# code path -- the validation genuinely tests what users run.
#
# Stage order mirrors docs/pipeline.svg:
#   reconcile -> QC (pre-model) -> GATE -> filter -> DESeq -> VST
#   -> batch view -> PCA + variance-vs-metadata -> QC (post-model)
#   -> DE (+shrink+annotate+volcano) -> enrichment -> report
# =============================================================================

#' Run the full pipeline on an in-memory count matrix + metadata.
#'
#' @param counts Integer matrix (genes x samples).
#' @param metadata data.frame; sample IDs in `sample_col`.
#' @param cfg Resolved config list (from load_config()).
#' @param outdir Output root.
#' @param run_timestamp Timestamp string for provenance/report.
#' @param sample_col Metadata column holding sample IDs (default "sample").
#' @param tx2gene Unused here (kept for signature symmetry with loaders).
#' @param force If TRUE, downgrade QC FAIL -> WARN and continue.
#' @param resume If TRUE (default), reuse cached intermediates from a prior run
#'   when inputs + config + code are unchanged; only recompute what changed.
#' @return list(qc, de, dds, report) -- used by the airway validation to assert.
run_pipeline_core <- function(counts, metadata, cfg, outdir, run_timestamp,
                              sample_col = "sample", tx2gene = NULL, force = FALSE,
                              pipeline_version = "dev", resume = TRUE) {
  # Determinism: fix the seed before any stochastic step (shrinkage, plotting).
  set.seed(cfg$seed)
  tested <- cfg$contrast[[1]]
  t0 <- Sys.time()   # wall-clock start, for the run-metrics panel

  # Headless-safe plotting on Linux/HPC. WHY: on Linux compute nodes without an
  # X11 display the default bitmap device can fail when saving PNGs; the Cairo
  # device renders off-screen. We scope this to Linux ONLY: macOS already renders
  # headless via its native 'quartz' device (and its cairo build may depend on an
  # absent X11), and Windows uses its own device -- forcing cairo there would hurt.
  if (Sys.info()[["sysname"]] == "Linux" && isTRUE(capabilities("cairo"))) {
    options(bitmapType = "cairo")
  }

  # --- Align counts <-> metadata (identity-checked; FAIL on mismatch) --------
  counts <- reconcile_samples(counts, metadata, sample_col)
  rownames(metadata) <- as.character(metadata[[sample_col]])
  metadata <- metadata[colnames(counts), , drop = FALSE]

  # --- Resumability: build the checkpoint cache ------------------------------
  # The base key captures everything run-wide that must invalidate ALL cached
  # stages when it changes: the count data, the sample sheet, the design/contrast,
  # the seed, the low-count-filter thresholds, and a fingerprint of the compute
  # code itself. Stage-specific inputs are added per stage below.
  cache <- new_cache(
    outdir, enable = isTRUE(resume), version = pipeline_version,
    base_deps = list(counts, metadata, cfg$design, cfg$contrast, cfg$seed,
                     cfg$qc$min_count, cfg$qc$min_replicates_per_group,
                     cache_code_fingerprint()))
  if (isTRUE(resume)) {
    log_info(sprintf("Resume enabled: checkpoint cache at %s (finished stages are reused).",
                     cache$dir))
  } else {
    log_info("Resume disabled (--no-resume): every stage computed fresh.")
  }

  # ===========================================================================
  # QC -- pre-model checks (the ones that can FAIL and stop the run)
  # ===========================================================================
  log_section("Quality control (pre-model)")
  qc <- qc_init()
  qc_check_sample_alignment(qc, counts, metadata)
  qc_check_raw_counts(qc, counts)
  qc_check_library_size(qc, counts, cfg)
  qc_check_genes_detected(qc, counts, cfg)
  qc_check_replicates(qc, metadata, tested, cfg)
  qc_check_design_identifiable(qc, metadata, cfg$design, tested)

  # Persist QC + render a QC-only report BEFORE gating, so a FAIL still leaves
  # the user with a readable red/amber/green report and the TSV (per spec).
  qc_write(qc, outdir)
  if (any(qc$results$status == "FAIL") && !isTRUE(force)) {
    qc_summary(qc)
    render_report(cfg, outdir, run_timestamp, character(0), forced = FALSE)
  }
  qc_gate(qc, force = force)   # stops with exit 1 on FAIL unless --force

  # ===========================================================================
  # Normalize: build, filter, fit
  # ===========================================================================
  log_section("Normalization & model fit")
  # The DESeq fit (size factors + dispersions + GLM) is the single most expensive
  # stage, so it is the primary checkpoint. We cache the fitted object together
  # with the pre/post-filter gene counts (needed by QC + metrics on resume).
  test_type <- cfg$de$test %||% "Wald"
  fit <- with_cache(
    cache, "deseq_fit", deps = list(test_type, cfg$de$reduced),
    compute = function() {
      dds0 <- build_dds(counts, metadata, cfg$design, cfg$contrast)
      flt0 <- filter_low_counts(dds0, cfg)
      list(dds = run_deseq(flt0$dds, test = test_type, reduced = cfg$de$reduced),
           n_before = flt0$n_before, n_after = flt0$n_after)
    })
  dds <- fit$dds
  flt <- list(n_before = fit$n_before, n_after = fit$n_after)
  dds_key <- key_of(fit)   # thread as `parent` so downstream cascades on refit
  # QC checks always run (cheap) so the report reflects this run regardless of cache.
  qc_check_low_count_filter(qc, flt$n_before, flt$n_after, cfg)
  qc_check_dispersion(qc, dds)

  # ===========================================================================
  # Explore: VST -> batch view -> PCA -> variance-vs-metadata
  # ===========================================================================
  log_section("Exploratory analysis")
  # VST is moderately expensive on large data; cache it, chained to the fit key.
  vsd <- with_cache(cache, "vst", deps = list(blind = TRUE), parent = dds_key,
                    compute = function() get_vst(dds, blind = TRUE))

  # Number of top-variance genes for PCA is config-driven (no magic number);
  # default to the DESeq2 convention of 500 if the field is absent.
  pca_ntop <- cfg$explore$pca_ntop %||% 500L
  pca <- compute_pca(vsd, ntop = pca_ntop)
  shape_by <- {
    bv <- identify_batch_vars(cfg$design, tested)
    bv <- bv[bv %in% colnames(metadata)]
    if (length(bv)) bv[1] else NULL
  }
  plot_pca(pca, as.data.frame(metadata), color_by = tested, shape_by = shape_by, outdir = outdir)
  variance_vs_metadata(pca, as.data.frame(metadata), outdir = outdir)

  # Batch-corrected VIEW (visualization only; DE still uses the raw-count design).
  batch_vars <- identify_batch_vars(cfg$design, tested)
  batch_vars <- batch_vars[batch_vars %in% colnames(metadata)]
  if (length(batch_vars)) {
    vsd_bc <- remove_batch_for_viz(vsd, as.data.frame(metadata), batch_vars, tested)
    pca_bc <- compute_pca(vsd_bc, ntop = pca_ntop)
    # Persist the corrected scores too, so the report can render an INTERACTIVE
    # batch-corrected PCA (matching the main one) in the Batch effects section.
    plot_pca(pca_bc, as.data.frame(metadata), color_by = tested, shape_by = shape_by,
             outdir = outdir, filename = "pca_batch_corrected.png", save_scores = TRUE,
             scores_file = "pca_scores_bc.tsv", pv_file = "pca_percentvar_bc.txt")
  }

  # Post-model outlier check (WARN only; never auto-removes samples).
  qc_check_cooks_outliers(qc, dds, cfg)

  # ===========================================================================
  # Differential expression -- one or more result specs (contrasts / coefficients)
  # ===========================================================================
  log_section("Differential expression")
  specs <- resolve_result_specs(cfg)
  log_info(sprintf("Extracting %d result(s): %s", length(specs),
                   paste(vapply(specs, `[[`, character(1), "name"), collapse = ", ")))
  de_all <- list()
  primary_de_key <- NULL
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    # Primary spec keeps the unsuffixed de_results.tsv / volcano.png names so the
    # report + airway validation are unchanged; extras get a sanitised suffix.
    suffix <- if (i == 1) "" else paste0("_", gsub("[^A-Za-z0-9]+", "_", spec$name))
    # Cache the shrunk + annotated results (the apeglm shrinkage and OrgDb lookup
    # are the cost). Thresholding (significant flag), the TSV, and the volcano are
    # cheap and always regenerated -- so changing padj/lfc cutoffs does NOT force a
    # re-shrink (padj itself depends on alpha, so that stays in the key).
    stage <- paste0("de_", gsub("[^A-Za-z0-9]+", "_", spec$name))
    de_df <- with_cache(
      cache, stage,
      deps = list(spec, cfg$de$shrink, cfg$de$padj_cutoff, cfg$organism),
      parent = dds_key,
      compute = function() {
        res <- run_results(dds, spec, cfg)
        res <- shrink_lfc(dds, spec, res, cfg$de$shrink)
        annotate_results(res, cfg$organism)
      })
    if (i == 1) primary_de_key <- key_of(de_df)
    detab_i <- build_de_table(de_df, cfg, outdir, suffix = suffix)
    plot_volcano(detab_i$table, cfg, outdir, suffix = suffix,
                 title = sprintf("Differential expression: %s", spec$name))
    detab_i$name <- spec$name
    de_all[[spec$name]] <- detab_i
  }
  detab <- de_all[[1]]   # primary contrast (back-compatible return + report featured view)

  # ===========================================================================
  # Functional enrichment (guarded) -- run on the PRIMARY contrast
  # ===========================================================================
  log_section("Functional enrichment")
  # ORA/GSEA (GSEA especially) are the cost here; cache them chained to the
  # primary DE key. run_enrichment also WRITES tsv/plot side effects, so on a
  # cache hit we verify those declared files still exist -- if any are missing
  # (e.g. outputs cleared but cache kept), we recompute rather than point at gaps.
  enr <- with_cache(
    cache, "enrichment",
    deps = list(cfg$enrichment, cfg$organism, cfg$de$padj_cutoff, cfg$de$lfc_cutoff),
    parent = primary_de_key,
    verify = function(e) {
      paths <- Filter(Negate(is.null), list(e$ora_path, e$gsea_path, e$dotplot))
      all(vapply(paths, file.exists, logical(1)))
    },
    compute = function() run_enrichment(detab$table, cfg, cfg$organism, outdir))

  # ===========================================================================
  # Run metrics & environment provenance
  # ===========================================================================
  # WHY: a report should be self-describing -- anyone reading it later can see how
  # long it took, which pipeline version + package versions produced it, and the
  # headline counts, without digging through the log. Written as a 2-column TSV
  # the report renders verbatim.
  sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NA_real_)
  sig <- detab$table$significant %in% c(TRUE, "TRUE")
  pkgver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "NA")
  metrics <- data.frame(
    metric = c(
      "Pipeline version", "Run timestamp", "Runtime",
      "Samples", "Genes (input)", "Genes (after filter)", "Genes filtered",
      "Test", "Contrasts", "Significant DEGs", "DEGs up / down", "Size factors (min-max)",
      "Organism", "Design", "Contrast",
      "R version", "Bioconductor", "DESeq2", "apeglm", "clusterProfiler",
      "Platform", "Host"),
    value = c(
      pipeline_version, run_timestamp,
      sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))),
      ncol(counts), flt$n_before, flt$n_after,
      sprintf("%.1f%%", 100 * (flt$n_before - flt$n_after) / flt$n_before),
      cfg$de$test %||% "Wald", length(specs),
      detab$n_sig, sprintf("%d / %d", sum(sig & detab$table$log2FoldChange > 0, na.rm = TRUE),
                                       sum(sig & detab$table$log2FoldChange < 0, na.rm = TRUE)),
      if (all(is.na(sf))) "NA" else sprintf("%.2f - %.2f", min(sf), max(sf)),
      cfg$organism, cfg$design,
      sprintf("%s: %s vs %s", cfg$contrast[[1]], cfg$contrast[[2]], cfg$contrast[[3]]),
      R.version.string, tryCatch(as.character(BiocManager::version()), error = function(e) "NA"),
      pkgver("DESeq2"), pkgver("apeglm"), pkgver("clusterProfiler"),
      R.version$platform, Sys.info()[["nodename"]]),
    stringsAsFactors = FALSE)
  utils::write.table(metrics, file.path(outdir, "qc", "run_metrics.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  log_info(sprintf("Run metrics written (runtime %s).", metrics$value[3]))

  # ===========================================================================
  # Final QC summary + report
  # ===========================================================================
  qc_write(qc, outdir)
  qc_summary(qc)
  log_section("Report")
  report_path <- render_report(cfg, outdir, run_timestamp, batch_vars,
                               forced = isTRUE(qc$forced))

  list(qc = qc, de = detab, de_all = de_all, specs = specs,
       dds = dds, enrichment = enr, report = report_path)
}
