# =============================================================================
# R/pipeline.R — the end-to-end orchestration of a single analysis
# -----------------------------------------------------------------------------
# One concern: given a raw count matrix + metadata + resolved config, run every
# stage in order and produce all artifacts. Shared by real runs (run_pipeline.R)
# and the airway self-test (tests/validate_airway.R) so both exercise the SAME
# code path — the validation genuinely tests what users run.
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
#' @return list(qc, de, dds, report) — used by the airway validation to assert.
run_pipeline_core <- function(counts, metadata, cfg, outdir, run_timestamp,
                              sample_col = "sample", tx2gene = NULL, force = FALSE) {
  # Determinism: fix the seed before any stochastic step (shrinkage, plotting).
  set.seed(cfg$seed)
  tested <- cfg$contrast[[1]]

  # --- Align counts <-> metadata (identity-checked; FAIL on mismatch) --------
  counts <- reconcile_samples(counts, metadata, sample_col)
  rownames(metadata) <- as.character(metadata[[sample_col]])
  metadata <- metadata[colnames(counts), , drop = FALSE]

  # ===========================================================================
  # QC — pre-model checks (the ones that can FAIL and stop the run)
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
  dds <- build_dds(counts, metadata, cfg$design, cfg$contrast)
  flt <- filter_low_counts(dds, cfg)
  qc_check_low_count_filter(qc, flt$n_before, flt$n_after, cfg)
  dds <- run_deseq(flt$dds)
  qc_check_dispersion(qc, dds)

  # ===========================================================================
  # Explore: VST -> batch view -> PCA -> variance-vs-metadata
  # ===========================================================================
  log_section("Exploratory analysis")
  vsd <- get_vst(dds, blind = TRUE)

  pca <- compute_pca(vsd)
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
    pca_bc <- compute_pca(vsd_bc)
    plot_pca(pca_bc, as.data.frame(metadata), color_by = tested, shape_by = shape_by,
             outdir = outdir, filename = "pca_batch_corrected.png", save_scores = FALSE)
  }

  # Post-model outlier check (WARN only; never auto-removes samples).
  qc_check_cooks_outliers(qc, dds, cfg)

  # ===========================================================================
  # Differential expression
  # ===========================================================================
  log_section("Differential expression")
  res <- run_results(dds, cfg$contrast, cfg)
  res <- shrink_lfc(dds, cfg$contrast, res, cfg$de$shrink)
  de_df <- annotate_results(res, cfg$organism)
  detab <- build_de_table(de_df, cfg, outdir)
  plot_volcano(detab$table, cfg, outdir)

  # ===========================================================================
  # Functional enrichment (guarded)
  # ===========================================================================
  log_section("Functional enrichment")
  enr <- run_enrichment(detab$table, cfg, cfg$organism, outdir)

  # ===========================================================================
  # Final QC summary + report
  # ===========================================================================
  qc_write(qc, outdir)
  qc_summary(qc)
  log_section("Report")
  report_path <- render_report(cfg, outdir, run_timestamp, batch_vars,
                               forced = isTRUE(qc$forced))

  list(qc = qc, de = detab, dds = dds, enrichment = enr, report = report_path)
}
