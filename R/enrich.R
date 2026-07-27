# =============================================================================
# R/enrich.R -- functional enrichment (ORA + GSEA), fully guarded
# -----------------------------------------------------------------------------
# One concern: given the DE table, ask which GO biological processes are
# over-represented among the significant genes (ORA) and which are coordinately
# shifted across the whole ranked list (GSEA).
#
# GUARDING PRINCIPLE: enrichment is the LAST, most fragile step (needs an OrgDb,
# enough mapped genes, non-degenerate statistics). Every failure mode here must
# WARN and return gracefully -- never abort a run that already produced valid DE
# results. A missing OrgDb, zero significant genes, or a solver hiccup all yield
# a logged warning and NULL, not a crash.
# =============================================================================

#' Run GO enrichment (ORA + GSEA) on a DE table.
#'
#' @param de_table DE table from build_de_table() (needs entrez, log2FoldChange,
#'   stat, padj, significant).
#' @param cfg Config (enrichment$enable, $ontology, $padj_cutoff).
#' @param organism "human" or "mouse".
#' @param outdir Output root; tables under <outdir>/tables/, plot under /plots/.
#' @return list(ora, gsea, ora_path, gsea_path, dotplot) -- any element may be NULL.
run_enrichment <- function(de_table, cfg, organism, outdir) {
  null_result <- list(ora = NULL, gsea = NULL, ora_path = NULL,
                      gsea_path = NULL, dotplot = NULL)

  if (!isTRUE(cfg$enrichment$enable)) {
    log_info("Enrichment disabled in config; skipping.")
    return(null_result)
  }
  # --- Dependency guards -----------------------------------------------------
  orgdb_pkg <- if (organism == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    log_warn("clusterProfiler not installed; skipping enrichment."); return(null_result)
  }
  if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
    log_warn(sprintf("%s not installed; skipping enrichment.", orgdb_pkg)); return(null_result)
  }
  orgdb <- getExportedValue(orgdb_pkg, orgdb_pkg)
  ont <- cfg$enrichment$ontology
  pcut <- cfg$enrichment$padj_cutoff %||% cfg$de$padj_cutoff

  tables_dir <- file.path(outdir, "tables"); dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  plots_dir  <- file.path(outdir, "plots");  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

  # Universe = all genes with a mapped Entrez id (the tested background).
  universe <- unique(stats::na.omit(as.character(de_table$entrez)))
  if (length(universe) < 10) {
    log_warn("Fewer than 10 genes have Entrez IDs; enrichment is not meaningful. Skipping.")
    return(null_result)
  }

  # ---------------------------------------------------------------------------
  # (1) ORA -- over-representation of GO terms among significant genes.
  # ---------------------------------------------------------------------------
  sig_entrez <- unique(stats::na.omit(as.character(de_table$entrez[de_table$significant])))
  ora <- NULL; ora_path <- NULL
  if (length(sig_entrez) < 5) {
    log_warn(sprintf("Only %d significant genes with Entrez IDs; skipping ORA.", length(sig_entrez)))
  } else {
    ora <- tryCatch(
      clusterProfiler::enrichGO(
        gene = sig_entrez, universe = universe, OrgDb = orgdb,
        keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH",
        pvalueCutoff = pcut, qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) { log_warn(sprintf("enrichGO failed: %s", conditionMessage(e))); NULL })
    if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
      ora_path <- file.path(tables_dir, "enrichment_ora.tsv")
      utils::write.table(as.data.frame(ora), ora_path, sep = "\t", quote = FALSE, row.names = FALSE)
      log_success(sprintf("ORA: %d enriched GO:%s terms -> %s",
                         nrow(as.data.frame(ora)), ont, ora_path))
    } else {
      log_info("ORA produced no significantly enriched terms.")
    }
  }

  # ---------------------------------------------------------------------------
  # (2) GSEA -- coordinated shifts across the full ranked gene list.
  # Rank by the Wald statistic (sign = direction, magnitude = confidence).
  # ---------------------------------------------------------------------------
  gsea <- NULL; gsea_path <- NULL
  rank_df <- de_table[!is.na(de_table$entrez), ]
  rank_df <- rank_df[!duplicated(rank_df$entrez), ]
  # Ranking metric: prefer the Wald statistic, but apeglm/ashr shrinkage DROPS the
  # `stat` column, so fall back to signed -log10(p) = sign(LFC) * -log10(pvalue),
  # which encodes both direction and confidence and is a standard GSEA ranking.
  metric <- rank_df$stat
  if (is.null(metric) || all(is.na(metric))) {
    metric <- sign(rank_df$log2FoldChange) *
      -log10(pmax(rank_df$pvalue, .Machine$double.xmin))
  }
  rank_df$metric <- metric
  rank_df <- rank_df[is.finite(rank_df$metric), ]
  if (nrow(rank_df) >= 20) {
    ranks <- sort(setNames(rank_df$metric, rank_df$entrez), decreasing = TRUE)
    gsea <- tryCatch(
      suppressWarnings(clusterProfiler::gseGO(
        geneList = ranks, OrgDb = orgdb, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = pcut, verbose = FALSE)),
      error = function(e) { log_warn(sprintf("gseGO failed: %s", conditionMessage(e))); NULL })
    if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
      gsea_path <- file.path(tables_dir, "enrichment_gsea.tsv")
      utils::write.table(as.data.frame(gsea), gsea_path, sep = "\t", quote = FALSE, row.names = FALSE)
      log_success(sprintf("GSEA: %d enriched GO:%s gene sets -> %s",
                         nrow(as.data.frame(gsea)), ont, gsea_path))
    } else {
      log_info("GSEA produced no significantly enriched gene sets.")
    }
  } else {
    log_warn(sprintf("Only %d ranked genes; skipping GSEA (needs >= 20).", nrow(rank_df)))
  }

  # ---------------------------------------------------------------------------
  # (3) Dotplot of the ORA result (most interpretable single figure).
  # ---------------------------------------------------------------------------
  dotplot_path <- NULL
  if (!is.null(ora) && nrow(as.data.frame(ora)) > 0 &&
      requireNamespace("enrichplot", quietly = TRUE)) {
    dotplot_path <- file.path(plots_dir, "enrichment_dotplot.png")
    ok <- tryCatch({
      p <- enrichplot::dotplot(ora, showCategory = 15) +
        ggplot2::ggtitle(sprintf("Over-represented GO:%s terms", ont))
      ggplot2::ggsave(dotplot_path, p, width = 8, height = 7, dpi = 150); TRUE
    }, error = function(e) { log_warn(sprintf("dotplot failed: %s", conditionMessage(e))); FALSE })
    if (!ok) dotplot_path <- NULL
  }

  list(ora = ora, gsea = gsea, ora_path = ora_path,
       gsea_path = gsea_path, dotplot = dotplot_path)
}

# NULL-coalescing helper (used for optional config values).
`%||%` <- function(a, b) if (is.null(a)) b else a
