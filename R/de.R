# =============================================================================
# R/de.R -- differential expression: results, LFC shrinkage, annotation, outputs
# -----------------------------------------------------------------------------
# One concern: turn a fitted DESeqDataSet into an annotated, thresholded,
# shrinkage-stabilized DE table plus a volcano plot.
#
# Why shrink log2 fold-changes: raw MLE fold-changes for low-count genes are
# wildly overdispersed and unreliable. apeglm shrinks them toward zero in
# proportion to their uncertainty, so ranking and visualization are dominated by
# genes with real, well-estimated effects rather than noisy low-count outliers.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

#' Extract results for the requested contrast.
#'
#' @param dds Fitted DESeqDataSet.
#' @param contrast [factor, numerator, denominator].
#' @param cfg Config (uses de$padj_cutoff for `alpha`, which tunes independent
#'   filtering to the significance level we actually use).
#' @return A DESeqResults object.
run_results <- function(dds, contrast, cfg) {
  res <- DESeq2::results(
    dds,
    contrast = c(contrast[[1]], contrast[[2]], contrast[[3]]),
    alpha = cfg$de$padj_cutoff)
  log_success(sprintf("Extracted results for %s: %s vs %s.",
                     contrast[[1]], contrast[[2]], contrast[[3]]))
  res
}

#' Apply log2 fold-change shrinkage.
#'
#' @param dds Fitted DESeqDataSet.
#' @param contrast [factor, numerator, denominator].
#' @param res The unshrunk DESeqResults (used as fallback / for ashr baseline).
#' @param method One of "apeglm","ashr","normal","none".
#' @return A DESeqResults object with shrunk log2FoldChange (or `res` if none).
#' apeglm shrinks a single coefficient and needs the coef NAME, not a contrast.
#' Because build_dds() set the denominator as the reference level, that coef is
#' "<factor>_<numerator>_vs_<denominator>". If it is not present (e.g. an
#' interaction design), we fall back to ashr, which accepts an explicit contrast.
shrink_lfc <- function(dds, contrast, res, method = "apeglm") {
  if (identical(method, "none")) {
    log_info("LFC shrinkage disabled (de$shrink: none).")
    return(res)
  }
  factor <- contrast[[1]]; num <- contrast[[2]]; denom <- contrast[[3]]
  coef_name <- sprintf("%s_%s_vs_%s", factor, num, denom)

  if (method == "apeglm") {
    if (coef_name %in% DESeq2::resultsNames(dds)) {
      log_info(sprintf("Shrinking LFCs with apeglm (coef '%s').", coef_name))
      return(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", res = res))
    }
    log_warn(sprintf(paste0("apeglm needs coef '%s' which is not in resultsNames(); ",
                            "falling back to ashr (supports arbitrary contrasts)."), coef_name))
    method <- "ashr"
  }
  if (method == "ashr") {
    if (!requireNamespace("ashr", quietly = TRUE)) {
      log_warn("ashr not installed; falling back to 'normal' shrinkage.")
      method <- "normal"
    } else {
      return(DESeq2::lfcShrink(dds, contrast = c(factor, num, denom), type = "ashr", res = res))
    }
  }
  # 'normal' shrinkage supports contrasts and ships with DESeq2.
  log_info("Shrinking LFCs with 'normal' estimator.")
  DESeq2::lfcShrink(dds, contrast = c(factor, num, denom), type = "normal", res = res)
}

#' Annotate result rows with gene SYMBOL and ENTREZ id via the organism OrgDb.
#'
#' @param res DESeqResults.
#' @param organism "human" or "mouse".
#' @return data.frame of results with added gene_id, symbol, entrez columns.
#' GUARDED: if the OrgDb is missing, we WARN and return the table with NA symbol
#' columns rather than crashing -- annotation is a convenience, not a correctness
#' requirement of the DE call itself.
annotate_results <- function(res, organism) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  # Ensembl gene IDs often carry a version suffix (ENSG00000123456.7); strip it
  # so the key matches the OrgDb ENSEMBL keytype.
  looks_ensembl <- mean(grepl("^ENS[A-Z]*G[0-9]+", df$gene_id)) > 0.5
  key <- if (looks_ensembl) sub("\\..*$", "", df$gene_id) else df$gene_id
  keytype <- if (looks_ensembl) "ENSEMBL" else "SYMBOL"

  orgdb_pkg <- if (organism == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
  if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
    log_warn(sprintf("%s not installed; skipping symbol/entrez annotation.", orgdb_pkg))
    df$symbol <- NA_character_; df$entrez <- NA_character_
    return(df)
  }
  orgdb <- getExportedValue(orgdb_pkg, orgdb_pkg)
  safe_map <- function(column) {
    tryCatch(
      suppressMessages(AnnotationDbi::mapIds(orgdb, keys = key, column = column,
                                             keytype = keytype, multiVals = "first")),
      error = function(e) { log_warn(sprintf("mapIds(%s) failed: %s", column, conditionMessage(e)));
                            setNames(rep(NA_character_, length(key)), key) })
  }
  # If already symbols, keep them; otherwise map from Ensembl.
  df$symbol <- if (keytype == "SYMBOL") df$gene_id else unname(safe_map("SYMBOL"))
  df$entrez <- unname(safe_map("ENTREZID"))
  df
}

#' Build the final, ordered, flagged DE table and write it to disk.
#'
#' @param df Annotated results data.frame.
#' @param cfg Config (padj_cutoff, lfc_cutoff).
#' @param outdir Output root; TSV written under <outdir>/tables/.
#' @return list(table = df with `significant` flag, path = tsv path, n_sig = int).
build_de_table <- function(df, cfg, outdir) {
  # A gene is "significant" when it clears BOTH the statistical (adjusted p) and
  # the effect-size (|log2FC|) thresholds. WHY both: statistical significance
  # alone can flag tiny, biologically irrelevant changes in high-power designs.
  df$significant <- !is.na(df$padj) &
    df$padj < cfg$de$padj_cutoff &
    abs(df$log2FoldChange) >= cfg$de$lfc_cutoff
  # Order by adjusted p-value (NA padj sink to the bottom).
  df <- df[order(df$padj, na.last = TRUE), ]

  # Reorder columns to a readable, documented layout.
  front <- intersect(c("gene_id", "symbol", "entrez", "baseMean",
                        "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                        "significant"), colnames(df))
  df <- df[, c(front, setdiff(colnames(df), front))]

  tables_dir <- file.path(outdir, "tables")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(tables_dir, "de_results.tsv")
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

  n_sig <- sum(df$significant, na.rm = TRUE)
  log_success(sprintf("DE table written to %s (%d significant genes at padj<%.3g, |LFC|>=%.2g).",
                     path, n_sig, cfg$de$padj_cutoff, cfg$de$lfc_cutoff))
  list(table = df, path = path, n_sig = n_sig)
}

#' Draw and save a volcano plot (EnhancedVolcano, with a ggplot fallback).
#'
#' @param df DE table (with log2FoldChange, padj, symbol).
#' @param cfg Config (thresholds define the guide lines).
#' @param outdir Output root; PNG under <outdir>/plots/.
#' @return Path to the PNG.
plot_volcano <- function(df, cfg, outdir) {
  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(plots_dir, "volcano.png")
  labs_vec <- ifelse(is.na(df$symbol) | df$symbol == "", df$gene_id, df$symbol)

  ok <- tryCatch({
    p <- EnhancedVolcano::EnhancedVolcano(
      df, lab = labs_vec, x = "log2FoldChange", y = "padj",
      pCutoff = cfg$de$padj_cutoff, FCcutoff = cfg$de$lfc_cutoff,
      title = "Differential expression", subtitle = NULL,
      ylab = bquote(~-Log[10] ~ italic(padj)))
    ggplot2::ggsave(path, p, width = 8, height = 7, dpi = 150)
    TRUE
  }, error = function(e) { log_warn(sprintf("EnhancedVolcano failed (%s); using ggplot fallback.",
                                            conditionMessage(e))); FALSE })

  if (!ok) {
    d <- df[!is.na(df$padj), ]
    d$neglog10 <- -log10(d$padj)
    d$dir <- ifelse(d$significant & d$log2FoldChange > 0, "up",
                    ifelse(d$significant & d$log2FoldChange < 0, "down", "ns"))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["log2FoldChange"]],
                                         y = .data[["neglog10"]],
                                         color = .data[["dir"]])) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_vline(xintercept = c(-1, 1) * cfg$de$lfc_cutoff, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(cfg$de$padj_cutoff), linetype = "dashed") +
      ggplot2::scale_color_manual(values = c(up = "firebrick", down = "steelblue", ns = "grey70")) +
      ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p", title = "Differential expression") +
      ggplot2::theme_bw(base_size = 13)
    ggplot2::ggsave(path, p, width = 8, height = 7, dpi = 150)
  }
  log_success(sprintf("Volcano plot saved to %s", path))
  path
}
