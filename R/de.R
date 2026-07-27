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

# Coefficient name implied by a [factor, numerator, denominator] contrast, given
# the denominator was set as the reference level in build_dds().
.coef_from_contrast <- function(ct) sprintf("%s_%s_vs_%s", ct[[1]], ct[[2]], ct[[3]])

#' Resolve the config into an ordered list of result specifications.
#'
#' @param cfg Config list.
#' @return A list of specs; each is list(name, contrast=[f,n,d] or NULL, coef=str or NULL).
#' The primary `de$contrast` is always first (keeps single-contrast behaviour and
#' file names unchanged). `de$contrasts` (optional) adds more results, each either
#' a plain [factor, numerator, denominator] vector OR a named entry with a `coef`
#' (a resultsNames() coefficient -- how you request an INTERACTION effect) or a
#' `contrast`. This is what lets one run emit several comparisons.
resolve_result_specs <- function(cfg) {
  specs <- list(list(name = .coef_from_contrast(cfg$contrast),
                     contrast = cfg$contrast, coef = NULL))
  extra <- cfg$de$contrasts
  if (!is.null(extra)) {
    for (item in extra) {
      if (is.null(names(item))) {
        # Plain [factor, numerator, denominator] vector.
        specs[[length(specs) + 1]] <- list(name = .coef_from_contrast(item),
                                            contrast = item, coef = NULL)
      } else if (!is.null(item$coef)) {
        nm <- if (!is.null(item$name)) item$name else item$coef
        specs[[length(specs) + 1]] <- list(name = nm, contrast = NULL, coef = item$coef)
      } else {
        ct <- item$contrast
        nm <- if (!is.null(item$name)) item$name else .coef_from_contrast(ct)
        specs[[length(specs) + 1]] <- list(name = nm, contrast = ct, coef = NULL)
      }
    }
  }
  # De-duplicate by name, keeping first occurrence (primary wins).
  specs[!duplicated(vapply(specs, `[[`, character(1), "name"))]
}

#' Extract results for one result spec (a contrast OR a named coefficient).
#'
#' @param dds Fitted DESeqDataSet.
#' @param spec A spec from resolve_result_specs().
#' @param cfg Config (uses de$padj_cutoff for `alpha`).
#' @return A DESeqResults object.
#' NOTE on LRT: when `dds` was fit with test="LRT", results() returns the LRT
#' p-values regardless of the contrast/name; the contrast/coef only selects which
#' log2FoldChange to report alongside those p-values.
run_results <- function(dds, spec, cfg) {
  if (!is.null(spec$coef)) {
    if (!spec$coef %in% DESeq2::resultsNames(dds)) {
      stop_pipeline(sprintf("Requested coefficient '%s' is not in resultsNames(): %s",
                            spec$coef, paste(DESeq2::resultsNames(dds), collapse = ", ")))
    }
    res <- DESeq2::results(dds, name = spec$coef, alpha = cfg$de$padj_cutoff)
    log_success(sprintf("Extracted results for coefficient '%s'.", spec$coef))
  } else {
    ct <- spec$contrast
    res <- DESeq2::results(dds, contrast = c(ct[[1]], ct[[2]], ct[[3]]),
                           alpha = cfg$de$padj_cutoff)
    log_success(sprintf("Extracted results for %s: %s vs %s.", ct[[1]], ct[[2]], ct[[3]]))
  }
  res
}

#' Apply log2 fold-change shrinkage for one result spec.
#'
#' @param dds Fitted DESeqDataSet.
#' @param spec A spec from resolve_result_specs().
#' @param res The unshrunk DESeqResults.
#' @param method One of "apeglm","ashr","normal","none".
#' @return A DESeqResults object with shrunk log2FoldChange (or `res` if none).
#' apeglm shrinks a single coefficient by NAME, so it is the natural choice for
#' both contrasts (coef derived from the reference level) and interaction coefs.
#' ashr/normal need an explicit contrast; a coef-only spec (interaction) therefore
#' cannot use them and is returned unshrunk with a warning.
shrink_lfc <- function(dds, spec, res, method = "apeglm") {
  if (identical(method, "none")) {
    log_info("LFC shrinkage disabled (de$shrink: none).")
    return(res)
  }
  coef_name <- if (!is.null(spec$coef)) spec$coef else .coef_from_contrast(spec$contrast)

  if (method == "apeglm") {
    if (coef_name %in% DESeq2::resultsNames(dds)) {
      log_info(sprintf("Shrinking LFCs with apeglm (coef '%s').", coef_name))
      return(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", res = res))
    }
    log_warn(sprintf("apeglm needs coef '%s' (not in resultsNames()); trying ashr.", coef_name))
    method <- "ashr"
  }
  # ashr/normal require a contrast vector; interaction (coef-only) specs can't use them.
  if (is.null(spec$contrast)) {
    log_warn(sprintf("Cannot apply '%s' shrinkage to coefficient-only spec '%s'; returning unshrunk LFCs.",
                     method, spec$name))
    return(res)
  }
  ct <- spec$contrast
  if (method == "ashr") {
    if (!requireNamespace("ashr", quietly = TRUE)) {
      log_warn("ashr not installed; falling back to 'normal' shrinkage.")
      method <- "normal"
    } else {
      return(DESeq2::lfcShrink(dds, contrast = c(ct[[1]], ct[[2]], ct[[3]]), type = "ashr", res = res))
    }
  }
  log_info("Shrinking LFCs with 'normal' estimator.")
  DESeq2::lfcShrink(dds, contrast = c(ct[[1]], ct[[2]], ct[[3]]), type = "normal", res = res)
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
#' @param suffix Optional filename suffix (e.g. "_myContrast") so multiple
#'   contrasts write to distinct files; "" keeps the primary de_results.tsv name.
#' @return list(table = df with `significant` flag, path = tsv path, n_sig = int).
build_de_table <- function(df, cfg, outdir, suffix = "") {
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
  path <- file.path(tables_dir, sprintf("de_results%s.tsv", suffix))
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
#' @param suffix Optional filename suffix for multi-contrast runs ("" = volcano.png).
#' @param title Plot title (defaults to "Differential expression").
#' @return Path to the PNG.
plot_volcano <- function(df, cfg, outdir, suffix = "", title = "Differential expression") {
  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(plots_dir, sprintf("volcano%s.png", suffix))
  labs_vec <- ifelse(is.na(df$symbol) | df$symbol == "", df$gene_id, df$symbol)

  ok <- tryCatch({
    p <- EnhancedVolcano::EnhancedVolcano(
      df, lab = labs_vec, x = "log2FoldChange", y = "padj",
      pCutoff = cfg$de$padj_cutoff, FCcutoff = cfg$de$lfc_cutoff,
      title = title, subtitle = NULL,
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
