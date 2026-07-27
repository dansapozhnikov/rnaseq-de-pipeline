#!/usr/bin/env Rscript
# Generate tiny, deterministic fixtures for the io_detect loaders. Run from the
# repo root: Rscript tests/fixtures/make_fixtures.R
# WHY committed generator + committed outputs: the outputs are what tests read,
# but the generator documents exactly how each fixture was shaped.
suppressPackageStartupMessages(library(SummarizedExperiment))

fx <- "tests/fixtures"
dir.create(fx, showWarnings = FALSE, recursive = TRUE)

genes   <- c("G1", "G2", "G3", "G4", "G5")
samples <- c("s1", "s2", "s3", "s4")
# Fixed integer count matrix (genes x samples).
m <- matrix(c(100, 0, 50, 200, 10,
              120, 2, 60, 180, 12,
              5,  300, 40, 20, 500,
              8,  280, 45, 25, 520),
            nrow = 5, dimnames = list(genes, samples))

# 1) Plain matrix (CSV).
df <- data.frame(gene = rownames(m), m, check.names = FALSE)
write.csv(df, file.path(fx, "matrix.csv"), row.names = FALSE, quote = FALSE)

# 2) featureCounts (leading '#' comment; annotation cols; Length kept).
fc_path <- file.path(fx, "featurecounts.txt")
cat("# Program:featureCounts v2.0.1; Command:...\n", file = fc_path)
hdr <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", samples)
cat(paste(hdr, collapse = "\t"), "\n", file = fc_path, append = TRUE)
for (i in seq_along(genes)) {
  row <- c(genes[i], "chr1", 100 * i, 100 * i + 50, "+", 1000 + i, m[i, ])
  cat(paste(row, collapse = "\t"), "\n", file = fc_path, append = TRUE)
}

# 3) STAR per-sample *ReadsPerGene.out.tab (4 summary rows + 4 count cols).
star_dir <- file.path(fx, "star"); dir.create(star_dir, showWarnings = FALSE)
for (s in samples) {
  f <- file.path(star_dir, paste0(s, "_ReadsPerGene.out.tab"))
  summ <- rbind(c("N_unmapped", 10, 10, 10),
                c("N_multimapping", 20, 20, 20),
                c("N_noFeature", 30, 30, 30),
                c("N_ambiguous", 5, 5, 5))
  body <- cbind(genes, m[, s], m[, s], m[, s])  # unstranded=col2 == our matrix
  write.table(rbind(summ, body), f, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# 4) HTSeq per-sample 2-col files with trailing __* rows.
htseq_dir <- file.path(fx, "htseq"); dir.create(htseq_dir, showWarnings = FALSE)
for (s in samples) {
  f <- file.path(htseq_dir, paste0(s, ".txt"))
  body <- cbind(genes, m[, s])
  tail <- rbind(c("__no_feature", 30), c("__ambiguous", 5),
                c("__too_low_aQual", 0), c("__not_aligned", 0),
                c("__alignment_not_unique", 0))
  write.table(rbind(body, tail), f, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# 5) Salmon: per-sample quant.sf at transcript level + tx2gene map.
#    Two transcripts per gene; gene counts must sum back to the matrix.
salmon_dir <- file.path(fx, "salmon"); dir.create(salmon_dir, showWarnings = FALSE)
tx  <- paste0(rep(genes, each = 2), c(".t1", ".t2"))
t2g <- data.frame(tx = tx, gene = rep(genes, each = 2))
write.csv(t2g, file.path(fx, "tx2gene.csv"), row.names = FALSE, quote = FALSE)
for (s in samples) {
  sdir <- file.path(salmon_dir, s); dir.create(sdir, showWarnings = FALSE)
  # Split each gene's count across its two transcripts.
  reads <- as.numeric(m[, s])
  num   <- as.vector(rbind(ceiling(reads / 2), floor(reads / 2)))
  sf <- data.frame(Name = tx, Length = 1000, EffectiveLength = 900,
                   TPM = round(num / sum(num) * 1e6, 4), NumReads = num)
  write.table(sf, file.path(sdir, "quant.sf"), sep = "\t",
              quote = FALSE, row.names = FALSE)
}

# 6) RangedSummarizedExperiment (.rds) carrying counts + colData.
coldata <- DataFrame(sample = samples,
                     condition = factor(c("control", "control", "treated", "treated")),
                     batch = factor(c("A", "B", "A", "B")),
                     row.names = samples)
se <- SummarizedExperiment(assays = list(counts = m), colData = coldata)
saveRDS(se, file.path(fx, "se.rds"))

# 7) A NORMALIZED (fractional) matrix that must be REJECTED by load_matrix.
mn <- m / rep(colSums(m), each = nrow(m)) * 1e6   # CPM-like, fractional
dfn <- data.frame(gene = rownames(mn), round(mn, 3), check.names = FALSE)
write.csv(dfn, file.path(fx, "normalized.csv"), row.names = FALSE, quote = FALSE)

# 8) Sample sheet matching s1..s4 (for reconcile tests + end-to-end).
meta <- data.frame(sample = samples,
                   condition = c("control", "control", "treated", "treated"),
                   batch = c("A", "B", "A", "B"))
write.csv(meta, file.path(fx, "sample_sheet.csv"), row.names = FALSE, quote = FALSE)

# 9) A mismatched sample sheet (wrong ID) to exercise reconcile FAIL.
meta_bad <- meta; meta_bad$sample[4] <- "s99"
write.csv(meta_bad, file.path(fx, "sample_sheet_bad.csv"), row.names = FALSE, quote = FALSE)

cat("Fixtures written to", fx, "\n")
