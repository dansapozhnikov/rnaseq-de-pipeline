# Unit tests for R/io_detect.R — one test per supported format plus the
# guardrails (normalized-data rejection, sample reconciliation).

# The reference matrix every fixture encodes (genes x samples), for equality checks.
ref <- matrix(c(100, 0, 50, 200, 10,
                120, 2, 60, 180, 12,
                5,  300, 40, 20, 500,
                8,  280, 45, 25, 520),
              nrow = 5,
              dimnames = list(c("G1","G2","G3","G4","G5"), c("s1","s2","s3","s4")))
storage.mode(ref) <- "integer"

test_that("detect_format identifies each format", {
  expect_equal(detect_format(fixture("matrix.csv")), "matrix")
  expect_equal(detect_format(fixture("featurecounts.txt")), "featurecounts")
  expect_equal(detect_format(fixture("star")), "star")
  expect_equal(detect_format(fixture("htseq")), "htseq")
  expect_equal(detect_format(fixture("salmon")), "salmon")
  expect_equal(detect_format(fixture("se.rds")), "rangedSE")
})

test_that("plain matrix loads to the reference integer matrix", {
  res <- load_matrix(fixture("matrix.csv"))
  expect_true(is.integer(res$counts))
  expect_equal(res$counts[, colnames(ref)], ref)
})

test_that("featureCounts drops annotation cols and keeps gene length", {
  res <- load_featurecounts(fixture("featurecounts.txt"))
  expect_equal(res$counts[, colnames(ref)], ref)
  expect_equal(nrow(res$feature_meta), 5L)
  expect_true(all(res$feature_meta$length > 1000))
})

test_that("STAR uses the unstranded column and drops summary rows", {
  res <- load_star(fixture("star"), strand_col = 2L)
  expect_equal(res$counts[rownames(ref), colnames(ref)], ref)
  expect_false(any(grepl("^N_", rownames(res$counts))))
})

test_that("HTSeq drops __* rows", {
  res <- load_htseq(fixture("htseq"))
  expect_equal(res$counts[rownames(ref), colnames(ref)], ref)
  expect_false(any(grepl("^__", rownames(res$counts))))
})

test_that("Salmon summarises transcripts to gene level and rounds", {
  res <- load_salmon(fixture("salmon"), tx2gene = fixture("tx2gene.csv"))
  expect_true(is.integer(res$counts))
  # Two transcripts split each gene's count via ceil/floor, so gene totals match.
  expect_equal(res$counts[rownames(ref), colnames(ref)], ref)
})

test_that("rangedSE returns counts AND colData", {
  res <- load_rangedSE(fixture("se.rds"))
  expect_equal(res$counts[rownames(ref), colnames(ref)], ref)
  expect_true(!is.null(res$coldata))
  expect_true("condition" %in% colnames(res$coldata))
})

test_that("normalized (fractional) input is REJECTED, not silently coerced", {
  expect_error(load_matrix(fixture("normalized.csv")), "normalized|raw")
})

test_that("reconcile_samples fails loudly and prints offending IDs", {
  res  <- load_matrix(fixture("matrix.csv"))
  meta_ok  <- read.csv(fixture("sample_sheet.csv"))
  meta_bad <- read.csv(fixture("sample_sheet_bad.csv"))
  # Good sheet reorders cleanly.
  out <- reconcile_samples(res$counts, meta_ok, "sample")
  expect_equal(colnames(out), as.character(meta_ok$sample))
  # Bad sheet (s99) triggers a mismatch error mentioning the IDs.
  expect_error(reconcile_samples(res$counts, meta_bad, "sample"), "s99|s4")
})
