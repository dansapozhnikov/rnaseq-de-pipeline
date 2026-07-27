# Unit tests for R/qc.R -- severity mapping, threshold logic, and the gate.

# Minimal config stub carrying just the QC thresholds the checks read.
cfg <- list(qc = list(
  min_library_size = 1e6, min_library_size_fail = 1e5,
  min_genes_detected = 3, min_replicates_per_group = 2,
  min_count = 10, max_pct_genes_filtered = 90, cooks_outlier = TRUE))

test_that("qc_check maps condition + severity to the right status", {
  qc <- qc_init()
  qc_check(qc, "a", "WARN", TRUE,  "ok")
  qc_check(qc, "b", "WARN", FALSE, "not ok")
  qc_check(qc, "c", "FAIL", FALSE, "bad")
  expect_equal(qc$results$status, c("PASS", "WARN", "FAIL"))
})

test_that("library size uses hard floor (FAIL) vs soft min (WARN)", {
  qc <- qc_init()
  # 'low' sums to 5e5: above the 1e5 hard floor (PASS) but below the 1e6 soft
  # minimum (WARN). 'ok' sums to 2e6 (fine). Exercises the two-threshold logic.
  counts <- matrix(c(rep(100000L, 5), rep(400000L, 5)), nrow = 5,
                   dimnames = list(paste0("g", 1:5), c("low", "ok")))
  storage.mode(counts) <- "integer"
  qc_check_library_size(qc, counts, cfg)
  st <- setNames(qc$results$status, qc$results$check)
  expect_equal(unname(st["library_size_floor"]), "PASS")  # 5e5 >= 1e5 floor
  expect_equal(unname(st["library_size_min"]),   "WARN")  # 5e5 <  1e6 soft min
})

test_that("replicates_per_group FAILs when a group has too few", {
  qc <- qc_init()
  meta_bad <- data.frame(sample = paste0("s", 1:3),
                         condition = c("control", "treated", "treated"))
  qc_check_replicates(qc, meta_bad, "condition", cfg)  # control has n=1
  expect_equal(qc$results$status[qc$results$check == "replicates_per_group"], "FAIL")

  qc2 <- qc_init()
  meta_ok <- data.frame(sample = paste0("s", 1:4),
                        condition = c("control", "control", "treated", "treated"))
  qc_check_replicates(qc2, meta_ok, "condition", cfg)
  expect_equal(qc2$results$status[qc2$results$check == "replicates_per_group"], "PASS")
})

test_that("design_identifiable detects a confounded covariate", {
  # batch is perfectly confounded with condition => rank-deficient => FAIL.
  meta_conf <- data.frame(sample = paste0("s", 1:4),
                          batch = c("A", "A", "B", "B"),
                          condition = c("control", "control", "treated", "treated"))
  qc <- qc_init()
  qc_check_design_identifiable(qc, meta_conf, "~ batch + condition", "condition")
  expect_equal(qc$results$status[qc$results$check == "design_identifiable"], "FAIL")

  # A crossed design (batch varies within condition) => full rank => PASS.
  meta_ok <- data.frame(sample = paste0("s", 1:4),
                        batch = c("A", "B", "A", "B"),
                        condition = c("control", "control", "treated", "treated"))
  qc2 <- qc_init()
  qc_check_design_identifiable(qc2, meta_ok, "~ batch + condition", "condition")
  expect_equal(qc2$results$status[qc2$results$check == "design_identifiable"], "PASS")
})

test_that("qc_gate stops on FAIL and downgrades under --force", {
  qc <- qc_init()
  qc_check(qc, "x", "FAIL", FALSE, "bad")
  expect_error(qc_gate(qc, force = FALSE), "QC gate")

  qc2 <- qc_init()
  qc_check(qc2, "x", "FAIL", FALSE, "bad")
  qc_gate(qc2, force = TRUE)
  expect_true(qc2$forced)
  expect_false(any(qc2$results$status == "FAIL"))  # downgraded to WARN
})

test_that("low_count_filter WARNs past the percentage threshold", {
  qc <- qc_init()
  qc_check_low_count_filter(qc, n_before = 1000, n_after = 50, cfg)  # 95% removed
  expect_equal(qc$results$status[qc$results$check == "low_count_filter"], "WARN")
})
