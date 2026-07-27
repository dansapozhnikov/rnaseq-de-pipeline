# Unit tests for R/de.R result-spec resolution (Wald/LRT/multi-contrast/interaction).

test_that("primary contrast yields a single spec with the right name", {
  cfg <- list(contrast = list("condition", "treated", "control"), de = list())
  specs <- resolve_result_specs(cfg)
  expect_length(specs, 1)
  expect_equal(specs[[1]]$name, "condition_treated_vs_control")
  expect_equal(specs[[1]]$contrast, cfg$contrast)
  expect_null(specs[[1]]$coef)
})

test_that("a plain [factor,num,denom] extra contrast is appended", {
  cfg <- list(contrast = list("condition", "treated", "control"),
              de = list(contrasts = list(list("condition", "other", "control"))))
  specs <- resolve_result_specs(cfg)
  expect_length(specs, 2)
  expect_equal(specs[[2]]$name, "condition_other_vs_control")
  expect_equal(specs[[2]]$contrast[[2]], "other")
})

test_that("a {name, coef} extra contrast becomes a coefficient spec (interaction)", {
  cfg <- list(contrast = list("condition", "treated", "control"),
              de = list(contrasts = list(list(name = "interaction",
                                              coef = "conditiontreated.batchB"))))
  specs <- resolve_result_specs(cfg)
  expect_length(specs, 2)
  expect_equal(specs[[2]]$name, "interaction")
  expect_equal(specs[[2]]$coef, "conditiontreated.batchB")
  expect_null(specs[[2]]$contrast)
})

test_that("duplicate spec names are de-duplicated (primary wins)", {
  cfg <- list(contrast = list("condition", "treated", "control"),
              de = list(contrasts = list(list("condition", "treated", "control"))))
  specs <- resolve_result_specs(cfg)
  expect_length(specs, 1)  # the extra duplicates the primary and is dropped
})

test_that("load_config rejects LRT without a reduced model", {
  tmp <- tempfile(fileext = ".yaml")
  writeLines(c(
    "seed: 1", "organism: human", "input_format: matrix",
    "design: \"~ batch + condition\"", "contrast: [condition, treated, control]",
    "qc: {min_library_size: 1, min_library_size_fail: 1, min_genes_detected: 1,",
    "     min_replicates_per_group: 2, min_count: 1, max_pct_genes_filtered: 90}",
    "de: {padj_cutoff: 0.05, lfc_cutoff: 1, shrink: apeglm, test: LRT}",
    "enrichment: {enable: false, ontology: BP}",
    "outputs: {outdir: results/}"), tmp)
  expect_error(load_config(tmp), "LRT.*reduced|reduced.*LRT")
})
