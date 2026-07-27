# Unit tests for R/cache.R -- the resumability checkpoint layer.

test_that("with_cache computes once, then reloads without recomputing", {
  td <- tempfile("cache_"); dir.create(td)
  cache <- new_cache(td, enable = TRUE, version = "test", base_deps = list(1))
  calls <- 0
  compute <- function() { calls <<- calls + 1; list(value = 42) }

  a <- with_cache(cache, "stage1", deps = list("x"), compute = compute)
  b <- with_cache(cache, "stage1", deps = list("x"), compute = compute)

  expect_equal(a$value, 42)
  expect_equal(b$value, 42)
  expect_equal(calls, 1)                       # second call served from cache
  expect_gt(length(list.files(file.path(td, "cache"), pattern = "^stage1\\.")), 0)
})

test_that("changing deps invalidates the cache (recompute)", {
  td <- tempfile("cache_"); dir.create(td)
  cache <- new_cache(td, enable = TRUE, version = "test", base_deps = list(1))
  calls <- 0
  compute <- function() { calls <<- calls + 1; calls }

  with_cache(cache, "s", deps = list(a = 1), compute = compute)
  with_cache(cache, "s", deps = list(a = 2), compute = compute)   # different dep
  expect_equal(calls, 2)
})

test_that("changing base_deps (inputs) invalidates the cache", {
  td <- tempfile("cache_"); dir.create(td)
  calls <- 0
  compute <- function() { calls <<- calls + 1; calls }
  c1 <- new_cache(td, enable = TRUE, version = "test", base_deps = list(counts = "A"))
  c2 <- new_cache(td, enable = TRUE, version = "test", base_deps = list(counts = "B"))
  with_cache(c1, "s", deps = list(1), compute = compute)
  with_cache(c2, "s", deps = list(1), compute = compute)
  expect_equal(calls, 2)
})

test_that("disabled cache always computes and writes nothing", {
  td <- tempfile("cache_"); dir.create(td)
  cache <- new_cache(td, enable = FALSE, version = "test", base_deps = list(1))
  calls <- 0
  compute <- function() { calls <<- calls + 1; calls }
  with_cache(cache, "s", deps = list(1), compute = compute)
  with_cache(cache, "s", deps = list(1), compute = compute)
  expect_equal(calls, 2)                       # no reuse when disabled
  expect_equal(length(list.files(file.path(td, "cache"))), 0)
})

test_that("verify=FALSE on a cache hit forces recompute", {
  td <- tempfile("cache_"); dir.create(td)
  cache <- new_cache(td, enable = TRUE, version = "test", base_deps = list(1))
  calls <- 0
  compute <- function() { calls <<- calls + 1; list(n = calls) }
  with_cache(cache, "s", deps = list(1), compute = compute)             # writes
  # A verify that always fails should discard the hit and recompute.
  with_cache(cache, "s", deps = list(1), compute = compute, verify = function(o) FALSE)
  expect_equal(calls, 2)
})

test_that("chaining via parent cascades invalidation", {
  td <- tempfile("cache_"); dir.create(td)
  cache <- new_cache(td, enable = TRUE, version = "test", base_deps = list(1))
  calls <- 0
  compute <- function() { calls <<- calls + 1; calls }
  with_cache(cache, "child", deps = list(1), parent = "keyA", compute = compute)
  with_cache(cache, "child", deps = list(1), parent = "keyB", compute = compute)  # parent changed
  expect_equal(calls, 2)
})

test_that("code fingerprint changes when a compute function changes", {
  # cache_code_fingerprint() reads the compute functions from the global env
  # (where the pipeline sources them), so mutate there and restore afterwards.
  ge <- globalenv()
  had <- exists("build_dds", envir = ge, inherits = FALSE)
  old <- if (had) get("build_dds", envir = ge) else NULL
  on.exit(if (had) assign("build_dds", old, envir = ge) else
            suppressWarnings(rm("build_dds", envir = ge)), add = TRUE)

  assign("build_dds", function() 1, envir = ge); fp1 <- cache_code_fingerprint()
  assign("build_dds", function() 2, envir = ge); fp2 <- cache_code_fingerprint()
  expect_false(identical(fp1, fp2))
})
