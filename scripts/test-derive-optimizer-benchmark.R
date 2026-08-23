#!/usr/bin/env Rscript

source(file.path("scripts", "derive-optimizer-benchmark.R"))

resource_summary_fixture <- function() {
  grid <- expand.grid(
    case = "fixture",
    profile = c(
      "mize-newton",
      "mize-newton-safeguarded",
      "mize-tn",
      "mize-tn-l-bfgs",
      "mize-bfgs",
      "mize-l-bfgs",
      "mize-sr1",
      "mize-cg-pr+",
      "mize-cg-hz+",
      "mize-cg-pr+-l-bfgs"
    ),
    line_search = c(
      "More-Thuente",
      "Rasmussen",
      "Schmidt",
      "Hager-Zhang"
    ),
    step0 = "default",
    callback_mode = c("combined", "separate"),
    rep = 1:3,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  metrics <- data.frame(
    final_f = 0,
    grad_norm = 0,
    iter = 1L,
    n_fn_call = 1L,
    n_gr_call = 1L,
    n_fg_call = 1L,
    n_hs_call = 0L,
    elapsed_sec = 0,
    objective_reduction_hit = TRUE,
    gradient_reduction_hit = TRUE,
    reference_gap_target_applicable = FALSE,
    reference_gap_hit = NA
  )
  cbind(grid, metrics[rep(1L, nrow(grid)), , drop = FALSE])
}

expect_resource_summary_error <- function(summary, message) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  utils::write.csv(summary, path, row.names = FALSE, na = "")
  error <- tryCatch(
    {
      benchmark_read_resource_summary(path)
      NULL
    },
    error = identity
  )
  stopifnot(inherits(error, "error"))
  stopifnot(identical(conditionMessage(error), message))
}

run_derivation_completeness_checks <- function() {
  fixture <- resource_summary_fixture()
  complete_path <- tempfile(fileext = ".csv")
  on.exit(unlink(complete_path))
  utils::write.csv(fixture, complete_path, row.names = FALSE, na = "")
  stopifnot(nrow(benchmark_read_resource_summary(complete_path)) == 240L)

  missing_cell <- with(
    fixture,
    profile == "mize-newton" &
      line_search == "More-Thuente" &
      callback_mode == "combined"
  )
  stopifnot(sum(missing_cell) == 3L)
  expect_resource_summary_error(
    fixture[!missing_cell, , drop = FALSE],
    "resource summary does not contain the complete Packet 4 cell grid per case"
  )
}

run_derivation_completeness_checks()

message("derive-optimizer-benchmark completeness checks passed")
