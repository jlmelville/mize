test_that("CG restarts when recurrence arithmetic is unusable", {
  # Isolate coefficient, candidate-vector, and slope overflow before line search.
  cases <- list(
    list(g = 1e200, old = 1e200, p = -1),
    list(g = 1e200, old = 1, p = -1),
    list(g = 1e150, old = 1, p = -1e200),
    list(g = 1e100, old = 1, p = -1e100)
  )
  for (update in list(fr_update, prfr_update)) {
    for (case in cases) {
      direction <- cg_direction(cg_update = update)
      direction$pm_old <- case$p
      opt <- list(
        cache = list(gr_curr = case$g, gr_curr_iter = 2, gr_old = case$old)
      )
      result <- direction$calculate(opt, list(), direction, 0, list(), 2)
      expect_equal(result$sub_stage$value, -case$g)
    }
  }
})

test_that("CG arithmetic fallback is usable or fails explicitly downstream", {
  fg <- list(fn = function(x) sum(x), gr = function(x) rep(1e200, length(x)))
  for (update in c(
    "FR",
    "CD",
    "DY",
    "HS",
    "HS+",
    "PR",
    "PR+",
    "LS",
    "HZ",
    "HZ+",
    "PRFR"
  )) {
    result <- mize(
      0,
      fg,
      method = "CG",
      cg_update = update,
      line_search = "const",
      step0 = 1e-200,
      max_iter = 3,
      abs_tol = NULL,
      rel_tol = NULL,
      step_tol = NULL
    )
    expect_equal(result$par, -3, info = update)
    expect_identical(result$terminate$what, "max_iter", info = update)
    result <- mize(0, fg, method = "CG", cg_update = update, max_iter = 3)
    expect_false(result$converged, info = update)
    expect_identical(result$terminate$what, "line_search_failed", info = update)
    expect_identical(
      result$terminate$val,
      "nonfinite_initial_point",
      info = update
    )
  }
})
