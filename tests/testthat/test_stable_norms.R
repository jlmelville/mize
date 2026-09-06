test_that("Euclidean norms avoid intermediate overflow and underflow", {
  expect_equal(norm2(c(3e200, 4e200)) / 1e200, 5)
  expect_equal(norm2(c(3e-200, 4e-200)) / 1e-200, 5)
  expect_equal(norm2(c(-3, 4)), 5)
  expect_identical(norm2(numeric()), 0)
  expect_identical(norm2(c(0, 0)), 0)
  expect_identical(norm2(Inf), Inf)
  expect_true(is.na(norm2(NA_real_)))
  expect_true(is.nan(norm2(NaN)))
  expect_identical(norm2(rep(.Machine$double.xmax, 2)), Inf)
})

test_that("normalization preserves its small-vector policy", {
  tiny <- c(x = 3e-200, y = 4e-200)
  expect_identical(normalize(tiny), tiny)
  expect_identical(normalize(0), 0)
  expect_identical(normalize(.Machine$double.eps / 2), .Machine$double.eps / 2)
  expect_equal(normalize(.Machine$double.eps), 1)
  expect_equal(normalize(c(3e200, 4e200)), c(0.6, 0.8))
})

test_that("public summaries retain representable extreme gradient norms", {
  for (scale in c(1e200, 1e-200)) {
    fg <- list(fn = function(x) sum(x), gr = function(x) c(3, 4) * scale)
    opt <- make_mize(
      method = "SD",
      line_search = "const",
      step0 = 1,
      par = c(0, 0),
      fg = fg
    )
    summary <- mize_step_summary(opt, c(0, 0), fg, calc_gr = TRUE)
    expect_equal(summary$g2n / scale, 5)
  }
})
