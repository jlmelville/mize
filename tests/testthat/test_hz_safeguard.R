test_that("HZ+ clips using both Euclidean norms in equation 8", {
  # Direct algebra isolates the safeguard from line-search trajectory changes.
  cases <- list(
    list(g = c(-100, 1), old = c(2, -1), p = c(-10, 0)),
    list(g = c(-0.08, -0.07), old = c(0.003, -0.004), p = c(8, -10))
  )
  for (case in cases) {
    y <- case$g - case$old
    py <- sum(case$p * y)
    beta <- sum((y - 2 * case$p * sum(y^2) / py) * case$g) / py
    lower <- -1 /
      (sqrt(sum(case$p^2)) *
        min(0.01, sqrt(sum(case$old^2))))
    expect_lt(beta, lower)
    expect_equal(hz_update(case$g, case$old, case$p, eps = 0), beta)
    expect_equal(hz_plus_update(case$g, case$old, case$p, eps = 0), lower)
  }
  # Also retain an unrestricted negative coefficient when it meets the bound.
  expect_equal(
    hz_plus_update(c(-1, 1), c(2, -1), c(-10, 0), eps = 0),
    hz_update(c(-1, 1), c(2, -1), c(-10, 0), eps = 0)
  )
})
