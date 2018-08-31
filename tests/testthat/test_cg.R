context("Conjugate Gradient")

test_that("CG with different updates", {
  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "fr", grad_tol = 1e-5)
  expect_equal(res$nf, 20)
  expect_equal(res$ng, 20)
  expect_equal(res$f, 3.117, tol = 1e-3)
  expect_equal(res$g2n, 25.835, tol = 1e-3)
  expect_equal(res$par, c(-0.600, 0.286), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "cd", grad_tol = 1e-5)
  expect_equal(res$nf, 20)
  expect_equal(res$ng, 20)
  expect_equal(res$f, 3.126, tol = 1e-3)
  expect_equal(res$g2n, 25.915, tol = 1e-3)
  expect_equal(res$par, c(-0.603, 0.289), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "dy", grad_tol = 1e-5)
  expect_equal(res$nf, 20)
  expect_equal(res$ng, 20)
  expect_equal(res$f, 3.117, tol = 1e-3)
  expect_equal(res$g2n, 25.845, tol = 1e-3)
  expect_equal(res$par, c(-0.600, 0.286), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "hs", grad_tol = 1e-5)
  expect_equal(res$nf, 21)
  expect_equal(res$ng, 21)
  expect_equal(res$f, 1.829, tol = 1e-3)
  expect_equal(res$g2n, 3.665, tol = 1e-3)
  expect_equal(res$par, c(-0.351, 0.118), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "hs+", grad_tol = 1e-5)
  expect_equal(res$nf, 21)
  expect_equal(res$ng, 21)
  expect_equal(res$f, 1.728, tol = 1e-3)
  expect_equal(res$g2n, 2.350, tol = 1e-3)
  expect_equal(res$par, c(-0.312, 0.106), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "pr", grad_tol = 1e-5)
  expect_equal(res$nf, 21)
  expect_equal(res$ng, 21)
  expect_equal(res$f, 1.552, tol = 1e-3)
  expect_equal(res$g2n, 3.163, tol = 1e-3)
  expect_equal(res$par, c(-0.245, 0.0547), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "pr+", grad_tol = 1e-5)
  expect_equal(res$nf, 21)
  expect_equal(res$ng, 21)
  expect_equal(res$f, 1.475, tol = 1e-3)
  expect_equal(res$g2n, 2.525, tol = 1e-3)
  expect_equal(res$par, c(-0.211, 0.054), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "ls", grad_tol = 1e-5)
  expect_equal(res$nf, 21)
  expect_equal(res$ng, 21)
  expect_equal(res$f, 1.828, tol = 1e-3)
  expect_equal(res$g2n, 3.666, tol = 1e-3)
  expect_equal(res$par, c(-0.351, 0.118), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "hz", grad_tol = 1e-5)
  expect_equal(res$nf, 19)
  expect_equal(res$ng, 19)
  expect_equal(res$f, 2.377, tol = 1e-3)
  expect_equal(res$g2n, 14.363, tol = 1e-3)
  expect_equal(res$par, c(-0.413, 0.232), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 5,
              cg_update = "hz+", grad_tol = 1e-5)
  expect_equal(res$nf, 19)
  expect_equal(res$ng, 19)
  expect_equal(res$f, 2.487, tol = 1e-3)
  expect_equal(res$g2n, 16.750, tol = 1e-3)
  expect_equal(res$par, c(-0.407, 0.237), tol = 1e-3)
})
