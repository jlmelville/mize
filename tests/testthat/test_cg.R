context("Conjugate Gradient")

test_that("CG with different updates", {
  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "fr", grad_tol = 1e-5)
  expect_equal(res$nf, 14)
  expect_equal(res$ng, 14)
  expect_equal(res$f, 3.504, tol = 1e-3)
  expect_equal(res$g2n, 20.226, tol = 1e-3)
  expect_equal(res$par, c(-0.819, 0.626), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "cd", grad_tol = 1e-5)
  expect_equal(res$nf, 14)
  expect_equal(res$ng, 14)
  expect_equal(res$f, 3.512, tol = 1e-3)
  expect_equal(res$g2n, 20.274, tol = 1e-3)
  expect_equal(res$par, c(-0.821, 0.629), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "dy", grad_tol = 1e-5)
  expect_equal(res$nf, 10)
  expect_equal(res$ng, 10)
  expect_equal(res$f, 4.112, tol = 1e-3)
  expect_equal(res$g2n, 2.579, tol = 1e-3)
  expect_equal(res$par, c(-1.028, 1.060), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "hs", grad_tol = 1e-5)
  expect_equal(res$nf, 13)
  expect_equal(res$ng, 13)
  expect_equal(res$f, 3.531, tol = 1e-3)
  expect_equal(res$g2n, 24.935, tol = 1e-3)
  expect_equal(res$par, c(-0.785, 0.558), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "hs+", grad_tol = 1e-5)
  expect_equal(res$nf, 13)
  expect_equal(res$ng, 13)
  expect_equal(res$f, 3.531, tol = 1e-3)
  expect_equal(res$g2n, 24.935, tol = 1e-3)
  expect_equal(res$par, c(-0.785, 0.558), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "pr", grad_tol = 1e-5)
  expect_equal(res$nf, 13)
  expect_equal(res$ng, 13)
  expect_equal(res$f, 3.455, tol = 1e-3)
  expect_equal(res$g2n, 24.32, tol = 1e-3)
  expect_equal(res$par, c(-0.765, 0.527), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "pr+", grad_tol = 1e-5)
  expect_equal(res$nf, 13)
  expect_equal(res$ng, 13)
  expect_equal(res$f, 3.455, tol = 1e-3)
  expect_equal(res$g2n, 24.32, tol = 1e-3)
  expect_equal(res$par, c(-0.765, 0.527), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "ls", grad_tol = 1e-5)
  expect_equal(res$nf, 9)
  expect_equal(res$ng, 9)
  expect_equal(res$f, 4.122, tol = 1e-3)
  expect_equal(res$g2n, 1.833, tol = 1e-3)
  expect_equal(res$par, c(-1.029, 1.066), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "hz", grad_tol = 1e-5)
  expect_equal(res$nf, 14)
  expect_equal(res$ng, 14)
  expect_equal(res$f, 3.484, tol = 1e-3)
  expect_equal(res$g2n, 19.505, tol = 1e-3)
  expect_equal(res$par, c(-0.817, 0.626), tol = 1e-3)

  res <- mize(rb0, rosenbrock_fg, method = "cg", max_iter = 3,
              cg_update = "hz+", grad_tol = 1e-5)
  expect_equal(res$nf, 14)
  expect_equal(res$ng, 14)
  expect_equal(res$f, 3.490, tol = 1e-3)
  expect_equal(res$g2n, 19.655, tol = 1e-3)
  expect_equal(res$par, c(-0.818, 0.626), tol = 1e-3)
})
