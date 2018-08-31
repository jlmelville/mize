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
