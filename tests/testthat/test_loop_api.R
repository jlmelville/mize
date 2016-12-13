context("External loop API")
test_that("steepest descent with constant step size", {

  opt <- make_mizer(method = "sd", line_search = "const", step0 = 0.0001)

  opt <- mizer_init(opt, rb0, rosenbrock_fg)
  par <- rb0
  for (iter in 1:3) {
    res <- mizer_step(opt, par, rosenbrock_fg, iter)
    par <- res$par
    opt <- res$opt
  }

  expect_equal(res$nf, 0)
  expect_equal(res$ng, 3)
  expect_equal(rosenbrock_fg$fn(par), 12.81, tol = 1e-3)
  expect_equal(norm2(rosenbrock_fg$gr(par)), 147.11, tol = 1e-3)
  expect_equal(res$par, c(-1.144, 1.023), tol = 1e-3)
})

test_that("can initialize in make_mizer if par and fg are to hand", {

  opt <- make_mizer(method = "sd", line_search = "const", step0 = 0.0001,
                    par = rb0, fg = rosenbrock_fg)

  par <- rb0
  for (iter in 1:3) {
    res <- mizer_step(opt, par, rosenbrock_fg, iter)
    par <- res$par
    opt <- res$opt
  }

  expect_equal(res$nf, 0)
  expect_equal(res$ng, 3)
  expect_equal(rosenbrock_fg$fn(par), 12.81, tol = 1e-3)
  expect_equal(norm2(rosenbrock_fg$gr(par)), 147.11, tol = 1e-3)
  expect_equal(res$par, c(-1.144, 1.023), tol = 1e-3)
})
