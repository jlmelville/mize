context("API tests")

# Repeat some of the basic tests, using the consumer API
test_that("steepest descent with constant step size", {
  res <- mizer(rb0, rosenbrock_fg, method = "SD", max_iter = 3,
               line_search = "const", step0 = 0.0001, grad_tol = 1e-5,
               check_conv_every = NULL)

  expect_equal(res$nf, 1)
  expect_equal(res$ng, 4)
  expect_equal(res$f, 12.81, tol = 1e-3)
  expect_equal(res$g2n, 147.11, tol = 1e-3)
  expect_equal(res$par, c(-1.144, 1.023), tol = 1e-3)
})

test_that("grad norm not returned (or calculated) if grad tol is NULL", {
  res <- mizer(rb0, rosenbrock_fg, method = "SD", max_iter = 3,
               line_search = "const", step0 = 0.0001, grad_tol = NULL,
               check_conv_every = NULL)

  expect_equal(res$nf, 1)
  expect_equal(res$ng, 3)
  expect_equal(res$f, 12.81, tol = 1e-3)
  expect_true(is.null(res$g2n))
  expect_equal(res$par, c(-1.144, 1.023), tol = 1e-3)
})


test_that("L-BFGS with More-Thuente LS", {
  # can abbreviate line search name and initializer
  res <- mizer(rb0, rosenbrock_fg, method = "L-BFGS", max_iter = 3,
               line_search = "mo", c1 = 5e-10, c2 = 1e-9, step0 = "s",
               ls_initializer = "q", scale_hess = FALSE, grad_tol = 1e-5)

  expect_equal(res$nf, 17)
  expect_equal(res$ng, 17)
  expect_equal(res$f, 3.53, tol = 1e-3)
  expect_equal(res$g2n, 24.98, tol = 1e-3)
  expect_equal(res$par, c(-0.785, 0.558), tol = 1e-3)
})

test_that("BFGS with More-Thuente LS", {
  res <- mizer(rb0, rosenbrock_fg, method = "BFGS", max_iter = 3,
               line_search = "more-thuente", c1 = 5e-10, c2 = 1e-9, step0 = "s",
               ls_initializer = "quad", scale_hess = FALSE, grad_tol = 1e-5)

  expect_equal(res$nf, 17)
  expect_equal(res$ng, 17)
  expect_equal(res$f, 3.53, tol = 1e-3)
  expect_equal(res$g2n, 24.98, tol = 1e-3)
  expect_equal(res$par, c(-0.785, 0.558), tol = 1e-3)
})

test_that("CG with Rasmussen LS", {
  res <- mizer(rb0, rosenbrock_fg, method = "CG",
               cg_update = "PR+",
               max_iter = 3,
               line_search = "ras", c1 = 5e-10, c2 = 1e-9, step0 = "r",
               ls_initializer = "slope", grad_tol = 1e-5)

  expect_equal(res$nf, 27)
  expect_equal(res$ng, 27)
  expect_equal(res$f, 3.53, tol = 1e-3)
  expect_equal(res$g2n, 24.98, tol = 1e-3)
  expect_equal(res$par, c(-0.785, 0.558), tol = 1e-3)
})

test_that("NAG with Rasmussen LS", {
  res <- mizer(rb0, rosenbrock_fg, method = "NAG",
               nest_convex_approx = FALSE, nest_q = 0, nest_burn_in = 0,
               max_iter = 3,
               line_search = "rasmussen", c1 = 5e-10, c2 = 1e-9,
               step0 = "rasmussen",
               ls_initializer = "slope", grad_tol = 1e-5)

  expect_equal(res$nf, 29)
  expect_equal(res$ng, 29)
  expect_equal(res$f, 3.56, tol = 1e-3)
  expect_equal(res$g2n, 7.2, tol = 1e-3)
  expect_equal(res$par, c(-0.869, 0.781), tol = 1e-3)
})

test_that("bold driver SD and classical momentum", {
  res <- mizer(rb0, rosenbrock_fg,
               method = "SD", norm_direction = TRUE,
               line_search = "bold",
               mom_type = "classical",
               mom_schedule = "ramp", mom_init = 0.1, mom_final = 0.3,
               max_iter = 3, grad_tol = 1e-5)

  expect_equal(res$nf, 11)
  expect_equal(res$ng, 4)
  expect_equal(res$f, 4.37, tol = 1e-3)
  expect_equal(res$g2n, 24.09, tol = 1e-3)
  expect_equal(res$par, c(-1.043, 1.043), tol = 1e-3)
})

test_that("bold driver SD and nesterov momentum", {
  res <- mizer(rb0, rosenbrock_fg,
               method = "SD", norm_direction = TRUE,
               line_search = "bold",
               mom_type = "nesterov",
               mom_schedule = "ramp", mom_init = 0.1, mom_final = 0.3,
               max_iter = 3, grad_tol = 1e-5)

  expect_equal(res$nf, 11)
  expect_equal(res$ng, 3)
  expect_equal(res$f, 4.09, tol = 1e-3)
  expect_equal(res$g2n, 2.344, tol = 1e-3)
  expect_equal(res$par, c(-1.019, 1.050), tol = 1e-3)
})

test_that("Delta bar delta adaptive learning rate and momentum", {
  res <- mizer(rb0, rosenbrock_fg,
               method = "DBD", norm_direction = TRUE,
               step0 = 0.1,
               mom_type = "constant",
               mom_schedule = 0.2,
               max_iter = 3, grad_tol = 1e-5, check_conv_every = NULL)

  expect_equal(res$nf, 1)
  expect_equal(res$ng, 3)
  expect_equal(res$f, 4.84, tol = 1e-3)
  expect_equal(res$g2n, 37.85, tol = 1e-3)
  expect_equal(res$par, c(-0.993, 1.079), tol = 1e-3)
})

test_that("Terminates semi-gracefully if function value is non-finite", {
  res <- mizer(rb0, rosenbrock_fg, "DBD", step0 = 1, check_conv_every = 1)
  expect_equal(res$terminate$what, "fn_inf")
  expect_equal(res$iter, 4)
})

test_that("Terminates semi-gracefully if gradient is non-finite", {
  # If we don't check convergence often enough, solution can diverge
  # in between checks. If NaN is detected in a gradient calculation, we
  # terminate early even if not on a convergence check iteration
  res <- mizer(rb0, rosenbrock_fg, "DBD", step0 = 1, check_conv_every = 10)
  expect_equal(res$terminate$what, "gr_inf")
  expect_equal(res$iter, 6)
})


