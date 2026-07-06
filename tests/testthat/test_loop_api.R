test_that("steepest descent with constant step size", {
  opt <- make_mize(method = "SD", line_search = "const", step0 = 0.0001)

  opt <- mize_init(opt, rb0, rosenbrock_fg)
  par <- rb0
  for (iter in 1:3) {
    res <- mize_step(opt, par, rosenbrock_fg)
    par <- res$par
    opt <- res$opt
  }

  expect_equal(res$nf, 0)
  expect_equal(res$ng, 3)
  expect_equal(rosenbrock_fg$fn(par), 12.81, tolerance = 1e-3)
  expect_equal(norm2(rosenbrock_fg$gr(par)), 147.11, tolerance = 1e-3)
  expect_equal(res$par, c(-1.144, 1.023), tolerance = 1e-3)
})

test_that("can initialize in make_mize if par and fg are to hand", {
  opt <- make_mize(
    method = "SD",
    line_search = "const",
    step0 = 0.0001,
    par = rb0,
    fg = rosenbrock_fg
  )

  par <- rb0
  for (iter in 1:3) {
    res <- mize_step(opt, par, rosenbrock_fg)
    par <- res$par
    opt <- res$opt
  }

  expect_equal(res$nf, 0)
  expect_equal(res$ng, 3)
  expect_equal(rosenbrock_fg$fn(par), 12.81, tolerance = 1e-3)
  expect_equal(norm2(rosenbrock_fg$gr(par)), 147.11, tolerance = 1e-3)
  expect_equal(res$par, c(-1.144, 1.023), tolerance = 1e-3)
})

test_that("reinitializing produces the same results", {
  opt <- make_mize(method = "BFGS", line_search = "more-thuente")

  opt <- mize_init(opt, rb0, rosen_no_hess)
  par <- rb0
  for (iter in 1:3) {
    res <- mize_step(opt, par, rosen_no_hess)
    par <- res$par
    opt <- res$opt
  }

  expect_equal(res$nf, 5)
  expect_equal(res$ng, 5)
  expect_equal(rosen_no_hess$fn(par), 4.28, tolerance = 1e-3)
  expect_equal(norm2(rosen_no_hess$gr(par)), 17.29, tolerance = 1e-3)
  expect_equal(res$par, c(-1.048, 1.070), tolerance = 1e-3)

  opt <- mize_init(opt, rb0, rosen_no_hess)
  par <- rb0
  for (iter in 1:3) {
    res <- mize_step(opt, par, rosen_no_hess)
    par <- res$par
    opt <- res$opt
  }

  # nf and ng are remembered
  expect_equal(res$nf, 10)
  expect_equal(res$ng, 10)
  expect_equal(rosen_no_hess$fn(par), 4.28, tolerance = 1e-3)
  expect_equal(norm2(rosen_no_hess$gr(par)), 17.29, tolerance = 1e-3)
  expect_equal(res$par, c(-1.048, 1.070), tolerance = 1e-3)
})

test_that("cache predicates require exact value fields", {
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old), add = TRUE)

  opt <- list(
    cache = list(
      fn_new_iter = 1,
      fn_curr_iter = 1,
      gr_curr_iter = 1
    )
  )

  expect_warning(fn_new <- has_fn_new(opt, 1), NA)
  expect_warning(fn_curr <- has_fn_curr(opt, 1), NA)
  expect_warning(gr_curr <- has_gr_curr(opt, 1), NA)

  expect_false(fn_new)
  expect_false(fn_curr)
  expect_false(gr_curr)

  opt$cache$fn_new <- 1
  opt$cache$fn_curr <- 2
  opt$cache$gr_curr <- c(1, 2)

  expect_warning(fn_new <- has_fn_new(opt, 1), NA)
  expect_warning(fn_curr <- has_fn_curr(opt, 1), NA)
  expect_warning(gr_curr <- has_gr_curr(opt, 1), NA)

  expect_true(fn_new)
  expect_true(fn_curr)
  expect_true(gr_curr)
})

test_that("stateful convergence reports status fields", {
  opt <- make_mize(
    method = "SD",
    line_search = "const",
    step0 = 0.0001,
    max_iter = 1
  )
  opt <- mize_init(opt, rb0, rosenbrock_fg)
  step <- mize_step(opt, rb0, rosenbrock_fg)
  step_info <- mize_step_summary(step$opt, step$par, rosenbrock_fg, rb0)
  opt <- check_mize_convergence(step_info)

  expect_true(opt$is_terminated)
  expect_equal(opt$terminate$what, "max_iter")
  expect_false(opt$converged)
  expect_equal(opt$status, "budget_exhausted")
  expect_true(grepl("max_iter", opt$message, fixed = TRUE))
})
