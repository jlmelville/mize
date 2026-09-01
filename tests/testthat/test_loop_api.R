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

test_that("mize_init explicit convergence arguments override factory values", {
  opt <- make_mize(
    max_iter = 100,
    max_fn = 101,
    max_gr = 102,
    max_fg = 103,
    abs_tol = 0.9,
    rel_tol = 0.8,
    grad_tol = 0.7,
    ginf_tol = 0.6,
    step_tol = 0.5
  )

  opt <- mize_init(
    opt,
    rb0,
    rosenbrock_fg,
    max_iter = 10,
    max_fn = 11,
    max_gr = 12,
    max_fg = 13,
    abs_tol = 0.4,
    rel_tol = 0.3,
    grad_tol = 0.2,
    ginf_tol = 0.1,
    step_tol = 0.05
  )

  expect_equal(
    opt$convergence,
    list(
      max_iter = 10,
      max_fn = 11,
      max_gr = 12,
      max_fg = 13,
      abs_tol = 0.4,
      rel_tol = 0.3,
      grad_tol = 0.2,
      ginf_tol = 0.1,
      step_tol = 0.05
    )
  )
})

test_that("mize_init distinguishes omitted, infinite, and null arguments", {
  opt <- make_mize(max_iter = 12, abs_tol = 0.2, rel_tol = 0.3)
  initialized <- mize_init(opt, rb0, rosenbrock_fg)

  expect_equal(initialized$convergence, opt$convergence)

  initialized <- mize_init(
    initialized,
    rb0,
    rosenbrock_fg,
    max_iter = Inf,
    abs_tol = NULL,
    rel_tol = NULL
  )

  expect_equal(initialized$convergence$max_iter, Inf)
  expect_null(initialized$convergence$abs_tol)
  expect_null(initialized$convergence$rel_tol)

  initialized <- mize_init(
    make_mize(abs_tol = 0.2, rel_tol = 0.3),
    rb0,
    rosenbrock_fg,
    abs_tol = 0.1
  )

  expect_equal(initialized$convergence$abs_tol, 0.1)
  expect_equal(initialized$convergence$rel_tol, 0.1)
})

test_that("reinitializing clears convergence and terminal state", {
  quadratic_fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )

  terminate_by_convergence <- function() {
    opt <- make_mize(
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      grad_tol = 2
    )
    opt <- mize_init(opt, c(1), quadratic_fg)
    step <- mize_step(opt, c(1), quadratic_fg)
    check_mize_convergence(
      mize_step_summary(step$opt, step$par, quadratic_fg, c(1))
    )
  }

  terminate_by_budget <- function() {
    opt <- make_mize(
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      max_iter = 1,
      abs_tol = 0
    )
    opt <- mize_init(opt, c(1), quadratic_fg)
    step <- mize_step(opt, c(1), quadratic_fg)
    check_mize_convergence(
      mize_step_summary(step$opt, step$par, quadratic_fg, c(1))
    )
  }

  terminate_by_nonfinite <- function() {
    nonfinite_fg <- list(
      fn = function(x) Inf,
      gr = function(x) 2 * x
    )
    opt <- make_mize(
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      abs_tol = 0
    )
    opt <- mize_init(opt, c(1), nonfinite_fg)
    step <- mize_step(opt, c(1), nonfinite_fg)
    check_mize_convergence(
      mize_step_summary(step$opt, step$par, nonfinite_fg, c(1))
    )
  }

  terminated_opts <- list(
    terminate_by_convergence(),
    terminate_by_budget(),
    terminate_by_nonfinite()
  )
  expect_equal(
    vapply(terminated_opts, function(opt) opt$status, character(1)),
    c("converged", "budget_exhausted", "failed")
  )

  for (opt in terminated_opts) {
    expect_true(opt$is_terminated)
    opt <- mize_init(opt, c(1), quadratic_fg)

    expect_false(opt$is_terminated)
    expect_null(opt$terminate)
    expect_null(opt$converged)
    expect_null(opt$status)
    expect_null(opt$message)
    expect_null(opt$ok)
    expect_null(opt$convergence$fn_new)

    step <- mize_step(opt, c(1), quadratic_fg)
    expect_equal(step$opt$iter, 1)
    expect_equal(step$par, 0.8)
  }
})

test_that("mize_step rejects an uninitialized optimizer without mutation", {
  opt <- make_mize(method = "SD", line_search = "constant", step0 = 0.1)
  opt_before <- opt
  par <- c(1)

  expect_error(
    mize_step(opt, par, list(fn = function(x) x^2, gr = function(x) 2 * x)),
    "mize_init"
  )
  expect_identical(opt, opt_before)
  expect_identical(par, c(1))
})

test_that("stateful max_iter detects a limit crossed between checks", {
  quadratic_fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    max_iter = 2,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, c(1), quadratic_fg)
  par <- c(1)
  for (i in seq_len(3)) {
    step <- mize_step(opt, par, quadratic_fg)
    opt <- step$opt
    par <- step$par
  }

  opt <- check_mize_convergence(
    mize_step_summary(opt, par, quadratic_fg, calc_fn = FALSE)
  )

  expect_true(opt$is_terminated)
  expect_equal(opt$terminate$what, "max_iter")
  expect_equal(opt$terminate$val, 2)
})

test_that("stateful convergence ignores stale function and gradient caches", {
  quadratic_fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    max_iter = Inf,
    abs_tol = 0.1,
    rel_tol = NULL,
    grad_tol = 0.1,
    step_tol = NULL
  )
  opt <- mize_init(opt, c(1), quadratic_fg)
  opt$convergence$fn_new <- 1
  opt$cache$fn_curr <- 1
  opt$cache$fn_curr_iter <- "invalid"
  opt$cache$gr_curr <- 0
  opt$cache$gr_curr_iter <- "invalid"

  summary <- mize_step_summary(
    opt,
    c(1),
    quadratic_fg,
    calc_fn = FALSE,
    calc_gr = FALSE
  )
  checked <- check_mize_convergence(summary)

  expect_null(checked$terminate)
  expect_false(checked$is_terminated)
  expect_equal(checked$convergence$fn_new, 1)
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

test_that("stateful DBD reports a later non-finite gradient", {
  gradient_calls <- 0L
  fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) {
      gradient_calls <<- gradient_calls + 1L
      if (gradient_calls >= 3L) rep(NaN, length(x)) else 2 * x
    }
  )
  par <- c(1, -1)
  opt <- make_mize(
    method = "DBD",
    step0 = 0.1,
    par = par,
    fg = fg,
    max_iter = 10,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = 1e-6,
    step_tol = NULL
  )

  while (!opt$is_terminated) {
    par_old <- par
    step <- mize_step(opt, par, fg)
    opt <- step$opt
    par <- step$par
    if (opt$is_terminated) {
      break
    }

    summary <- mize_step_summary(opt, par, fg, par_old)
    opt <- summary$opt
    if (!opt$is_terminated) {
      opt <- check_mize_convergence(summary)
    }
  }

  expect_equal(gradient_calls, 3L)
  expect_equal(opt$terminate$what, "gr_inf")
  expect_false(opt$converged)
  expect_equal(opt$status, "failed")
})
