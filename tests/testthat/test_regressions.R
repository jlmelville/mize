test_that("Newton direction uses inverse Hessian functions", {
  quadratic_fg <- list(
    fn = function(x) x[1]^2 + 3 * x[2]^2,
    gr = function(x) c(2 * x[1], 6 * x[2])
  )
  inverse_hessians <- list(
    function(x) c(0.5, 1 / 6),
    function(x) diag(c(0.5, 1 / 6))
  )

  for (hi in inverse_hessians) {
    fg <- quadratic_fg
    fg$hi <- hi

    res <- mize(
      c(2, -3),
      fg,
      method = "NEWTON",
      line_search = "constant",
      step0 = 1,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL
    )

    expect_equal(as.numeric(res$par), c(0, 0), tolerance = 1e-12)
  }
})

test_that("make_mize validates public configuration values", {
  cases <- list(
    list(args = list(memory = 0), pattern = "memory must be > 0"),
    list(
      args = list(nest_q = -0.1),
      pattern = "nest_q must be between 0 and 1"
    ),
    list(
      args = list(nest_burn_in = -1),
      pattern = "nest_burn_in must be non-negative"
    ),
    list(args = list(step_up = 0), pattern = "step_up must be positive"),
    list(
      args = list(step_down = -0.1),
      pattern = "step_down must be between 0 and 1"
    ),
    list(
      args = list(dbd_weight = -0.1),
      pattern = "dbd_weight must be between 0 and 1"
    ),
    list(args = list(c1 = -0.1), pattern = "c1 must be between 0 and 1"),
    list(
      args = list(c1 = 0.5, c2 = 0.4),
      pattern = "c2 must be between c1 and 1"
    ),
    list(
      args = list(ls_max_fn = -1),
      pattern = "ls_max_fn must be non-negative"
    ),
    list(
      args = list(ls_max_gr = -1),
      pattern = "ls_max_gr must be non-negative"
    ),
    list(
      args = list(ls_max_fg = -1),
      pattern = "ls_max_fg must be non-negative"
    ),
    list(
      args = list(ls_max_alpha_mult = 0),
      pattern = "ls_max_alpha_mult must be positive"
    ),
    list(
      args = list(ls_max_alpha = 0),
      pattern = "ls_max_alpha must be positive"
    ),
    list(
      args = list(restart_wait = 0),
      pattern = "restart_wait must be a positive integer"
    ),
    list(
      args = list(step_next_init = 0),
      pattern = "numeric argument for step_next_init must be positive"
    )
  )

  for (case in cases) {
    expect_error(do.call(make_mize, case$args), case$pattern)
  }
})

test_that("BFGS update skips bad curvature pairs", {
  hm <- matrix(c(2, 0.25, 0.25, 1), nrow = 2)

  expect_equal(
    bfgs_update(hm, sm = c(1, 0), ym = c(-1, 0), eps = 0),
    hm
  )
  expect_equal(
    bfgs_update(hm, sm = c(1, 0), ym = c(1e-12, 1), eps = 0),
    hm
  )
})

test_that("L-BFGS memory trimming keeps only the newest updates", {
  state <- list(
    memory = 1,
    sms = list(c(1, 0)),
    yms = list(c(2, 0)),
    rhos = list(0.5)
  )

  state <- lbfgs_memory_update(state, ym = c(0, 2), sm = c(0, 1), eps = 0)

  expect_length(state$sms, 1)
  expect_equal(state$sms[[1]], c(0, 1))
  expect_equal(state$yms[[1]], c(0, 2))
  expect_equal(state$rhos[[1]], 0.5)

  state <- list(
    memory = 2,
    sms = list(c(1, 0), c(0, 1)),
    yms = list(c(2, 0), c(0, 2)),
    rhos = list(0.5, 0.5)
  )

  state <- lbfgs_memory_update(state, ym = c(2, 2), sm = c(1, 1), eps = 0)

  expect_length(state$sms, 2)
  expect_equal(state$sms[[1]], c(0, 1))
  expect_equal(state$sms[[2]], c(1, 1))
  expect_equal(state$yms[[1]], c(0, 2))
  expect_equal(state$yms[[2]], c(2, 2))
  expect_equal(state$rhos[[1]], 0.5)
  expect_equal(state$rhos[[2]], 0.25)
})

test_that("L-BFGS memory update skips bad curvature pairs", {
  state <- list(
    memory = 2,
    sms = list(c(1, 0)),
    yms = list(c(2, 0)),
    rhos = list(0.5)
  )

  expect_equal(
    lbfgs_memory_update(state, ym = c(-1, 0), sm = c(1, 0), eps = 0),
    state
  )
  expect_equal(
    lbfgs_memory_update(state, ym = c(1e-12, 1), sm = c(1, 0), eps = 0),
    state
  )
})

test_that("L-BFGS supports memory of one through the public API", {
  res <- mize(
    rb0,
    rosen_no_hess,
    method = "L-BFGS",
    memory = 1,
    max_iter = 3,
    line_search = "mo",
    c1 = 5e-10,
    c2 = 1e-9,
    step0 = "scipy",
    step_next_init = "q",
    scale_hess = FALSE,
    grad_tol = 1e-5
  )

  expect_true(is.finite(res$f))
})

test_that("store_progress records infinity gradient norm", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "SD",
    max_iter = 1,
    line_search = "const",
    step0 = 0.0001,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = 1e-5,
    store_progress = TRUE
  )

  expect_true("ginfn" %in% names(res$progress))
  expect_false("ginf" %in% names(res$progress))
  expect_true(all(is.finite(res$progress$ginfn)))
})

test_that("Delta-Bar-Delta numeric step0 initializes epsilon vector", {
  opt <- make_mize(method = "DBD", step0 = 0.25)
  opt <- mize_init(opt, rb0, rosenbrock_fg)

  expect_equal(
    opt$stages$gradient_descent$step_size$value,
    rep(0.25, length(rb0))
  )
})

test_that("function-valued momentum schedules may take iteration only", {
  seen <- c()
  schedule <- function(iter) {
    seen <<- c(seen, iter)
    0.5
  }

  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "MOM",
    line_search = "const",
    step0 = 0.001,
    mom_schedule = schedule,
    max_iter = 3,
    check_conv_every = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    store_progress = TRUE
  )

  expect_equal(seen, c(2, 3))
  expect_equal(tail(res$progress$mu, 1), 0.5)
})

test_that("function-valued momentum schedules may still take max_iter", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "MOM",
    line_search = "const",
    step0 = 0.001,
    mom_schedule = function(iter, max_iter) iter / max_iter,
    use_init_mom = TRUE,
    max_iter = 2,
    check_conv_every = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    store_progress = TRUE
  )

  expect_equal(tail(res$progress$mu, 1), 1)
})

test_that("append_depends handles optimizer, stage, and sub-stage dependencies", {
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(),
        step_size = constant_step_size(value = 1)
      )
    )
  )

  opt <- append_depends(opt, new_depends = "opt_dep")
  expect_equal(opt$depends, "opt_dep")

  opt <- append_depends(
    opt,
    stage_type = "gradient_descent",
    new_depends = "stage_dep"
  )
  expect_true("stage_dep" %in% opt$stages$gradient_descent$depends)

  opt <- append_depends(
    opt,
    stage_type = "gradient_descent",
    sub_stage_type = "direction",
    new_depends = "direction_dep"
  )
  expect_true(
    "direction_dep" %in% opt$stages$gradient_descent$direction$depends
  )
})

test_that("representative optimizers register lifecycle dependency hooks", {
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = list(fn = function(x) x[1]^2, gr = function(x) 2 * x)
  )
  expect_true("gradient" %in% names(opt$hooks$gradient_descent$before))

  opt <- make_mize(
    method = "L-BFGS",
    line_search = "constant",
    step0 = 0.1,
    par = rb0,
    fg = rosen_no_hess
  )
  after_step_hooks <- names(opt$hooks$step$after)
  expect_true("gradient_old" %in% after_step_hooks)
  expect_true("update_old" %in% after_step_hooks)

  opt <- make_mize(
    method = "MOM",
    line_search = "constant",
    step0 = 0.001,
    mom_schedule = 0.9,
    restart = "fn",
    par = rb0,
    fg = rosenbrock_fg
  )
  expect_true("fn_curr" %in% names(opt$hooks$step$before))
  expect_true("update_fn_cache" %in% names(opt$hooks$step$after))
  expect_true("save_cache_on_failure" %in% names(opt$hooks$step$after))
  expect_true("fn_new" %in% names(opt$hooks$validation$before))
  expect_true("validate_fn" %in% names(opt$hooks$validation$during))
})

test_that("missing lifecycle dependencies fail loudly", {
  opt <- make_mize(method = "SD", line_search = "constant", step0 = 0.1)
  opt$stages$gradient_descent$depends <- c(
    opt$stages$gradient_descent$depends,
    "gradinet"
  )

  expect_error(
    mize_init(opt, c(1), list(fn = function(x) x[1]^2, gr = function(x) 2 * x)),
    "Missing lifecycle dependency 'gradinet'"
  )
})

test_that("CG helper updates default to the unpreconditioned gradient", {
  gm <- c(1, 2)
  gm_old <- c(2, 1)
  pm_old <- -gm_old

  expect_equal(
    hz_plus_update(gm, gm_old, pm_old),
    hz_plus_update(gm, gm_old, pm_old, wm = gm)
  )
  expect_equal(
    prfr_update(gm, gm_old, pm_old),
    prfr_update(gm, gm_old, pm_old, wm = gm)
  )
})

test_that("relative function convergence handles zero baselines", {
  opt <- list()

  expect_equal(
    check_fn_conv(opt, 1, 0, 0, abs_tol = NULL, rel_tol = 1e-8),
    list(what = "rel_tol", val = 0)
  )
  expect_null(check_fn_conv(opt, 1, 0, 1e-8, abs_tol = NULL, rel_tol = 1e-8))
  expect_null(check_fn_conv(opt, 1, 1e-8, 0, abs_tol = NULL, rel_tol = 1e-8))
})

test_that("best gradient restore does not cache an infinite function value", {
  double_well <- list(
    fn = function(x) {
      (x[1]^2 - 1)^2
    },
    gr = function(x) {
      4 * x[1] * (x[1]^2 - 1)
    }
  )

  res <- mize(
    c(1.2),
    double_well,
    method = "SD",
    line_search = "constant",
    step0 = 0.3,
    max_iter = 3,
    check_conv_every = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = 0
  )

  expect_true(is.finite(res$f))
  expect_equal(res$f, double_well$fn(res$par))
  expect_equal(res$par, 1.028033, tolerance = 1e-6)
})

test_that("mize_step_summary can force gradient norm calculation", {
  quadratic_fg <- list(
    fn = function(x) {
      sum(x^2)
    },
    gr = function(x) {
      2 * x
    }
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1, -2),
    fg = quadratic_fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL
  )
  step <- mize_step(opt, c(1, -2), quadratic_fg)

  default_summary <- mize_step_summary(step$opt, step$par, quadratic_fg)
  expect_null(default_summary$g2n)
  expect_null(default_summary$ginfn)
  expect_equal(default_summary$ng, step$ng)

  forced_summary <- mize_step_summary(
    step$opt,
    step$par,
    quadratic_fg,
    calc_gr = TRUE
  )
  grad <- quadratic_fg$gr(step$par)
  expect_equal(forced_summary$g2n, norm2(grad))
  expect_equal(forced_summary$ginfn, norm_inf(grad))
  expect_equal(forced_summary$ng, step$ng + 1)
})

test_that("verbose requires convergence checks to define logging cadence", {
  expect_error(
    mize(
      rb0,
      rosenbrock_fg,
      method = "SD",
      line_search = "constant",
      step0 = 0.0001,
      max_iter = 1,
      check_conv_every = NULL,
      verbose = TRUE
    ),
    "check_conv_every must be non-NULL if verbose is TRUE"
  )
})
