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
