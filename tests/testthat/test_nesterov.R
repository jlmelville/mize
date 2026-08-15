# From the table in https://jlmelville.github.io/mize/nesterov.html
# Uses R code, but nothing from mize
test_that("classical momentum with constant step size", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "MOM",
    line_search = "constant",
    step0 = 0.001,
    mom_schedule = 0.95,
    max_iter = 5,
    store_progress = TRUE,
    use_init_mom = TRUE
  )

  expect_equal(
    res$progress$f,
    c(24.20, 5.35, 25.54, 22.49, 4.32, 18.83),
    tolerance = 1e-3
  )
})

# Use Nesterovized momentum (useful if method is e.g. DBD)
test_that("nesterov momentum with constant step size", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "MOM",
    line_search = "constant",
    step0 = 0.001,
    mom_type = "nesterov",
    mom_schedule = 0.95,
    max_iter = 5,
    store_progress = TRUE,
    use_init_mom = TRUE
  )

  expect_equal(
    res$progress$f,
    c(24.20, 34.96, 7.04, 5.02, 3.85, 3.75),
    tolerance = 1e-3
  )
})

# Use NAG method directly
test_that("NAG with constant step size", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "NAG",
    line_search = "constant",
    step0 = 0.001,
    mom_schedule = 0.95,
    max_iter = 5,
    store_progress = TRUE,
    use_init_mom = TRUE
  )

  expect_equal(
    res$progress$f,
    c(24.20, 34.96, 7.04, 5.02, 3.85, 3.75),
    tolerance = 1e-3
  )
})

test_that("classical momentum  using initial momentum should make no difference", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "MOM",
    line_search = "constant",
    step0 = 0.001,
    mom_schedule = 0.95,
    max_iter = 5,
    store_progress = TRUE,
    use_init_mom = FALSE
  )

  expect_equal(
    res$progress$f,
    c(24.20, 5.35, 25.54, 22.49, 4.32, 18.83),
    tolerance = 1e-3
  )
})

# Use Nesterovized momentum
test_that("nesterov momentum with constant step size force gradient descent first iteration", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "MOM",
    line_search = "constant",
    step0 = 0.001,
    mom_type = "nesterov",
    mom_schedule = 0.95,
    max_iter = 5,
    store_progress = TRUE,
    use_init_mom = FALSE
  )

  expect_equal(
    res$progress$f,
    c(24.20, 5.35, 5.26, 4.13, 4.10, 4.07),
    tolerance = 1e-3
  )
})

# Use NAG method directly
test_that("NAG with constant step size force gradient descent first iteration", {
  res <- mize(
    rb0,
    rosenbrock_fg,
    method = "NAG",
    line_search = "constant",
    step0 = 0.001,
    mom_schedule = 0.95,
    max_iter = 5,
    store_progress = TRUE,
    use_init_mom = FALSE
  )

  expect_equal(
    res$progress$f,
    c(24.20, 5.35, 5.26, 4.13, 4.10, 4.07),
    tolerance = 1e-3
  )
})

test_that("Nesterov-first summaries count every gradient callback", {
  calls <- 0
  fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) {
      calls <<- calls + 1
      2 * x
    }
  )
  opt <- make_opt(make_stages(
    momentum_stage(
      direction = momentum_direction(),
      step_size = constant_step_size(0.5)
    ),
    gradient_stage(
      direction = sd_direction(),
      step_size = constant_step_size(0.1)
    )
  ))
  opt <- mize_init(
    opt,
    c(1),
    fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  first <- mize_step_summary(opt, c(1), fg, calc_gr = TRUE)
  second <- mize_step_summary(first$opt, c(1), fg, calc_gr = TRUE)

  expect_equal(calls, 2)
  expect_equal(first$ng, 1)
  expect_equal(second$ng, 2)
  expect_false(has_gr_curr(second$opt, 1))
})
