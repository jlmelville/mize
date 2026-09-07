test_that("Bold Driver no allowance is a failed search in both public routes", {
  fg <- list(fn = function(x) sum(x^2), gr = function(x) 2 * x)
  for (store in c(FALSE, TRUE)) {
    result <- mize(
      1,
      fg,
      method = "SD",
      line_search = "bold",
      ls_max_fn = 0,
      max_iter = 1,
      store_progress = store
    )
    expect_false(result$converged)
    expect_identical(result$terminate$what, "line_search_failed")
    expect_identical(result$terminate$val, "budget_exhausted")
    if (store) expect_equal(tail(result$progress$alpha, 1), 0)
  }
  opt <- make_mize(method = "SD", line_search = "bold", ls_max_fn = 0)
  opt <- mize_init(opt, 1, fg)
  step <- mize_step(opt, 1, fg)
  summary <- mize_step_summary(step$opt, step$par, fg, par_old = 1)
  opt <- check_mize_convergence(summary)
  expect_identical(opt$terminate$what, "line_search_failed")
  expect_equal(summary$alpha, 0)
})

test_that("Bold Driver reports completed alpha and preserves proposal schedule", {
  fg <- list(fn = function(x) x^2 / 4, gr = function(x) x / 2)
  opt <- make_mize(method = "SD", line_search = "bold")
  par <- 2
  opt <- mize_init(opt, par, fg)
  for (i in 1:4) {
    old <- par
    step <- mize_step(opt, par, fg)
    opt <- step$opt
    par <- step$par
    summary <- mize_step_summary(opt, par, fg, par_old = old)
    expected_alpha <- 1.1^(i - 1)
    expect_equal(par, old * (1 - expected_alpha / 2))
    expect_equal(summary$alpha, expected_alpha)
    expect_equal(summary$alpha_init, expected_alpha)
    expect_equal(opt$stages$gradient_descent$step_size$value, 1.1^i)
  }
})

test_that("Bold Driver zero gradient alpha allows a later momentum update", {
  fg <- list(
    fn = function(x) x^2 / 4,
    gr = function(x) if (x == 1) 0 else x / 2
  )
  opt <- make_mize(method = "MOM", line_search = "bold", mom_schedule = 0.9)
  opt <- mize_init(opt, 2, fg)
  first <- mize_step(opt, 2, fg)
  expect_equal(first$par, 1)
  second <- mize_step(first$opt, first$par, fg)
  summary <- mize_step_summary(second$opt, second$par, fg, par_old = first$par)
  expect_equal(summary$alpha, 0)
  expect_equal(summary$step, 0.9)
  expect_equal(second$par, 0.1)
})

test_that("Bold Driver honors step0 and adapts it through both public routes", {
  fg <- list(fn = function(x) x^2 / 4, gr = function(x) x / 2)
  for (step0 in list(NULL, 0.01, 5)) {
    # On this quadratic, 5 is rejected and halved to 2.5; smaller trials pass.
    first_alpha <- if (is.null(step0)) 1 else min(step0, 2.5)
    opt <- make_mize(method = "SD", line_search = "bold", step0 = step0)
    par <- 2
    opt <- mize_init(opt, par, fg)
    for (i in 1:4) {
      old <- par
      step <- mize_step(opt, par, fg)
      opt <- step$opt
      par <- step$par
      summary <- mize_step_summary(opt, par, fg, par_old = old)
      expected_alpha <- first_alpha * 1.1^(i - 1)
      expected_init <- if (i == 1 && !is.null(step0)) step0 else expected_alpha
      expect_equal(par, old * (1 - expected_alpha / 2))
      expect_equal(summary$alpha, expected_alpha)
      expect_equal(summary$alpha_init, expected_init)
    }

    result <- mize(
      2,
      fg,
      method = "SD",
      line_search = "bold",
      step0 = step0,
      max_iter = 4,
      store_progress = TRUE
    )
    expect_equal(result$par, par)
    expect_equal(tail(result$progress$alpha, 4), first_alpha * 1.1^(0:3))
  }
})

test_that("Bold Driver validates step0 before calling the objective", {
  invalid_steps <- list(
    0,
    -1,
    Inf,
    -Inf,
    NA_real_,
    NaN,
    numeric(),
    c(1, 2),
    "rasmussen",
    TRUE,
    1 + 1i,
    list(1),
    matrix(1)
  )
  fg <- list(
    fn = function(x) stop("objective called"),
    gr = function(x) stop("gradient called")
  )
  for (step0 in invalid_steps) {
    expect_error(
      make_mize(method = "SD", line_search = "bold", step0 = step0),
      "step0 must be a positive finite numeric scalar"
    )
    expect_error(
      mize(2, fg, method = "SD", line_search = "bold", step0 = step0),
      "step0 must be a positive finite numeric scalar"
    )
  }
})
