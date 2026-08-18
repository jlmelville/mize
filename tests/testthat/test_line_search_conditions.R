condition_quadratic_fg <- function(center = 10) {
  list(
    fn = function(x) {
      (x - center)^2
    },
    gr = function(x) {
      2 * (x - center)
    }
  )
}

condition_setup <- function(fg, x, pv = -fg$gr(x) / abs(fg$gr(x))) {
  step0 <- make_step0(fg, x, pv)
  list(
    pv = pv,
    step0 = step0,
    phi = make_phi_alpha(x, fg, pv, calc_gradient_default = TRUE)
  )
}

test_that("Schmidt Armijo controls reject invalid types and ranges", {
  expect_error(
    new_schmidt_armijo_search(armijo_constant = NA_real_),
    "armijo_constant"
  )
  expect_error(new_schmidt_armijo_search(step_down = "half"), "step_down")
  expect_error(new_schmidt_armijo_search(step_down = 2), "step_down")
  expect_error(
    new_schmidt_armijo_search(max_evaluations = -1),
    "max_evaluations"
  )
  expect_error(
    new_schmidt_armijo_search(parameter_tolerance = Inf),
    "parameter_tolerance"
  )
})

test_that("Armijo backtracking accepts steps with sufficient decrease", {
  cases <- list(
    list(
      name = "f1 small initial step",
      fg = f1,
      x = 0,
      alpha = 1e-1,
      c1 = 0.001
    ),
    list(name = "f2 large initial step", fg = f2, x = 0, alpha = 1e1, c1 = 0.1),
    list(
      name = "f4 very large initial step",
      fg = f4,
      x = 0,
      alpha = 1e3,
      c1 = 0.001
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- new_schmidt_armijo_search(
      armijo_constant = case$c1,
      max_evaluations = 10000
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = case$alpha,
      pm = setup$pv
    )

    expect_true(res$step$alpha > 0, info = case$name)
    expect_true(
      armijo_ok_step(setup$step0, res$step, case$c1),
      info = case$name
    )
  }
})

test_that("Schmidt Armijo accepts internal phi results without parameters", {
  step0 <- list(alpha = 0, f = 1, df = -1, d = -1)
  phi <- function(alpha, calc_gradient = TRUE) {
    list(alpha = alpha, f = 0, df = -1, d = -1)
  }

  res <- new_schmidt_armijo_search(
    armijo_constant = 0.1,
    max_evaluations = 1
  )(
    phi = phi,
    step0 = step0,
    alpha = 0.5,
    pm = 1
  )

  expect_equal(res$step$alpha, 0.5)
  expect_equal(res$step$f, 0)
  expect_equal(res$nfn, 1)
  expect_equal(res$ngr, 1)
})

test_that("Schmidt cubic Armijo backs off a nonfinite objective", {
  evaluated <- numeric()
  step0 <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    if (alpha >= 1) {
      return(list(alpha = alpha, f = Inf, df = Inf, d = Inf, par = alpha))
    }
    list(alpha = alpha, f = 0, df = 0, d = 0, par = alpha)
  }

  result <- new_schmidt_armijo_search(
    armijo_constant = 0.1,
    max_evaluations = 2
  )(
    phi = phi,
    step0 = step0,
    alpha = 1,
    pm = 1
  )

  expect_equal(evaluated, c(1, 0.5))
  expect_equal(result$step$alpha, 0.5)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})

test_that("Schmidt cubic Armijo uses values when a trial gradient is nonfinite", {
  evaluated <- numeric()
  step0 <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    if (alpha == 1) {
      return(list(alpha = alpha, f = 2, df = Inf, d = Inf, par = alpha))
    }
    list(alpha = alpha, f = 0, df = 0, d = 0, par = alpha)
  }

  result <- new_schmidt_armijo_search(
    armijo_constant = 0.1,
    max_evaluations = 2
  )(
    phi = phi,
    step0 = step0,
    alpha = 1,
    pm = 1
  )

  expect_equal(evaluated, c(1, 0.25))
  expect_equal(result$step$alpha, 0.25)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})

test_that("More-Thuente successful steps satisfy strong Wolfe conditions", {
  cases <- list(
    list(
      name = "f1 small initial step",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1
    ),
    list(
      name = "f2 small initial step",
      fg = f2,
      x = 0,
      alpha = 1e-3,
      c1 = 0.1,
      c2 = 0.1
    ),
    list(
      name = "f5 large initial step",
      fg = f5,
      x = 0,
      alpha = 1e1,
      c1 = 0.001,
      c2 = 0.001
    ),
    list(
      name = "function modification",
      fg = f6,
      x = 1,
      alpha = 1,
      c1 = 0.1,
      c2 = 0.9
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    search <- new_wolfe_line_search(
      more_thuente_core,
      armijo_constant = case$c1,
      curvature_constant = case$c2,
      max_evaluations = 10000,
      method_policy = new_more_thuente_policy()
    )
    res <- search(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = case$alpha,
      pm = setup$pv
    )
    conditions <- new_line_condition_policy(case$c1, case$c2)

    expect_true(
      conditions$wolfe(setup$step0, res$step),
      info = case$name
    )
  }
})

test_that("Rasmussen and Schmidt successful steps satisfy strong Wolfe conditions", {
  cases <- list(
    list(
      name = "rasmussen",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1,
      run = function(setup, case) {
        new_rasmussen_wolfe_search(
          armijo_constant = case$c1,
          curvature_constant = case$c2,
          max_evaluations = 10000
        )(
          phi = setup$phi,
          step0 = setup$step0,
          alpha = case$alpha,
          pm = setup$pv
        )
      }
    ),
    list(
      name = "schmidt",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1,
      run = function(setup, case) {
        new_schmidt_wolfe_search(
          armijo_constant = case$c1,
          curvature_constant = case$c2,
          max_evaluations = 10000
        )(
          phi = setup$phi,
          step0 = setup$step0,
          alpha = case$alpha,
          pm = setup$pv
        )
      }
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- case$run(setup, case)
    conditions <- new_line_condition_policy(case$c1, case$c2)

    expect_true(
      conditions$wolfe(setup$step0, res$step),
      info = case$name
    )
  }
})

test_that("weak Wolfe configuration does not require strong curvature", {
  c1 <- 1e-4
  c2 <- 0.1
  alpha <- 12.5
  setup <- condition_setup(condition_quadratic_fg(), x = 0, pv = 1)
  weak_conditions <- new_line_condition_policy(
    c1,
    c2,
    strong_curvature = FALSE
  )
  strong_conditions <- new_line_condition_policy(c1, c2)

  results <- list(
    `more-thuente` = new_wolfe_line_search(
      more_thuente_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 100,
      strong_curvature = FALSE,
      method_policy = new_more_thuente_policy()
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    ),
    rasmussen = new_rasmussen_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 100,
      strong_curvature = FALSE
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    ),
    schmidt = new_schmidt_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 100,
      strong_curvature = FALSE
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    )
  )

  for (name in names(results)) {
    step <- results[[name]]$step

    expect_equal(step$alpha, alpha, info = name)
    expect_true(weak_conditions$wolfe(setup$step0, step), info = name)
    expect_false(strong_conditions$wolfe(setup$step0, step), info = name)
  }
})

test_that("Hager-Zhang returned steps satisfy approximate weak Wolfe conditions", {
  cases <- list(
    list(
      name = "f1 small initial step",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1
    ),
    list(
      name = "f2 large initial step",
      fg = f2,
      x = 0,
      alpha = 1e3,
      c1 = 0.1,
      c2 = 0.1
    ),
    list(
      name = "f5 large initial step",
      fg = f5,
      x = 0,
      alpha = 1e1,
      c1 = 0.001,
      c2 = 0.001
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    conditions <- new_line_condition_policy(
      armijo_constant = case$c1,
      curvature_constant = case$c2,
      approximate_armijo = TRUE,
      strong_curvature = FALSE,
      approximation_tolerance = 1e-6
    )
    search <- make_hager_zhang_search(
      armijo_constant = case$c1,
      curvature_constant = case$c2,
      max_evaluations = 100,
      strong_curvature = FALSE,
      approximate_armijo = TRUE
    )
    res <- search(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = case$alpha,
      pm = setup$pv
    )

    expect_true(
      conditions$wolfe(setup$step0, res$step),
      info = case$name
    )
  }
})

test_that("Hager-Zhang initial exhaustion distinguishes acceptance from fallback", {
  quartic_start <- list(alpha = 0, f = 1, df = 4, d = -16)
  quartic_trial <- list(alpha = 1, f = 81, df = -108, d = 432)

  expect_true(curvature_ok_step(quartic_start, quartic_trial, c2 = 0.5))
  expect_false(strong_curvature_ok_step(quartic_start, quartic_trial, c2 = 0.5))
  expect_false(armijo_ok_step(quartic_start, quartic_trial, c1 = 0.1))
  weak_conditions <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = FALSE,
    approximation_tolerance = 1e-6
  )
  expect_false(
    weak_conditions$wolfe(quartic_start, quartic_trial)
  )

  linear_start <- list(alpha = 0, f = 1, df = -1, d = -1)
  linear_trial <- list(alpha = 1, f = 0, df = -1, d = -1)

  expect_true(armijo_ok_step(linear_start, linear_trial, c1 = 0.1))
  expect_false(curvature_ok_step(linear_start, linear_trial, c2 = 0.5))
  expect_false(weak_conditions$wolfe(linear_start, linear_trial))
})

test_that("Hager-Zhang repairs a high-value negative-slope bound", {
  make_step <- function(alpha, f, d) {
    list(alpha = alpha, f = f, d = d, df = d)
  }

  initial_step <- make_step(alpha = 0, f = 0, d = -1)
  bracket <- make_hager_zhang_bracket(
    lower_endpoint = make_step(alpha = 1, f = -1, d = -1),
    upper_endpoint = make_step(alpha = 8, f = 2, d = 1)
  )
  trial_step <- make_step(alpha = 8, f = 2, d = -1)
  evaluate <- function(alpha, calc_gradient = TRUE) {
    if (alpha > 4) {
      make_step(alpha, f = 1, d = -1)
    } else if (alpha < 3) {
      make_step(alpha, f = -1, d = -1)
    } else {
      make_step(alpha, f = -0.5, d = 1)
    }
  }
  evaluator <- new_line_evaluator(evaluate, max_evaluations = 5)
  condition_policy <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = TRUE
  )

  result <- update_hager_zhang_bracket(
    bracket = bracket,
    trial_step = trial_step,
    initial_step = initial_step,
    evaluator = evaluator,
    approximate_decrease_tolerance = 0,
    condition_policy = condition_policy,
    method_policy = make_hager_zhang_policy()
  )

  expect_true(result$succeeded)
  expect_identical(environment(evaluator)$evaluation_count, 3L)
  expect_equal(
    vapply(result$bracket, `[[`, numeric(1), "alpha"),
    c(lower_endpoint = 2.75, upper_endpoint = 3.625)
  )
  expect_equal(
    vapply(result$bracket, `[[`, numeric(1), "d"),
    c(lower_endpoint = -1, upper_endpoint = 1)
  )
  expect_true(result$bracket$lower_endpoint$f <= initial_step$f)
})

test_that("Hager-Zhang retains bracket repair completed at budget exhaustion", {
  make_step <- function(alpha, f, d) {
    list(alpha = alpha, f = f, d = d, df = d)
  }
  initial_step <- make_step(alpha = 0, f = 1, d = -1)
  bracket <- make_hager_zhang_bracket(
    lower_endpoint = make_step(alpha = 1, f = 1 + 5e-7, d = -1),
    upper_endpoint = make_step(alpha = 5, f = 2, d = 4)
  )
  trial_step <- make_step(alpha = 1.8, f = 2, d = -1)
  evaluator <- new_line_evaluator(
    function(alpha, calc_gradient = TRUE) {
      expect_equal(alpha, 1.4)
      make_step(alpha, f = 0, d = -1)
    },
    max_evaluations = 1
  )
  condition_policy <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.1,
    approximate_armijo = TRUE,
    strong_curvature = FALSE
  )

  result <- refine_hager_zhang_bracket_with_secants(
    bracket = bracket,
    trial_step = trial_step,
    initial_step = initial_step,
    evaluator = evaluator,
    approximate_decrease_tolerance = 1e-6,
    condition_policy = condition_policy,
    method_policy = make_hager_zhang_policy()
  )

  expect_false(result$succeeded)
  expect_identical(result$termination_reason, "budget_exhausted")
  expect_identical(environment(evaluator)$evaluation_count, 1L)
  expect_equal(
    vapply(result$bracket, `[[`, numeric(1), "alpha"),
    c(lower_endpoint = 1.4, upper_endpoint = 1.8)
  )
})

test_that("Hager-Zhang bracket updates have explicit trial classifications", {
  make_step <- function(alpha, f, d) {
    list(alpha = alpha, f = f, d = d, df = d)
  }
  initial_step <- make_step(alpha = 0, f = 0, d = -1)
  bracket <- make_hager_zhang_bracket(
    lower_endpoint = make_step(alpha = 1, f = -1, d = -1),
    upper_endpoint = make_step(alpha = 8, f = 1, d = 1)
  )
  cases <- list(
    nonfinite_trial = make_step(alpha = 4, f = Inf, d = NaN),
    outside_bracket = make_step(alpha = 9, f = -1, d = -1),
    upper_endpoint = make_step(alpha = 4, f = -0.5, d = 0),
    lower_endpoint = make_step(alpha = 4, f = -0.5, d = -1),
    needs_bisection = make_step(alpha = 4, f = 0.5, d = -1)
  )

  for (classification in names(cases)) {
    expect_identical(
      classify_hager_zhang_trial(
        bracket,
        cases[[classification]],
        initial_step,
        approximate_decrease_tolerance = 0
      ),
      classification
    )
  }
})

test_that("line searches return the initial step when evaluation budgets are exhausted", {
  c1 <- 1e-4
  c2 <- 0.1
  alpha <- 12.5
  setup <- condition_setup(condition_quadratic_fg(), x = 0, pv = 1)
  weak_conditions <- new_line_condition_policy(
    c1,
    c2,
    strong_curvature = FALSE
  )

  wolfe_searches <- list(
    `more-thuente` = new_wolfe_line_search(
      more_thuente_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      method_policy = new_more_thuente_policy(
        alpha_max = Inf,
        safeguard_cubic = FALSE
      )
    ),
    rasmussen = new_rasmussen_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      interior_fraction = 0.1,
      expansion_factor = 3,
      relative_interval_tolerance = 1e-6
    ),
    schmidt = new_schmidt_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0
    ),
    `hager-zhang` = make_hager_zhang_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0
    )
  )

  for (name in names(wolfe_searches)) {
    res <- wolfe_searches[[name]](
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    )

    expect_equal(res$nfn, 0, info = name)
    expect_equal(res$ngr, 0, info = name)
    expect_equal(res$step$alpha, 0, info = name)
    expect_false(weak_conditions$wolfe(setup$step0, res$step), info = name)
  }

  armijo_res <- new_schmidt_armijo_search(
    armijo_constant = c1,
    max_evaluations = 0
  )(
    phi = setup$phi,
    step0 = setup$step0,
    alpha = alpha,
    pm = setup$pv
  )

  expect_equal(armijo_res$nfn, 0)
  expect_equal(armijo_res$ngr, 0)
  expect_equal(armijo_res$step$alpha, 0)
})

test_that("finite-value guard backs off non-finite line-search evaluations", {
  fg <- list(
    fn = function(x) {
      if (x > 1) {
        return(Inf)
      }
      (x - 0.25)^2
    },
    gr = function(x) {
      if (x > 1) {
        return(Inf)
      }
      2 * (x - 0.25)
    }
  )
  setup <- condition_setup(fg, x = 0, pv = 1)

  res <- find_finite(setup$phi, alpha = 4, min_alpha = 0, max_fn = 4)

  expect_true(res$ok)
  expect_lte(res$step$alpha, 1)
  expect_true(step_is_finite(res$step))
})

test_that("finite-value recovery rejects a non-finite initial alpha", {
  phi <- function(alpha, calc_gradient = TRUE) {
    stop("phi should not be called")
  }

  res <- find_finite(phi, alpha = Inf, min_alpha = 0, max_fn = Inf)

  expect_null(res$step)
  expect_identical(res$nfn, 0L)
  expect_false(res$ok)
})

test_that("finite-value recovery stops without representable progress", {
  smallest_positive <- .Machine$double.xmin * .Machine$double.eps
  evaluated_alphas <- numeric()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated_alphas <<- c(evaluated_alphas, alpha)
    if (length(evaluated_alphas) > 1L) {
      stop("recovery repeated an endpoint")
    }
    list(alpha = alpha, f = Inf, df = Inf, d = Inf, par = alpha)
  }

  res <- find_finite(
    phi,
    alpha = smallest_positive,
    min_alpha = 0,
    max_fn = Inf
  )

  expect_equal(evaluated_alphas, smallest_positive)
  expect_equal(res$step$alpha, smallest_positive)
  expect_identical(res$nfn, 1L)
  expect_false(res$ok)
})
