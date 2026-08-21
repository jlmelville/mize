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
  initial_point <- make_initial_line_point(fg, x, pv)
  list(
    pv = pv,
    initial_point = initial_point,
    evaluate_line = make_line_function(x, fg, pv, calc_gradient_default = TRUE)
  )
}

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
    res <- make_schmidt_armijo_search(
      armijo_constant = case$c1,
      max_evaluations = 10000
    )(
      evaluate_line = setup$evaluate_line,
      initial_point = setup$initial_point,
      initial_alpha = case$alpha,
      search_direction = setup$pv
    )

    expect_true(res$line_point$alpha > 0, info = case$name)
    expect_true(
      line_point_satisfies_armijo(setup$initial_point, res$line_point, case$c1),
      info = case$name
    )
  }
})

test_that("Schmidt Armijo accepts internal phi results without parameters", {
  initial_point <- list(alpha = 0, value = 1, gradient = -1, slope = -1)
  phi <- function(alpha, calc_gradient = TRUE) {
    list(alpha = alpha, value = 0, gradient = -1, slope = -1)
  }

  res <- make_schmidt_armijo_search(
    armijo_constant = 0.1,
    max_evaluations = 1
  )(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = 0.5,
    search_direction = 1
  )

  expect_equal(res$line_point$alpha, 0.5)
  expect_equal(res$line_point$value, 0)
  expect_equal(res$function_evaluations, 1)
  expect_equal(res$gradient_evaluations, 1)
})

test_that("Schmidt cubic Armijo backs off a nonfinite objective", {
  evaluated <- numeric()
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    if (alpha >= 1) {
      return(list(
        alpha = alpha,
        value = Inf,
        gradient = Inf,
        slope = Inf,
        parameters = alpha
      ))
    }
    list(alpha = alpha, value = 0, gradient = 0, slope = 0, parameters = alpha)
  }

  result <- make_schmidt_armijo_search(
    armijo_constant = 0.1,
    max_evaluations = 2
  )(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = 1
  )

  expect_equal(evaluated, c(1, 0.5))
  expect_equal(result$line_point$alpha, 0.5)
  expect_equal(result$function_evaluations, 2)
  expect_equal(result$gradient_evaluations, 2)
})

test_that("Schmidt cubic Armijo uses values when a trial gradient is nonfinite", {
  evaluated <- numeric()
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    if (alpha == 1) {
      return(list(
        alpha = alpha,
        value = 2,
        gradient = Inf,
        slope = Inf,
        parameters = alpha
      ))
    }
    list(alpha = alpha, value = 0, gradient = 0, slope = 0, parameters = alpha)
  }

  result <- make_schmidt_armijo_search(
    armijo_constant = 0.1,
    max_evaluations = 2
  )(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = 1
  )

  expect_equal(evaluated, c(1, 0.25))
  expect_equal(result$line_point$alpha, 0.25)
  expect_equal(result$function_evaluations, 2)
  expect_equal(result$gradient_evaluations, 2)
})

run_public_wolfe_step <- function(
  line_search,
  fg,
  initial_parameters,
  initial_alpha,
  armijo_constant,
  curvature_constant,
  strong_curvature,
  approximate_armijo
) {
  initial_gradient <- fg$gr(initial_parameters)
  search_direction <- -initial_gradient / norm2(initial_gradient)
  initial_point <- list(
    alpha = 0,
    value = fg$fn(initial_parameters),
    gradient = initial_gradient,
    slope = dot(initial_gradient, search_direction),
    parameters = initial_parameters
  )
  optimizer <- make_mize(
    method = "SD",
    line_search = line_search,
    norm_direction = TRUE,
    step0 = initial_alpha,
    c1 = armijo_constant,
    c2 = curvature_constant,
    strong_curvature = strong_curvature,
    approx_armijo = approximate_armijo,
    ls_max_fn = 100,
    ls_max_gr = 100,
    ls_max_fg = 200,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  optimizer <- mize_init(optimizer, initial_parameters, fg)
  result <- mize_step(optimizer, initial_parameters, fg)
  alpha <- as.numeric(
    (result$par - initial_parameters) / search_direction
  )
  trial_point <- list(
    alpha = alpha,
    value = result$f,
    gradient = result$g,
    slope = dot(result$g, search_direction),
    parameters = result$par
  )

  list(
    result = result,
    initial_point = initial_point,
    trial_point = trial_point
  )
}

test_that("public Wolfe controls reach every line-search backend", {
  approximate_only_fg <- list(
    fn = function(x) -1.25 * x^3 + 3.125 * x^2 - 2 * x + 1,
    gr = function(x) -3.75 * x^2 + 6.25 * x - 2
  )
  cases <- list(
    exact_strong = list(
      fg = condition_quadratic_fg(center = 1),
      initial_alpha = 0.5,
      c1 = 0.1,
      c2 = 0.5,
      strong_curvature = TRUE,
      approximate_armijo = FALSE,
      violates_strong = FALSE,
      violates_exact_armijo = FALSE
    ),
    weak_only = list(
      fg = condition_quadratic_fg(center = 10),
      initial_alpha = 12.5,
      c1 = 1e-4,
      c2 = 0.1,
      strong_curvature = FALSE,
      approximate_armijo = FALSE,
      violates_strong = TRUE,
      violates_exact_armijo = FALSE
    ),
    approximate_only = list(
      fg = approximate_only_fg,
      initial_alpha = 1,
      c1 = 0.1,
      c2 = 0.5,
      strong_curvature = TRUE,
      approximate_armijo = TRUE,
      violates_strong = FALSE,
      violates_exact_armijo = TRUE
    )
  )
  line_searches <- c(
    more_thuente = "more-thuente",
    rasmussen = "rasmussen",
    schmidt = "schmidt",
    hager_zhang = "hager-zhang"
  )

  for (line_search_name in names(line_searches)) {
    for (case_name in names(cases)) {
      case <- cases[[case_name]]
      info <- paste(line_search_name, case_name)
      step <- run_public_wolfe_step(
        line_search = line_searches[[line_search_name]],
        fg = case$fg,
        initial_parameters = 0,
        initial_alpha = case$initial_alpha,
        armijo_constant = case$c1,
        curvature_constant = case$c2,
        strong_curvature = case$strong_curvature,
        approximate_armijo = case$approximate_armijo
      )
      selected_conditions <- make_line_condition_policy(
        armijo_constant = case$c1,
        curvature_constant = case$c2,
        approximate_armijo = case$approximate_armijo,
        strong_curvature = case$strong_curvature
      )

      expect_equal(
        step$trial_point$alpha,
        case$initial_alpha,
        info = info
      )
      expect_true(
        selected_conditions$wolfe(step$initial_point, step$trial_point),
        info = info
      )
      if (case$violates_strong) {
        strong_conditions <- make_line_condition_policy(
          case$c1,
          case$c2,
          strong_curvature = TRUE
        )
        expect_false(
          strong_conditions$curvature(step$initial_point, step$trial_point),
          info = info
        )
      }
      if (case$violates_exact_armijo) {
        expect_false(
          selected_conditions$exact_armijo(
            step$initial_point,
            step$trial_point
          ),
          info = info
        )
      }
    }
  }
})

test_that("public Wolfe backends agree for separate and combined callbacks", {
  objective <- condition_quadratic_fg(center = 1)
  interfaces <- list(
    separate = objective,
    combined = c(
      objective,
      list(fg = function(x) list(fn = objective$fn(x), gr = objective$gr(x)))
    )
  )
  line_searches <- c("more-thuente", "rasmussen", "schmidt", "hager-zhang")

  for (line_search in line_searches) {
    steps <- lapply(interfaces, function(fg) {
      run_public_wolfe_step(
        line_search = line_search,
        fg = fg,
        initial_parameters = 0,
        initial_alpha = 0.5,
        armijo_constant = 0.1,
        curvature_constant = 0.5,
        strong_curvature = TRUE,
        approximate_armijo = FALSE
      )
    })

    expect_equal(
      steps$combined$trial_point,
      steps$separate$trial_point,
      info = line_search
    )
    expect_equal(
      c(fn = steps$combined$result$nf, gr = steps$combined$result$ng),
      c(fn = steps$separate$result$nf, gr = steps$separate$result$ng),
      info = line_search
    )
  }
})

test_that("Hager-Zhang initial exhaustion distinguishes acceptance from fallback", {
  quartic_start <- list(alpha = 0, value = 1, gradient = 4, slope = -16)
  quartic_trial <- list(alpha = 1, value = 81, gradient = -108, slope = 432)

  expect_true(line_point_satisfies_weak_curvature(
    quartic_start,
    quartic_trial,
    c2 = 0.5
  ))
  expect_false(line_point_satisfies_armijo(
    quartic_start,
    quartic_trial,
    c1 = 0.1
  ))
  weak_conditions <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = FALSE,
    approximation_tolerance = 1e-6
  )
  strong_conditions <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = TRUE
  )
  expect_false(strong_conditions$curvature(quartic_start, quartic_trial))
  expect_false(
    weak_conditions$wolfe(quartic_start, quartic_trial)
  )

  linear_start <- list(alpha = 0, value = 1, gradient = -1, slope = -1)
  linear_trial <- list(alpha = 1, value = 0, gradient = -1, slope = -1)

  expect_true(line_point_satisfies_armijo(linear_start, linear_trial, c1 = 0.1))
  expect_false(line_point_satisfies_weak_curvature(
    linear_start,
    linear_trial,
    c2 = 0.5
  ))
  expect_false(weak_conditions$wolfe(linear_start, linear_trial))
})

test_that("Hager-Zhang repairs a high-value negative-slope bound", {
  make_line_point <- function(alpha, value, slope) {
    list(
      alpha = alpha,
      value = value,
      slope = slope,
      gradient = slope,
      parameters = alpha
    )
  }

  initial_point <- make_line_point(alpha = 0, value = 0, slope = -1)
  bracket <- make_hager_zhang_bracket(
    lower_endpoint = make_line_point(alpha = 1, value = -1, slope = -1),
    upper_endpoint = make_line_point(alpha = 8, value = 2, slope = 1)
  )
  trial_point <- make_line_point(alpha = 8, value = 2, slope = -1)
  evaluate <- function(alpha, calc_gradient = TRUE) {
    if (alpha > 4) {
      make_line_point(alpha, value = 1, slope = -1)
    } else if (alpha < 3) {
      make_line_point(alpha, value = -1, slope = -1)
    } else {
      make_line_point(alpha, value = -0.5, slope = 1)
    }
  }
  evaluator <- make_line_evaluator(evaluate, max_evaluations = 5)
  condition_policy <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = TRUE
  )

  result <- update_hager_zhang_bracket(
    bracket = bracket,
    trial_point = trial_point,
    initial_point = initial_point,
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
    vapply(result$bracket, `[[`, numeric(1), "slope"),
    c(lower_endpoint = -1, upper_endpoint = 1)
  )
  expect_true(result$bracket$lower_endpoint$value <= initial_point$value)
})

test_that("Hager-Zhang retains bracket repair completed at budget exhaustion", {
  make_line_point <- function(alpha, value, slope) {
    list(
      alpha = alpha,
      value = value,
      slope = slope,
      gradient = slope,
      parameters = alpha
    )
  }
  initial_point <- make_line_point(alpha = 0, value = 1, slope = -1)
  bracket <- make_hager_zhang_bracket(
    lower_endpoint = make_line_point(alpha = 1, value = 1 + 5e-7, slope = -1),
    upper_endpoint = make_line_point(alpha = 5, value = 2, slope = 4)
  )
  trial_point <- make_line_point(alpha = 1.8, value = 2, slope = -1)
  evaluator <- make_line_evaluator(
    function(alpha, calc_gradient = TRUE) {
      expect_equal(alpha, 1.4)
      make_line_point(alpha, value = 0, slope = -1)
    },
    max_evaluations = 1
  )
  condition_policy <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.1,
    approximate_armijo = TRUE,
    strong_curvature = FALSE
  )
  result <- refine_hager_zhang_bracket_with_secants(
    bracket = bracket,
    trial_point = trial_point,
    initial_point = initial_point,
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
  make_line_point <- function(alpha, value, slope) {
    list(
      alpha = alpha,
      value = value,
      slope = slope,
      gradient = slope,
      parameters = alpha
    )
  }
  initial_point <- make_line_point(alpha = 0, value = 0, slope = -1)
  bracket <- make_hager_zhang_bracket(
    lower_endpoint = make_line_point(alpha = 1, value = -1, slope = -1),
    upper_endpoint = make_line_point(alpha = 8, value = 1, slope = 1)
  )
  cases <- list(
    nonfinite_trial = make_line_point(alpha = 4, value = Inf, slope = NaN),
    outside_bracket = make_line_point(alpha = 9, value = -1, slope = -1),
    upper_endpoint = make_line_point(alpha = 4, value = -0.5, slope = 0),
    lower_endpoint = make_line_point(alpha = 4, value = -0.5, slope = -1),
    needs_bisection = make_line_point(alpha = 4, value = 0.5, slope = -1)
  )

  for (classification in names(cases)) {
    expect_identical(
      classify_hager_zhang_trial(
        bracket,
        cases[[classification]],
        initial_point,
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
  weak_conditions <- make_line_condition_policy(
    c1,
    c2,
    strong_curvature = FALSE
  )

  wolfe_searches <- list(
    `more-thuente` = make_wolfe_line_search(
      more_thuente_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      method_policy = make_more_thuente_policy(
        alpha_max = Inf,
        safeguard_cubic = FALSE
      )
    ),
    rasmussen = make_rasmussen_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      interior_fraction = 0.1,
      expansion_factor = 3,
      relative_interval_tolerance = 1e-6
    ),
    schmidt = make_schmidt_wolfe_search(
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
      evaluate_line = setup$evaluate_line,
      initial_point = setup$initial_point,
      initial_alpha = alpha,
      search_direction = setup$pv
    )

    expect_equal(res$function_evaluations, 0, info = name)
    expect_equal(res$gradient_evaluations, 0, info = name)
    expect_equal(res$line_point$alpha, 0, info = name)
    expect_false(
      weak_conditions$wolfe(setup$initial_point, res$line_point),
      info = name
    )
  }

  armijo_res <- make_schmidt_armijo_search(
    armijo_constant = c1,
    max_evaluations = 0
  )(
    evaluate_line = setup$evaluate_line,
    initial_point = setup$initial_point,
    initial_alpha = alpha,
    search_direction = setup$pv
  )

  expect_equal(armijo_res$function_evaluations, 0)
  expect_equal(armijo_res$gradient_evaluations, 0)
  expect_equal(armijo_res$line_point$alpha, 0)
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

  evaluator <- make_line_evaluator(
    setup$evaluate_line,
    initial_point = setup$initial_point,
    max_evaluations = 4
  )
  res <- recover_finite_line_point(
    evaluator,
    alpha = 4,
    min_alpha = 0,
    max_evaluations = 4
  )

  expect_true(res$succeeded)
  expect_lte(res$line_point$alpha, 1)
  expect_true(wolfe_trial_point_is_usable(res$line_point))
})

test_that("finite-value recovery rejects a non-finite initial alpha", {
  phi <- function(alpha, calc_gradient = TRUE) {
    stop("phi should not be called")
  }

  evaluator <- make_line_evaluator(
    phi,
    initial_point = list(
      alpha = 0,
      value = 0,
      gradient = -1,
      slope = -1,
      parameters = 0
    ),
    max_evaluations = Inf
  )
  res <- recover_finite_line_point(
    evaluator,
    alpha = Inf,
    min_alpha = 0,
    max_evaluations = Inf
  )

  expect_null(res$line_point)
  expect_identical(environment(evaluator)$evaluation_count, 0L)
  expect_false(res$succeeded)
})

test_that("finite-value recovery stops without representable progress", {
  smallest_positive <- .Machine$double.xmin * .Machine$double.eps
  evaluated_alphas <- numeric()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated_alphas <<- c(evaluated_alphas, alpha)
    if (length(evaluated_alphas) > 1L) {
      stop("recovery repeated an endpoint")
    }
    list(
      alpha = alpha,
      value = Inf,
      gradient = Inf,
      slope = Inf,
      parameters = alpha
    )
  }

  evaluator <- make_line_evaluator(
    phi,
    initial_point = list(
      alpha = 0,
      value = 0,
      gradient = -1,
      slope = -1,
      parameters = 0
    ),
    max_evaluations = Inf
  )
  res <- recover_finite_line_point(
    evaluator,
    alpha = smallest_positive,
    min_alpha = 0,
    max_evaluations = Inf
  )

  expect_equal(evaluated_alphas, smallest_positive)
  expect_equal(res$line_point$alpha, smallest_positive)
  expect_identical(environment(evaluator)$evaluation_count, 1L)
  expect_false(res$succeeded)
})

test_that("weak Wolfe curvature includes the exact boundary", {
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  boundary_step <- list(
    alpha = 0.25,
    value = 0.8125,
    gradient = -0.5,
    slope = -0.5
  )
  conditions <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = FALSE
  )

  expect_true(weak_curvature_condition_is_met(-1, -0.5, 0.5))
  expect_true(conditions$wolfe(initial_point, boundary_step))

  evaluated_alphas <- numeric()
  evaluate <- function(alpha, calc_gradient = TRUE) {
    evaluated_alphas <<- c(evaluated_alphas, alpha)
    list(
      alpha = alpha,
      value = 1 - alpha + alpha^2,
      gradient = -1 + 2 * alpha,
      slope = -1 + 2 * alpha,
      parameters = alpha
    )
  }
  result <- make_hager_zhang_search(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    max_evaluations = 20,
    strong_curvature = FALSE,
    approximate_armijo = TRUE
  )(
    evaluate_line = evaluate,
    initial_point = initial_point,
    initial_alpha = 0.25,
    search_direction = 1
  )

  expect_equal(evaluated_alphas, 0.25)
  expect_equal(result$line_point$alpha, 0.25)
  expect_identical(result$function_evaluations, 1L)
  expect_identical(result$gradient_evaluations, 1L)
  expect_identical(result$termination_reason, "wolfe")
})
