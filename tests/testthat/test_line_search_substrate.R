test_that("line evaluator validates, recovers, accounts, and tracks decrease", {
  initial_point <- list(
    alpha = 0,
    value = 4,
    gradient = -4,
    slope = -4,
    parameters = 0
  )
  evaluated <- numeric()
  evaluate_line <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    list(
      alpha = alpha,
      value = if (alpha > 2) Inf else (alpha - 1)^2,
      gradient = 2 * (alpha - 1),
      slope = 2 * (alpha - 1),
      parameters = alpha
    )
  }
  evaluator <- make_line_evaluator(
    evaluate_line,
    initial_point,
    max_evaluations = 3
  )
  evaluator_state <- environment(evaluator)

  recovered <- recover_finite_line_point(
    evaluator,
    alpha = 4,
    min_alpha = 0,
    max_evaluations = Inf
  )
  expect_true(recovered$succeeded)
  expect_equal(recovered$line_point$alpha, 2)
  expect_equal(evaluator_state$evaluation_count, 2)

  recover_finite_line_point(evaluator, alpha = 1, max_evaluations = 1)
  expect_equal(evaluator_state$evaluation_count, 3)
  expect_equal(evaluator_state$best_decreasing_point$alpha, 1)
  expect_equal(evaluated, c(4, 2, 1))
  expect_error(evaluator(0.5), "budget exhausted")
})

test_that("line evaluator rejects malformed callback results after accounting", {
  valid_initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  malformed_points <- list(
    not_a_list = 1,
    missing_slope = list(alpha = 1, value = 0),
    missing_gradient = list(alpha = 1, value = 0, slope = 0, parameters = 1),
    missing_parameters = list(alpha = 1, value = 0, gradient = 0, slope = 0),
    vector_value = list(alpha = 1, value = c(0, 1), slope = 0),
    character_gradient = list(
      alpha = 1,
      value = 0,
      gradient = "bad",
      slope = 0,
      parameters = 1
    ),
    character_parameters = list(
      alpha = 1,
      value = 0,
      gradient = 0,
      slope = 0,
      parameters = "bad"
    )
  )

  for (name in names(malformed_points)) {
    evaluator <- make_line_evaluator(
      function(alpha, calc_gradient = TRUE) malformed_points[[name]],
      valid_initial_point,
      max_evaluations = 1
    )

    expect_error(
      evaluator(1),
      "point",
      info = name
    )
    expect_equal(environment(evaluator)$evaluation_count, 1, info = name)
  }
})

test_that("finite recovery rejects nonfinite numeric line points", {
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  cases <- list(
    nonfinite_gradient = list(gradient = NaN, parameters = 1),
    nonfinite_parameters = list(gradient = 0, parameters = NaN)
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    evaluator <- make_line_evaluator(
      function(alpha, calc_gradient = TRUE) {
        list(
          alpha = alpha,
          value = 0,
          gradient = case$gradient,
          slope = 0,
          parameters = case$parameters
        )
      },
      initial_point = initial_point,
      max_evaluations = 1
    )
    recovered <- recover_finite_line_point(
      evaluator,
      alpha = 1,
      max_evaluations = 1
    )

    expect_false(recovered$succeeded, info = case_name)
    expect_identical(
      environment(evaluator)$evaluation_count,
      1L,
      info = case_name
    )
  }
})

test_that("Wolfe searches reject incomplete trial points", {
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  incomplete_points <- list(
    missing_gradient = list(alpha = 1, value = 0, slope = 0, parameters = 1),
    missing_parameters = list(alpha = 1, value = 0, gradient = 0, slope = 0)
  )

  for (name in names(incomplete_points)) {
    search <- make_rasmussen_wolfe_search(
      armijo_constant = 0.1,
      curvature_constant = 0.5,
      max_evaluations = 1
    )
    expect_error(
      search(
        evaluate_line = function(alpha, calc_gradient = TRUE) {
          incomplete_points[[name]]
        },
        initial_point = initial_point,
        initial_alpha = 1,
        search_direction = 1
      ),
      "Wolfe trial point",
      info = name
    )
  }
})

test_that("shared wrapper cannot report Wolfe for an incomplete core point", {
  search <- make_wolfe_line_search(
    core = function(
      evaluator,
      initial_point,
      initial_alpha,
      condition_policy,
      search_direction,
      method_policy
    ) {
      make_line_search_core_result(
        "wolfe",
        list(alpha = initial_alpha, value = 0, slope = 0)
      )
    },
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    max_evaluations = 1
  )
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )

  result <- search(
    evaluate_line = function(alpha, calc_gradient = TRUE) {
      stop("line callback should not be called")
    },
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = 1
  )

  expect_identical(result$termination_reason, "invalid_wolfe_point")
  expect_identical(result$line_point, initial_point)
  expect_identical(result$function_evaluations, 0L)
  expect_identical(result$gradient_evaluations, 0L)
})

test_that("condition policy selects exact, approximate, and curvature rules", {
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -2,
    slope = -2,
    parameters = 0
  )
  approximate_trial <- list(
    alpha = 1,
    value = 0.9,
    gradient = 1,
    slope = 1,
    parameters = 1
  )
  weak_only_trial <- list(
    alpha = 1,
    value = 0.5,
    gradient = 2,
    slope = 2,
    parameters = 1
  )

  exact <- make_line_condition_policy(0.1, 0.5)
  approximate <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = TRUE,
    approximation_tolerance = 0.2
  )
  expect_false(exact$selected_armijo(initial_point, approximate_trial))
  expect_true(approximate$selected_armijo(initial_point, approximate_trial))
  expect_false(approximate$exact_armijo(initial_point, approximate_trial))
  expect_true(approximate$curvature(initial_point, approximate_trial))
  expect_true(approximate$wolfe(initial_point, approximate_trial))

  weak <- make_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    strong_curvature = FALSE
  )
  expect_true(weak$curvature(initial_point, weak_only_trial))
  expect_false(exact$curvature(initial_point, weak_only_trial))
  expect_true(weak$wolfe(initial_point, weak_only_trial))
})

test_that("Wolfe finalization keeps the best usable evaluated decrease", {
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  evaluate_line <- function(alpha, calc_gradient = TRUE) {
    list(
      alpha = alpha,
      value = if (alpha == 2) 0.25 else Inf,
      gradient = if (alpha == 2) -0.5 else Inf,
      slope = if (alpha == 2) -0.5 else Inf,
      parameters = alpha
    )
  }
  evaluator <- make_line_evaluator(
    evaluate_line,
    initial_point,
    max_evaluations = 2
  )
  recover_finite_line_point(evaluator, alpha = 2, max_evaluations = 1)
  failed_point <- recover_finite_line_point(
    evaluator,
    alpha = 3,
    max_evaluations = 1
  )$line_point

  result <- finalize_wolfe_line_search_result(
    accepted_point = failed_point,
    evaluator = evaluator,
    termination_reason = "nonfinite_recovery"
  )

  expect_identical(
    names(result),
    c(
      "line_point",
      "function_evaluations",
      "gradient_evaluations",
      "outcome"
    )
  )
  expect_equal(result$line_point$alpha, 2)
  expect_equal(result$function_evaluations, 2)
  expect_equal(result$gradient_evaluations, 2)
  expect_identical(result$outcome, "improving_fallback")
})

test_that("shared Wolfe core boundary preserves result and callback types", {
  observed <- NULL
  core <- function(
    evaluator,
    initial_point,
    initial_alpha,
    condition_policy,
    search_direction,
    method_policy
  ) {
    observed <<- list(
      initial_point = initial_point,
      initial_alpha = initial_alpha,
      search_direction = search_direction,
      method_policy = method_policy
    )
    make_line_search_core_result("wolfe", evaluator(initial_alpha))
  }
  search <- make_wolfe_line_search(
    core = core,
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    max_evaluations = 2,
    method_policy = list(method = "test")
  )
  evaluate_line <- function(alpha, calc_gradient = TRUE) {
    list(alpha = alpha, value = 0, gradient = 0, slope = 0, parameters = alpha)
  }
  result <- search(
    evaluate_line = evaluate_line,
    initial_point = list(
      alpha = 0,
      value = 1,
      gradient = -1,
      slope = -1,
      parameters = 0
    ),
    initial_alpha = 1,
    remaining_function_evaluations = 1,
    remaining_gradient_evaluations = 1,
    remaining_combined_evaluations = 2,
    search_direction = 3
  )

  expect_identical(
    names(result),
    c(
      "line_point",
      "function_evaluations",
      "gradient_evaluations",
      "outcome",
      "termination_reason"
    )
  )
  expect_identical(
    names(result$line_point),
    c("alpha", "value", "gradient", "slope", "parameters")
  )
  expect_equal(result$line_point$alpha, 1)
  expect_equal(result$function_evaluations, 1)
  expect_equal(result$gradient_evaluations, 1)
  expect_identical(result$outcome, "wolfe")
  expect_equal(observed$initial_point$alpha, 0)
  expect_equal(observed$initial_alpha, 1)
  expect_equal(observed$search_direction, 3)
  expect_identical(observed$method_policy, list(method = "test"))
  expect_identical(result$termination_reason, "wolfe")
})

test_that("line-search backends share one explicit callable protocol", {
  expected <- c(
    "evaluate_line",
    "initial_point",
    "initial_alpha",
    "remaining_function_evaluations",
    "remaining_gradient_evaluations",
    "remaining_combined_evaluations",
    "search_direction"
  )
  searches <- list(
    more_thuente = make_wolfe_line_search(
      core = more_thuente_core,
      armijo_constant = 0.05,
      curvature_constant = 0.1,
      method_policy = make_more_thuente_policy()
    ),
    rasmussen = make_rasmussen_wolfe_search(0.05, 0.1),
    schmidt = make_schmidt_wolfe_search(0.05, 0.1),
    hager_zhang = make_hager_zhang_search(0.05, 0.1),
    schmidt_armijo = make_schmidt_armijo_search()
  )

  for (name in names(searches)) {
    expect_identical(names(formals(searches[[name]])), expected, info = name)
  }
})

test_that("line parameter projection uses exact numeric identity", {
  unnamed <- c(1, 0)
  named <- c(first = 1, second = -0)

  expect_equal(project_line_parameters(c(1, 2), 0.5, c(-2, 4)), c(0, 4))
  expect_true(line_parameters_are_identical(unnamed, named))
  expect_true(line_parameters_are_identical(c(0, -0), c(-0, 0)))
  expect_false(line_parameters_are_identical(c(1, 2), c(1, 2, 3)))
  expect_false(line_parameters_are_identical(c(1, 2), c(1, Inf)))
})

test_that("bracket admission reconstructs a truthful colliding point", {
  direction <- -2^-53
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -0.5 / direction,
    slope = -0.5,
    parameters = 1
  )
  lower_point <- list(
    alpha = 1,
    value = 0.9,
    gradient = -0.4 / direction,
    slope = -0.4,
    parameters = project_line_parameters(1, 1, direction)
  )
  upper_point <- list(
    alpha = 2,
    value = 0.8,
    gradient = 0,
    slope = 0,
    parameters = project_line_parameters(1, 2, direction)
  )
  conditions <- make_line_condition_policy(0.25, 0.5)

  expect_false(conditions$wolfe(initial_point, upper_point))
  admission <- admit_bracketed_line_alpha(
    alpha = 1.5,
    first_endpoint = lower_point,
    second_endpoint = upper_point,
    initial_parameters = 1,
    search_direction = direction
  )

  expect_identical(admission$status, "collision")
  expect_identical(admission$matching_endpoint, "second")
  expect_equal(admission$condition_point$alpha, 1.5)
  expect_equal(admission$condition_point$parameters, upper_point$parameters)
  expect_equal(admission$condition_point$value, upper_point$value)
  expect_equal(admission$condition_point$gradient, upper_point$gradient)
  expect_equal(admission$condition_point$slope, upper_point$slope)
  expect_true(conditions$wolfe(initial_point, admission$condition_point))
  expect_null(admission$replacement)
  expect_null(find_novel_bracketed_line_parameters(
    lower_point,
    upper_point,
    1,
    direction
  ))
})

test_that("bracket admission finds a third realized parameter vector", {
  direction <- -2^-52
  lower_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1 / direction,
    slope = -1,
    parameters = 1
  )
  upper_point <- list(
    alpha = 2,
    value = 2,
    gradient = 1 / direction,
    slope = 1,
    parameters = project_line_parameters(1, 2, direction)
  )

  admission <- admit_bracketed_line_alpha(
    alpha = 0.25,
    first_endpoint = lower_point,
    second_endpoint = upper_point,
    initial_parameters = 1,
    search_direction = direction
  )
  replacement <- find_novel_bracketed_line_parameters(
    lower_point,
    upper_point,
    1,
    direction
  )

  expect_identical(admission$status, "collision")
  expect_equal(replacement$alpha, 1)
  expect_false(line_parameters_are_identical(
    replacement$parameters,
    lower_point$parameters
  ))
  expect_false(line_parameters_are_identical(
    replacement$parameters,
    upper_point$parameters
  ))
})

test_that("a novel parameter vector remains evaluable at an equal objective", {
  initial_parameters <- 1
  search_direction <- -2^-52
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1 / search_direction,
    slope = -1,
    parameters = initial_parameters
  )
  upper_point <- list(
    alpha = 2,
    value = 1,
    gradient = 1 / search_direction,
    slope = 1,
    parameters = project_line_parameters(
      initial_parameters,
      2,
      search_direction
    )
  )
  callback_parameters <- numeric()
  evaluator <- make_line_evaluator(
    function(alpha, calc_gradient = TRUE) {
      parameters <- project_line_parameters(
        initial_parameters,
        alpha,
        search_direction
      )
      callback_parameters <<- c(callback_parameters, parameters)
      list(
        alpha = alpha,
        value = 1,
        gradient = 0,
        slope = 0,
        parameters = parameters
      )
    },
    initial_point = initial_point,
    max_evaluations = 1
  )

  recovery <- recover_bracketed_line_point(
    evaluator = evaluator,
    alpha = 1,
    min_alpha = 0,
    first_endpoint = initial_point,
    second_endpoint = upper_point,
    initial_point = initial_point,
    condition_policy = make_line_condition_policy(0.1, 0.5),
    search_direction = search_direction,
    max_evaluations = 1
  )

  expect_true(recovery$succeeded)
  expect_false(recovery$accepted)
  expect_identical(environment(evaluator)$evaluation_count, 1L)
  expect_equal(callback_parameters, recovery$line_point$parameters)
  expect_equal(recovery$line_point$value, initial_point$value)
})

test_that("bracketed nonfinite recovery neither repeats nor changes its reason", {
  initial_parameters <- 1
  search_direction <- -2^-52
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1 / search_direction,
    slope = -1,
    parameters = initial_parameters
  )
  upper_point <- list(
    alpha = 4,
    value = 2,
    gradient = 1 / search_direction,
    slope = 1,
    parameters = project_line_parameters(
      initial_parameters,
      4,
      search_direction
    )
  )
  callback_parameters <- numeric()
  evaluator <- make_line_evaluator(
    function(alpha, calc_gradient = TRUE) {
      parameters <- project_line_parameters(
        initial_parameters,
        alpha,
        search_direction
      )
      callback_parameters <<- c(callback_parameters, parameters)
      list(
        alpha = alpha,
        value = Inf,
        gradient = Inf,
        slope = Inf,
        parameters = parameters
      )
    },
    initial_point = initial_point,
    max_evaluations = Inf
  )

  recovery <- recover_bracketed_line_point(
    evaluator = evaluator,
    alpha = 2,
    min_alpha = 0,
    first_endpoint = initial_point,
    second_endpoint = upper_point,
    initial_point = initial_point,
    condition_policy = make_line_condition_policy(0.1, 0.5),
    search_direction = search_direction,
    max_evaluations = Inf
  )

  expect_false(recovery$succeeded)
  expect_null(recovery$termination_reason)
  expect_gt(length(callback_parameters), 0L)
  expect_identical(
    length(unique(callback_parameters)),
    length(callback_parameters)
  )
})

test_that("bracket admission preserves legacy behavior for inconsistent maps", {
  initial_parameters <- 1
  direction <- -2^-53
  lower_point <- list(
    alpha = 0,
    value = 1,
    gradient = 1,
    slope = -1,
    parameters = initial_parameters
  )
  inconsistent_upper <- list(
    alpha = 2,
    value = 2,
    gradient = 1,
    slope = 1,
    parameters = 123
  )

  cases <- list(
    inconsistent_endpoint = list(
      parameters = initial_parameters,
      direction = direction,
      endpoint = inconsistent_upper
    ),
    missing_direction = list(
      parameters = initial_parameters,
      direction = NULL,
      endpoint = lower_point
    ),
    nonfinite_parameters = list(
      parameters = Inf,
      direction = direction,
      endpoint = inconsistent_upper
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    admission <- admit_bracketed_line_alpha(
      alpha = 1,
      first_endpoint = lower_point,
      second_endpoint = case$endpoint,
      initial_parameters = case$parameters,
      search_direction = case$direction
    )
    expect_identical(admission$status, "unchecked", info = case_name)
  }
})

test_that("the fast admission path falls back for nonconformable endpoints", {
  initial_parameters <- c(1, 2)
  search_direction <- c(-0.25, 0)
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = c(1, 0),
    slope = -0.25,
    parameters = initial_parameters
  )
  inconsistent_endpoint <- list(
    alpha = 2,
    value = 2,
    gradient = 0,
    slope = 0,
    parameters = 0.5
  )
  evaluator <- make_line_evaluator(
    function(alpha, calc_gradient = TRUE) {
      list(
        alpha = alpha,
        value = 0.5,
        gradient = c(0, 0),
        slope = 0,
        parameters = initial_parameters + alpha * search_direction
      )
    },
    max_evaluations = 1
  )
  initialize_line_evaluator_parameter_witness(
    evaluator,
    initial_point,
    search_direction
  )

  recovery <- recover_bracketed_line_point(
    evaluator = evaluator,
    alpha = 1,
    min_alpha = 0,
    first_endpoint = initial_point,
    second_endpoint = inconsistent_endpoint,
    initial_point = initial_point,
    condition_policy = make_line_condition_policy(0.1, 0.5),
    search_direction = search_direction,
    max_evaluations = 1
  )

  expect_true(recovery$succeeded)
  expect_false(recovery$accepted)
  expect_equal(recovery$line_point$alpha, 1)
  expect_identical(environment(evaluator)$evaluation_count, 1L)
})
