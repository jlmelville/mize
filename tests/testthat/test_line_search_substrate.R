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
  expect_equal(recovered$function_evaluations, 2)
  expect_equal(evaluator_state$evaluation_count, 2)

  evaluator(1)
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

test_that("finite recovery requires a completely usable line point", {
  cases <- list(
    missing_gradient = list(gradient = NULL, parameters = 1),
    missing_parameters = list(gradient = 0, parameters = NULL),
    nonfinite_gradient = list(gradient = NaN, parameters = 1),
    nonfinite_parameters = list(gradient = 0, parameters = NaN)
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    recovered <- recover_finite_line_point(
      function(alpha, calc_gradient = TRUE) {
        list(
          alpha = alpha,
          value = 0,
          gradient = case$gradient,
          slope = 0,
          parameters = case$parameters
        )
      },
      alpha = 1,
      max_evaluations = 1
    )

    expect_false(recovered$succeeded, info = case_name)
    expect_identical(recovered$function_evaluations, 1L, info = case_name)
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
  expect_false(exact$armijo(initial_point, approximate_trial))
  expect_true(approximate$armijo(initial_point, approximate_trial))
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

  expect_error(
    make_line_condition_policy(0.6, 0.5),
    "curvature_constant"
  )
  expect_error(
    make_line_condition_policy(0.1, 0.5, approximate_armijo = NA),
    "approximate_armijo"
  )
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
  evaluator(2)
  failed_point <- evaluator(3)

  result <- finalize_wolfe_line_search_result(
    accepted_point = failed_point,
    evaluator = evaluator,
    termination_reason = "nonfinite_recovery"
  )

  expect_identical(
    names(result),
    c("line_point", "function_evaluations", "gradient_evaluations")
  )
  expect_equal(result$line_point$alpha, 2)
  expect_equal(result$function_evaluations, 2)
  expect_equal(result$gradient_evaluations, 2)
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

test_that("debug output includes the private line-search termination reason", {
  objective <- list(
    fn = function(parameters) (parameters - 1)^2,
    gr = function(parameters) 2 * (parameters - 1)
  )
  optimizer <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(),
        step_size = hager_zhang_ls(
          initial_step_length = 0.5,
          debug = TRUE
        )
      ),
      verbose = FALSE
    )
  )
  optimizer$count_res_fg <- FALSE

  messages <- character()
  withCallingHandlers(
    opt_loop(
      optimizer,
      par = 0,
      fg = objective,
      max_iter = 1,
      store_progress = FALSE,
      verbose = FALSE,
      grad_tol = NULL
    ),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl(
    "hager-zhang line search terminated: wolfe",
    messages,
    fixed = TRUE
  )))
})
