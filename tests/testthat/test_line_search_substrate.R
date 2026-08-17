test_that("line evaluator validates, recovers, accounts, and tracks decrease", {
  initial_step <- list(alpha = 0, f = 4, df = -4, d = -4, par = 0)
  evaluated <- numeric()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    list(
      alpha = alpha,
      f = if (alpha > 2) Inf else (alpha - 1)^2,
      df = 2 * (alpha - 1),
      d = 2 * (alpha - 1),
      par = alpha
    )
  }
  evaluator <- new_line_evaluator(
    phi,
    initial_step,
    max_evaluations = 3
  )
  evaluator_state <- environment(evaluator)

  recovered <- find_finite(
    evaluator,
    alpha = 4,
    min_alpha = 0,
    max_fn = Inf
  )
  expect_true(recovered$ok)
  expect_equal(recovered$step$alpha, 2)
  expect_equal(recovered$nfn, 2)
  expect_equal(evaluator_state$evaluation_count, 2)

  evaluator(1)
  expect_equal(evaluator_state$evaluation_count, 3)
  expect_equal(evaluator_state$best_decrease$alpha, 1)
  expect_equal(evaluated, c(4, 2, 1))
  expect_error(evaluator(0.5), "budget exhausted")
})

test_that("line evaluator rejects malformed callback results after accounting", {
  valid_initial_step <- list(
    alpha = 0,
    f = 1,
    df = -1,
    d = -1,
    par = 0
  )
  malformed_steps <- list(
    not_a_list = 1,
    missing_slope = list(alpha = 1, f = 0),
    vector_value = list(alpha = 1, f = c(0, 1), d = 0),
    character_gradient = list(alpha = 1, f = 0, df = "bad", d = 0),
    character_parameters = list(alpha = 1, f = 0, d = 0, par = "bad")
  )

  for (name in names(malformed_steps)) {
    evaluator <- new_line_evaluator(
      function(alpha, calc_gradient = TRUE) malformed_steps[[name]],
      valid_initial_step,
      max_evaluations = 1
    )

    expect_error(
      evaluator(1),
      "line step",
      info = name
    )
    expect_equal(environment(evaluator)$evaluation_count, 1, info = name)
  }
})

test_that("condition policy selects exact, approximate, and curvature rules", {
  initial_step <- list(alpha = 0, f = 1, df = -2, d = -2, par = 0)
  approximate_trial <- list(alpha = 1, f = 0.9, df = 1, d = 1, par = 1)
  weak_only_trial <- list(alpha = 1, f = 0.5, df = 2, d = 2, par = 1)

  exact <- new_line_condition_policy(0.1, 0.5)
  approximate <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = TRUE,
    approximation_tolerance = 0.2
  )
  expect_false(exact$armijo(initial_step, approximate_trial))
  expect_true(approximate$armijo(initial_step, approximate_trial))
  expect_true(approximate$curvature(initial_step, approximate_trial))
  expect_true(approximate$wolfe(initial_step, approximate_trial))

  weak <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    strong_curvature = FALSE
  )
  expect_true(weak$curvature(initial_step, weak_only_trial))
  expect_false(exact$curvature(initial_step, weak_only_trial))
  expect_true(weak$wolfe(initial_step, weak_only_trial))

  expect_error(
    new_line_condition_policy(0.6, 0.5),
    "curvature_constant"
  )
  expect_error(
    new_line_condition_policy(0.1, 0.5, approximate_armijo = NA),
    "approximate_armijo"
  )
})

test_that("Wolfe finalization keeps the best usable evaluated decrease", {
  initial_step <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)
  phi <- function(alpha, calc_gradient = TRUE) {
    list(
      alpha = alpha,
      f = if (alpha == 2) 0.25 else Inf,
      df = if (alpha == 2) -0.5 else Inf,
      d = if (alpha == 2) -0.5 else Inf,
      par = alpha
    )
  }
  evaluator <- new_line_evaluator(
    phi,
    initial_step,
    max_evaluations = 2
  )
  evaluator(2)
  failed_candidate <- evaluator(3)

  result <- finalize_wolfe_line_search_result(
    candidate = failed_candidate,
    evaluator = evaluator,
    termination_reason = "nonfinite_recovery"
  )

  expect_identical(names(result), c("step", "nfn", "ngr"))
  expect_equal(result$step$alpha, 2)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})

test_that("shared Wolfe core boundary preserves result and callback types", {
  observed <- NULL
  core <- function(
    evaluator,
    initial_step,
    initial_alpha,
    condition_policy,
    search_direction,
    method_policy
  ) {
    observed <<- list(
      initial_step = initial_step,
      initial_alpha = initial_alpha,
      search_direction = search_direction,
      method_policy = method_policy
    )
    list(
      candidate = evaluator(initial_alpha),
      termination_reason = "wolfe"
    )
  }
  search <- new_wolfe_line_search(
    core = core,
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    max_evaluations = 2,
    method_policy = list(method = "test")
  )
  phi <- function(alpha, calc_gradient = TRUE) {
    list(alpha = alpha, f = 0, df = 0, d = 0, par = alpha)
  }
  result <- search(
    phi = phi,
    step0 = list(alpha = 0, f = 1, df = -1, d = -1, par = 0),
    alpha = 1,
    total_max_fn = 1,
    total_max_gr = 1,
    total_max_fg = 2,
    pm = 3
  )

  expect_identical(names(result), c("step", "nfn", "ngr"))
  expect_identical(names(result$step), c("alpha", "f", "df", "d", "par"))
  expect_equal(result$step$alpha, 1)
  expect_equal(result$nfn, 1)
  expect_equal(result$ngr, 1)
  expect_equal(observed$initial_step$alpha, 0)
  expect_equal(observed$initial_alpha, 1)
  expect_equal(observed$search_direction, 3)
  expect_identical(observed$method_policy, list(method = "test"))
})

test_that("supported Wolfe cores share one private signature", {
  expected <- c(
    "evaluator",
    "initial_step",
    "initial_alpha",
    "condition_policy",
    "search_direction",
    "method_policy"
  )

  expect_identical(names(formals(more_thuente_core)), expected)
  expect_identical(names(formals(run_bracket_zoom)), expected)
})
