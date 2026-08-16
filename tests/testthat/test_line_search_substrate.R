test_that("line points have one validated canonical representation", {
  point <- new_line_point(
    alpha = 0.5,
    value = 1.25,
    gradient = c(1, 2),
    slope = -3,
    parameters = c(4, 5)
  )

  expect_identical(
    names(point),
    c("alpha", "value", "gradient", "slope", "parameters")
  )
  expect_type(point, "list")
  expect_error(new_line_point(c(1, 2), 1, 1, 1, 1), "alpha")
  expect_error(new_line_point(1, c(1, 2), 1, 1, 1), "value")
  expect_error(new_line_point(1, 1, "bad", 1, 1), "gradient")
  expect_error(new_line_point(1, 1, 1, c(1, 2), 1), "slope")
  expect_error(new_line_point(1, 1, 1, 1, "bad"), "parameters")
})

test_that("line evaluator owns recovery, accounting, and best decrease", {
  initial_point <- new_line_point(0, 4, -4, -4, 0)
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
  evaluator <- new_line_evaluator(phi, initial_point, max_evaluations = 3)

  recovered <- recover_finite_line_point(
    evaluator,
    alpha = 4,
    minimum_alpha = 0
  )
  expect_true(recovered$ok)
  expect_equal(recovered$point$alpha, 2)
  expect_equal(recovered$evaluations, 2)
  expect_equal(remaining_line_evaluations(evaluator), 1)

  evaluate_line_point(evaluator, 1)
  expect_equal(line_evaluator_evaluations(evaluator), 3)
  expect_equal(remaining_line_evaluations(evaluator), 0)
  expect_equal(line_evaluator_best_decrease(evaluator)$alpha, 1)
  expect_equal(evaluated, c(4, 2, 1))
  expect_error(evaluate_line_point(evaluator, 0.5), "budget exhausted")
})

test_that("condition policy reports named exact, approximate, and curvature results", {
  initial_point <- new_line_point(0, 1, -2, -2, 0)
  approximate_trial <- new_line_point(1, 0.9, 1, 1, 1)
  weak_only_trial <- new_line_point(1, 0.5, 2, 2, 1)

  approximate <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    approximate_armijo = TRUE,
    strong_curvature = TRUE,
    approximation_tolerance = 0.2
  )$classify(initial_point, approximate_trial)
  expect_identical(
    names(approximate),
    c(
      "exact_armijo",
      "approximate_armijo",
      "armijo",
      "weak_curvature",
      "strong_curvature",
      "curvature",
      "wolfe"
    )
  )
  expect_false(approximate$exact_armijo)
  expect_true(approximate$approximate_armijo)
  expect_true(approximate$armijo)
  expect_true(approximate$strong_curvature)
  expect_true(approximate$wolfe)

  weak <- new_line_condition_policy(
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    strong_curvature = FALSE
  )$classify(initial_point, weak_only_trial)
  expect_true(weak$weak_curvature)
  expect_false(weak$strong_curvature)
  expect_true(weak$curvature)
  expect_true(weak$wolfe)

  expect_error(
    new_line_condition_policy(0.6, 0.5),
    "curvature_constant"
  )
  expect_error(
    new_line_condition_policy(0.1, 0.5, approximate_armijo = NA),
    "approximate_armijo"
  )
})

test_that("line result finalizer separates acceptance, fallback, and failure", {
  initial_point <- new_line_point(0, 1, -1, -1, 0)
  policy <- new_line_condition_policy(0.1, 0.5)
  phi <- function(alpha, calc_gradient = TRUE) {
    points <- list(
      `1` = list(alpha = 1, f = 0.4, df = -1, d = -1, par = 1),
      `2` = list(alpha = 2, f = 0.2, df = -2, d = -2, par = 2),
      `3` = list(alpha = 3, f = Inf, df = Inf, d = Inf, par = 3)
    )
    points[[as.character(alpha)]]
  }
  evaluator <- new_line_evaluator(phi, initial_point, max_evaluations = 3)
  evaluate_line_point(evaluator, 1)
  evaluate_line_point(evaluator, 2)
  nonfinite <- evaluate_line_point(evaluator, 3)

  fallback <- finalize_line_search_result(
    candidate = nonfinite,
    initial_point = initial_point,
    evaluator = evaluator,
    condition_policy = policy,
    termination_reason = "nonfinite_recovery"
  )
  expect_false(fallback$accepted)
  expect_identical(fallback$selection, "strict_decrease")
  expect_identical(fallback$termination_reason, "nonfinite_recovery")
  expect_equal(fallback$point$alpha, 2)
  expect_equal(fallback$evaluations, 3)
  expect_true(fallback$recovered_nonfinite)

  accepted_evaluator <- new_line_evaluator(
    function(alpha, calc_gradient = TRUE) {
      list(alpha = alpha, f = 0, df = 0, d = 0, par = alpha)
    },
    initial_point,
    max_evaluations = 1,
    initial_step = as_legacy_line_step(initial_point)
  )
  accepted_point <- evaluate_line_point(accepted_evaluator, 1)
  accepted <- finalize_line_search_result(
    candidate = accepted_point,
    initial_point = initial_point,
    evaluator = accepted_evaluator,
    condition_policy = policy,
    termination_reason = "wolfe"
  )
  expect_true(accepted$accepted)
  expect_identical(accepted$selection, "wolfe")
  expect_identical(accepted$termination_reason, "wolfe")
  expect_true(accepted$conditions$wolfe)

  unchanged <- finalize_line_search_result(
    candidate = NULL,
    initial_point = initial_point,
    evaluator = new_line_evaluator(phi, initial_point, max_evaluations = 0),
    condition_policy = policy,
    termination_reason = "budget_exhausted"
  )
  expect_false(unchanged$accepted)
  expect_identical(unchanged$selection, "unchanged_start")
  expect_identical(unchanged$termination_reason, "budget_exhausted")
  expect_equal(unchanged$point, initial_point)

  progress <- finalize_line_search_result(
    candidate = NULL,
    initial_point = initial_point,
    evaluator = new_line_evaluator(
      phi,
      initial_point,
      max_evaluations = 1,
      initial_step = as_legacy_line_step(initial_point)
    ),
    condition_policy = policy,
    termination_reason = "progress_failure"
  )
  expect_identical(progress$selection, "unchanged_start")
  expect_identical(progress$termination_reason, "progress_failure")
})

test_that("production finalization retains safe selection without diagnostics", {
  initial_point <- new_line_point(0, 1, -1, -1, 0)
  initial_step <- as_legacy_line_step(initial_point)
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
    initial_point,
    max_evaluations = 2,
    initial_step = initial_step
  )
  evaluate_legacy_line_step(evaluator, 2)
  failed_candidate <- evaluate_legacy_line_step(evaluator, 3)

  result <- finalize_legacy_line_search_result(
    candidate = failed_candidate,
    evaluator = evaluator,
    termination_reason = "nonfinite_recovery"
  )

  expect_identical(names(result), c("step", "nfn", "ngr"))
  expect_equal(result$step$alpha, 2)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})

test_that("shared Wolfe core shape retains the public private-boundary result", {
  observed <- NULL
  core <- function(
    evaluator,
    initial_point,
    initial_alpha,
    condition_policy,
    direction,
    method_policy
  ) {
    observed <<- list(
      initial_point = initial_point,
      initial_alpha = initial_alpha,
      direction = direction,
      method_policy = method_policy
    )
    list(
      candidate = evaluate_line_point(evaluator, initial_alpha),
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
  expect_equal(observed$initial_alpha, 1)
  expect_equal(observed$direction, 3)
  expect_identical(observed$method_policy, list(method = "test"))
})

test_that("supported Wolfe cores share one private signature", {
  expected <- c(
    "evaluator",
    "initial_point",
    "initial_alpha",
    "condition_policy",
    "direction",
    "method_policy"
  )

  expect_identical(names(formals(more_thuente_core)), expected)
  expect_identical(names(formals(rasmussen_core)), expected)
  expect_identical(names(formals(schmidt_core)), expected)
})
