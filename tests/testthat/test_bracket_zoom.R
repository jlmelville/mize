test_that("Schmidt brackets when the expanded objective stops decreasing", {
  initial_point <- list(alpha = 0, value = 1, gradient = -1, slope = -1)
  previous_point <- list(alpha = 1, value = 0.5, gradient = -1, slope = -1)
  trial_point <- list(alpha = 2, value = 0.5, gradient = 0, slope = 0)
  expansion_state <- list(
    initial_point = initial_point,
    previous_point = previous_point,
    iteration = 1L
  )
  conditions <- make_line_condition_policy(0.05, 0.1)

  result <- classify_schmidt_expansion(
    expansion_state,
    trial_point,
    conditions
  )

  expect_true(conditions$wolfe(initial_point, trial_point))
  expect_false(result$accepted)
  expect_equal(result$bracket, list(previous_point, trial_point))
})

test_that("bracket-and-zoom rejects non-descent directions before callbacks", {
  searches <- list(
    rasmussen = make_rasmussen_wolfe_search(0.05, 0.1),
    schmidt = make_schmidt_wolfe_search(0.05, 0.1)
  )
  phi <- function(alpha, calc_gradient = TRUE) {
    stop("phi should not be called")
  }

  for (name in names(searches)) {
    search <- searches[[name]]
    non_descent_step <- list(
      alpha = 0,
      value = 1,
      gradient = 1,
      slope = 1,
      parameters = 0
    )
    result <- search(
      phi,
      non_descent_step,
      initial_alpha = 1,
      search_direction = 1
    )
    expect_equal(result$line_point, non_descent_step, info = name)
    expect_identical(result$function_evaluations, 0L, info = name)
    expect_identical(result$gradient_evaluations, 0L, info = name)
  }
})

test_that("bracket-and-zoom expansion never evaluates an infinite proposal", {
  searches <- list(
    rasmussen = make_rasmussen_wolfe_search(0.05, 0.1),
    schmidt = make_schmidt_wolfe_search(0.05, 0.1)
  )
  initial_point <- list(
    alpha = 0,
    value = 0,
    gradient = -1,
    slope = -1,
    parameters = 0
  )

  for (name in names(searches)) {
    evaluated_alphas <- numeric()
    phi <- function(alpha, calc_gradient = TRUE) {
      if (!is.finite(alpha)) {
        stop("non-finite alpha reached the callback")
      }
      evaluated_alphas <<- c(evaluated_alphas, alpha)
      list(
        alpha = alpha,
        value = -alpha,
        gradient = -1,
        slope = -1,
        parameters = alpha
      )
    }

    result <- searches[[name]](
      phi,
      initial_point = initial_point,
      initial_alpha = 1e308,
      search_direction = 1
    )

    expect_equal(
      evaluated_alphas,
      c(1e308, .Machine$double.xmax),
      info = name
    )
    expect_equal(result$line_point$alpha, .Machine$double.xmax, info = name)
    expect_identical(result$function_evaluations, 2L, info = name)
    expect_identical(result$gradient_evaluations, 2L, info = name)
  }
})

test_that("Rasmussen zoom stops when no interior alpha is representable", {
  smallest_positive <- .Machine$double.xmin * .Machine$double.eps
  evaluated_alphas <- numeric()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated_alphas <<- c(evaluated_alphas, alpha)
    if (length(evaluated_alphas) > 1L) {
      stop("zoom repeated without representable progress")
    }
    list(alpha = alpha, value = 2, gradient = 1, slope = 1, parameters = alpha)
  }
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )

  result <- make_rasmussen_wolfe_search(0.05, 0.1)(
    phi,
    initial_point = initial_point,
    initial_alpha = smallest_positive,
    search_direction = 1
  )

  expect_equal(evaluated_alphas, smallest_positive)
  expect_equal(result$line_point, initial_point)
  expect_identical(result$function_evaluations, 1L)
  expect_identical(result$gradient_evaluations, 1L)
})

test_that("Schmidt zoom proposals stay safeguarded inside their bracket", {
  zoom_state <- initialize_schmidt_zoom_state(list(
    list(alpha = 0, value = 1, gradient = -1, slope = -1),
    list(alpha = 1, value = 2, gradient = 1, slope = 1)
  ))
  zoom_state$previous_proposal_not_safely_interior <- TRUE

  proposal <- propose_schmidt_zoom_alpha(
    zoom_state,
    interior_fraction = 0.1
  )

  expect_gte(proposal$alpha, 0.1)
  expect_lte(proposal$alpha, 0.9)
})

test_that("Schmidt zoom updates preserve the lower-value role", {
  conditions <- make_line_condition_policy(0.1, 0.5)
  initial_point <- list(alpha = 0, value = 1, gradient = -1, slope = -1)
  trial_point <- list(alpha = 1, value = 0.5, gradient = -1, slope = -1)
  other_point <- list(alpha = 2, value = 2, gradient = 1, slope = 1)
  zoom_state <- initialize_schmidt_zoom_state(list(other_point, initial_point))
  expect_equal(zoom_state$endpoints, list(other_point, initial_point))

  old_width <- abs(
    zoom_state$endpoints[[2L]]$alpha - zoom_state$endpoints[[1L]]$alpha
  )

  zoom_state <- update_schmidt_zoom(
    zoom_state,
    trial_point,
    initial_point,
    conditions
  )
  new_width <- abs(
    zoom_state$endpoints[[2L]]$alpha - zoom_state$endpoints[[1L]]$alpha
  )
  endpoint_values <- vapply(zoom_state$endpoints, `[[`, numeric(1L), "value")

  expect_lt(new_width, old_width)
  expect_equal(min(endpoint_values), trial_point$value)
})

test_that("Schmidt zoom preserves endpoint identity across value ties", {
  conditions <- make_line_condition_policy(0.1, 0.5)
  initial_point <- list(alpha = 0, value = 1, gradient = -1, slope = -1)
  best_point <- list(alpha = 2, value = 0, gradient = -1, slope = -1)
  tied_trial <- list(alpha = 1, value = 0, gradient = -1, slope = -1)
  zoom_state <- initialize_schmidt_zoom_state(list(initial_point, best_point))

  zoom_state <- update_schmidt_zoom(
    zoom_state,
    tied_trial,
    initial_point,
    conditions
  )

  expect_equal(zoom_state$endpoints, list(tied_trial, best_point))
  expect_identical(
    which.min(vapply(zoom_state$endpoints, `[[`, numeric(1L), "value")),
    1L
  )
})

test_that("bracket-and-zoom reports method-specific tolerance reasons", {
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
      value = 2,
      gradient = 1,
      slope = 1,
      parameters = alpha
    )
  }
  searches <- list(
    rasmussen = list(
      search = make_rasmussen_wolfe_search(
        armijo_constant = 0.05,
        curvature_constant = 0.1,
        relative_interval_tolerance = 10
      ),
      reason = "relative_interval_tolerance"
    ),
    schmidt = list(
      search = make_wolfe_line_search(
        core = run_bracket_zoom,
        armijo_constant = 0.05,
        curvature_constant = 0.1,
        method_policy = make_schmidt_wolfe_policy(parameter_tolerance = 10)
      ),
      reason = "parameter_tolerance"
    )
  )

  for (name in names(searches)) {
    case <- searches[[name]]
    result <- case$search(
      evaluate_line = evaluate_line,
      initial_point = initial_point,
      initial_alpha = 1,
      search_direction = 1
    )

    expect_identical(result$termination_reason, case$reason, info = name)
    expect_equal(result$line_point, initial_point, info = name)
  }
})
