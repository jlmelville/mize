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

make_parameter_resolution_searches <- function(
  armijo_constant,
  curvature_constant,
  max_evaluations = Inf
) {
  list(
    more_thuente = make_wolfe_line_search(
      more_thuente_core,
      armijo_constant = armijo_constant,
      curvature_constant = curvature_constant,
      max_evaluations = max_evaluations,
      method_policy = make_more_thuente_policy()
    ),
    rasmussen = make_rasmussen_wolfe_search(
      armijo_constant,
      curvature_constant,
      max_evaluations = max_evaluations
    ),
    schmidt = make_schmidt_wolfe_search(
      armijo_constant,
      curvature_constant,
      max_evaluations = max_evaluations
    ),
    hager_zhang = make_hager_zhang_search(
      armijo_constant,
      curvature_constant,
      max_evaluations = max_evaluations,
      strong_curvature = TRUE,
      approximate_armijo = FALSE
    )
  )
}

test_that("bracket-and-zoom accepts a colliding alpha without a callback", {
  initial_parameters <- 1
  search_direction <- -2^-53
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -0.5 / search_direction,
    slope = -0.5,
    parameters = initial_parameters
  )
  endpoint_parameters <- project_line_parameters(
    initial_parameters,
    2,
    search_direction
  )
  method_policies <- list(
    rasmussen = make_rasmussen_wolfe_policy(),
    schmidt = make_schmidt_wolfe_policy()
  )

  for (method_name in names(method_policies)) {
    method_policy <- method_policies[[method_name]]
    method_policy$classify_expansion <- function(
      expansion_state,
      trial_point,
      condition_policy
    ) {
      bracket <- if (expansion_state$iteration == 1L) {
        list(expansion_state$previous_point, trial_point)
      } else {
        NULL
      }
      list(accepted = FALSE, bracket = bracket)
    }
    method_policy$propose_expansion <- function(
      expansion_state,
      trial_point
    ) {
      2
    }
    method_policy$propose_zoom <- function(zoom_state, initial_point) {
      list(alpha = 1.5, state = zoom_state)
    }
    callback_alphas <- numeric()
    evaluate <- function(alpha, calc_gradient = TRUE) {
      callback_alphas <<- c(callback_alphas, alpha)
      parameters <- project_line_parameters(
        initial_parameters,
        alpha,
        search_direction
      )
      if (alpha == 1) {
        return(list(
          alpha = alpha,
          value = 0.9,
          gradient = -0.4 / search_direction,
          slope = -0.4,
          parameters = parameters
        ))
      }
      list(
        alpha = alpha,
        value = 0.8,
        gradient = 0,
        slope = 0,
        parameters = parameters
      )
    }
    search <- make_wolfe_line_search(
      run_bracket_zoom,
      armijo_constant = 0.25,
      curvature_constant = 0.5,
      max_evaluations = 4,
      method_policy = method_policy
    )

    result <- search(
      evaluate_line = evaluate,
      initial_point = initial_point,
      initial_alpha = 1,
      search_direction = search_direction
    )

    expect_identical(
      result$termination_reason,
      "wolfe",
      info = method_name
    )
    expect_identical(result$outcome, "wolfe", info = method_name)
    expect_identical(result$function_evaluations, 2L, info = method_name)
    expect_identical(result$gradient_evaluations, 2L, info = method_name)
    expect_equal(callback_alphas, c(1, 2), info = method_name)
    expect_equal(result$line_point$alpha, 1.5, info = method_name)
    expect_equal(
      result$line_point$parameters,
      endpoint_parameters,
      info = method_name
    )
  }
})

test_that("parameter-exhausted Wolfe brackets stop without repeated callbacks", {
  initial_parameters <- 1
  search_direction <- -2^-54
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -0.5 / search_direction,
    slope = -0.5,
    parameters = initial_parameters
  )
  endpoint_parameters <- initial_parameters + 2 * search_direction
  searches <- make_parameter_resolution_searches(0.25, 0.5)

  for (search_name in names(searches)) {
    callback_parameters <- numeric()
    evaluate <- function(alpha, calc_gradient = TRUE) {
      parameters <- initial_parameters + alpha * search_direction
      callback_parameters <<- c(callback_parameters, parameters)
      list(
        alpha = alpha,
        value = if (parameters == initial_parameters) 1 else 1.1,
        gradient = if (parameters == initial_parameters) {
          initial_point$gradient
        } else {
          0
        },
        slope = if (parameters == initial_parameters) -0.5 else 0,
        parameters = parameters
      )
    }

    result <- searches[[search_name]](
      evaluate_line = evaluate,
      initial_point = initial_point,
      initial_alpha = 2,
      search_direction = search_direction
    )

    expect_identical(
      result$termination_reason,
      "rounding_stagnation",
      info = search_name
    )
    expect_identical(result$outcome, "no_step", info = search_name)
    expect_identical(result$line_point, initial_point, info = search_name)
    expect_identical(result$function_evaluations, 1L, info = search_name)
    expect_identical(result$gradient_evaluations, 1L, info = search_name)
    expect_equal(callback_parameters, endpoint_parameters, info = search_name)
  }
})

test_that("bracket-and-zoom recovers a novel vector after a collision", {
  initial_parameters <- 1
  search_direction <- -2^-52
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1 / search_direction,
    slope = -1,
    parameters = initial_parameters
  )
  method_policy <- make_rasmussen_wolfe_policy()
  method_policy$classify_expansion <- function(
    expansion_state,
    trial_point,
    condition_policy
  ) {
    list(accepted = FALSE, bracket = list(initial_point, trial_point))
  }
  method_policy$propose_zoom <- function(zoom_state, initial_point) {
    list(alpha = 0.25, state = zoom_state)
  }
  evaluated_parameters <- numeric()
  evaluate <- function(alpha, calc_gradient = TRUE) {
    parameters <- initial_parameters + alpha * search_direction
    evaluated_parameters <<- c(evaluated_parameters, parameters)
    if (alpha == 2) {
      return(list(
        alpha = alpha,
        value = 2,
        gradient = 1 / search_direction,
        slope = 1,
        parameters = parameters
      ))
    }
    list(
      alpha = alpha,
      value = 0.5,
      gradient = 0,
      slope = 0,
      parameters = parameters
    )
  }
  search <- make_wolfe_line_search(
    run_bracket_zoom,
    armijo_constant = 0.1,
    curvature_constant = 0.5,
    max_evaluations = 3,
    method_policy = method_policy
  )

  result <- search(
    evaluate,
    initial_point,
    initial_alpha = 2,
    search_direction = search_direction
  )

  expect_identical(result$termination_reason, "wolfe")
  expect_equal(result$line_point$alpha, 1)
  expect_identical(result$function_evaluations, 2L)
  expect_equal(
    evaluated_parameters,
    c(
      initial_parameters + 2 * search_direction,
      initial_parameters + search_direction
    )
  )
})
