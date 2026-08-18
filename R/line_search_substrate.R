# Common private substrate for Wolfe line searches.

validate_line_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop(name, " must be a numeric scalar")
  }
}

validate_line_evaluation_limit <- function(value, name) {
  validate_line_scalar(value, name)
  if (
    is.na(value) ||
      value < 0 ||
      (is.finite(value) && value != floor(value))
  ) {
    stop(name, " must be a nonnegative whole number or Inf")
  }
}

validate_line_point_structure <- function(line_point) {
  if (!is.list(line_point)) {
    stop("line point must be a list")
  }
  for (name in c("alpha", "value", "slope")) {
    if (
      !is.numeric(line_point[[name]]) ||
        length(line_point[[name]]) != 1L
    ) {
      stop("line point ", name, " must be a numeric scalar")
    }
  }
  if (!is.null(line_point$gradient) && !is.numeric(line_point$gradient)) {
    stop("line point gradient must be numeric or NULL")
  }
  if (!is.null(line_point$parameters) && !is.numeric(line_point$parameters)) {
    stop("line point parameters must be numeric or NULL")
  }
  line_point
}

validate_initial_line_point <- function(line_point) {
  line_point <- validate_line_point_structure(line_point)
  if (is.null(line_point$gradient)) {
    stop("initial line point gradient must be numeric")
  }
  line_point
}

validate_wolfe_trial_point <- function(line_point) {
  line_point <- validate_line_point_structure(line_point)
  if (is.null(line_point$gradient)) {
    stop("Wolfe trial point gradient must be numeric")
  }
  if (is.null(line_point$parameters)) {
    stop("Wolfe trial point parameters must be numeric")
  }
  line_point
}

initial_line_point_is_usable <- function(line_point) {
  isTRUE(is.finite(line_point$alpha)) &&
    isTRUE(is.finite(line_point$value)) &&
    isTRUE(is.finite(line_point$slope)) &&
    !is.null(line_point$gradient) &&
    all(is.finite(line_point$gradient)) &&
    (is.null(line_point$parameters) || all(is.finite(line_point$parameters)))
}

wolfe_trial_point_is_usable <- function(line_point) {
  initial_line_point_is_usable(line_point) &&
    !is.null(line_point$parameters) &&
    all(is.finite(line_point$parameters))
}

make_line_evaluator <- function(
  evaluate_line,
  initial_point = NULL,
  max_evaluations = Inf
) {
  if (!is.function(evaluate_line)) {
    stop("evaluate_line must be a function")
  }
  validate_line_evaluation_limit(max_evaluations, "max_evaluations")
  if (!is.null(initial_point)) {
    initial_point <- validate_initial_line_point(initial_point)
  }

  evaluation_count <- 0L
  best_decreasing_point <- NULL
  evaluator <- function(alpha, calc_gradient = TRUE) {
    if (evaluation_count >= max_evaluations) {
      stop("line evaluator budget exhausted")
    }
    line_point <- evaluate_line(alpha, calc_gradient = calc_gradient)
    evaluation_count <<- evaluation_count + 1L
    line_point <- validate_wolfe_trial_point(line_point)

    if (
      !is.null(initial_point) &&
        wolfe_trial_point_is_usable(line_point) &&
        line_point$value < initial_point$value &&
        (is.null(best_decreasing_point) ||
          line_point$value < best_decreasing_point$value)
    ) {
      best_decreasing_point <<- line_point
    }
    line_point
  }
  attr(evaluator, "line_evaluator") <- evaluator
  evaluator
}

make_line_condition_policy <- function(
  armijo_constant,
  curvature_constant,
  approximate_armijo = FALSE,
  strong_curvature = TRUE,
  approximation_tolerance = 1e-6
) {
  validate_line_scalar(armijo_constant, "armijo_constant")
  validate_line_scalar(curvature_constant, "curvature_constant")
  validate_line_scalar(approximation_tolerance, "approximation_tolerance")
  if (
    is.na(armijo_constant) ||
      armijo_constant < 0 ||
      armijo_constant > 1
  ) {
    stop("armijo_constant must be between zero and one")
  }
  if (
    is.na(curvature_constant) ||
      curvature_constant < armijo_constant ||
      curvature_constant > 1
  ) {
    stop(
      "curvature_constant must be between armijo_constant and one"
    )
  }
  if (is.na(approximation_tolerance) || approximation_tolerance < 0) {
    stop("approximation_tolerance must be nonnegative")
  }
  if (
    !is.logical(approximate_armijo) ||
      length(approximate_armijo) != 1L ||
      is.na(approximate_armijo)
  ) {
    stop("approximate_armijo must be TRUE or FALSE")
  }
  if (
    !is.logical(strong_curvature) ||
      length(strong_curvature) != 1L ||
      is.na(strong_curvature)
  ) {
    stop("strong_curvature must be TRUE or FALSE")
  }
  exact_armijo <- function(initial_point, trial_point, ...) {
    trial_point$value <=
      initial_point$value +
        armijo_constant * trial_point$alpha * initial_point$slope
  }
  selected_armijo <- if (approximate_armijo) {
    function(initial_point, trial_point, ...) {
      exact_armijo(initial_point, trial_point) ||
        (trial_point$value <=
          initial_point$value +
            approximation_tolerance * abs(initial_point$value) &&
          (2 * armijo_constant - 1) * initial_point$slope >= trial_point$slope)
    }
  } else {
    exact_armijo
  }

  selected_curvature <- if (strong_curvature) {
    function(initial_point, trial_point, ...) {
      abs(trial_point$slope) <= -curvature_constant * initial_point$slope
    }
  } else {
    function(initial_point, trial_point, ...) {
      trial_point$slope >= curvature_constant * initial_point$slope
    }
  }

  list(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    armijo = selected_armijo,
    curvature = selected_curvature,
    wolfe = function(initial_point, trial_point, ...) {
      selected_armijo(initial_point, trial_point) &&
        selected_curvature(initial_point, trial_point)
    }
  )
}

finalize_wolfe_line_search_result <- function(
  accepted_point,
  evaluator,
  termination_reason
) {
  evaluator_state <- environment(evaluator)
  accepted_point_is_usable <- !is.null(accepted_point) &&
    wolfe_trial_point_is_usable(accepted_point)

  if (accepted_point_is_usable && identical(termination_reason, "wolfe")) {
    selected_point <- accepted_point
  } else if (!is.null(evaluator_state$best_decreasing_point)) {
    selected_point <- evaluator_state$best_decreasing_point
  } else {
    selected_point <- evaluator_state$initial_point
  }

  evaluations <- evaluator_state$evaluation_count
  list(
    line_point = selected_point,
    function_evaluations = evaluations,
    gradient_evaluations = evaluations
  )
}

make_line_search_core_result <- function(
  termination_reason,
  accepted_point = NULL
) {
  list(
    accepted_point = accepted_point,
    termination_reason = termination_reason
  )
}

make_wolfe_line_search <- function(
  core,
  armijo_constant,
  curvature_constant,
  max_evaluations = Inf,
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE,
  method_policy = list()
) {
  validate_line_evaluation_limit(max_evaluations, "max_evaluations")
  evaluator <- make_line_evaluator(
    function(alpha, calc_gradient = TRUE) NULL,
    max_evaluations = 0
  )
  evaluator_state <- environment(evaluator)
  condition_policy <- make_line_condition_policy(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    approximation_tolerance = approximation_tolerance
  )

  function(
    evaluate_line,
    initial_point,
    initial_alpha,
    remaining_function_evaluations = Inf,
    remaining_gradient_evaluations = Inf,
    remaining_combined_evaluations = Inf,
    search_direction = NULL
  ) {
    initial_point <- validate_initial_line_point(initial_point)
    validate_line_scalar(initial_alpha, "initial_alpha")
    if (
      is.na(initial_alpha) ||
        !is.finite(initial_alpha) ||
        initial_alpha <= 0
    ) {
      stop("initial_alpha must be positive and finite")
    }
    evaluation_limit <- min(
      max_evaluations,
      remaining_function_evaluations,
      remaining_gradient_evaluations,
      floor(remaining_combined_evaluations / 2)
    )
    evaluator_state$evaluate_line <- evaluate_line
    evaluator_state$max_evaluations <- evaluation_limit
    evaluator_state$initial_point <- initial_point
    evaluator_state$evaluation_count <- 0L
    evaluator_state$best_decreasing_point <- NULL

    if (evaluation_limit <= 0) {
      core_result <- make_line_search_core_result("budget_exhausted")
    } else if (!initial_line_point_is_usable(initial_point)) {
      core_result <- make_line_search_core_result("nonfinite_initial_point")
    } else {
      core_result <- core(
        evaluator = evaluator,
        initial_point = initial_point,
        initial_alpha = initial_alpha,
        condition_policy = condition_policy,
        search_direction = search_direction,
        method_policy = method_policy
      )
    }

    if (
      identical(core_result$termination_reason, "wolfe") &&
        !wolfe_trial_point_is_usable(core_result$accepted_point)
    ) {
      core_result <- make_line_search_core_result("invalid_wolfe_point")
    }

    finalized_result <- finalize_wolfe_line_search_result(
      accepted_point = core_result$accepted_point,
      evaluator = evaluator,
      termination_reason = core_result$termination_reason
    )
    finalized_result$termination_reason <- core_result$termination_reason
    finalized_result
  }
}
