# Common private substrate for Wolfe line searches.

validate_line_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop(name, " must be a numeric scalar")
  }
}

new_line_point <- function(
  alpha,
  value,
  gradient = NULL,
  slope,
  parameters = NULL
) {
  validate_line_scalar(alpha, "alpha")
  validate_line_scalar(value, "value")
  validate_line_scalar(slope, "slope")
  if (!is.null(gradient) && !is.numeric(gradient)) {
    stop("gradient must be numeric or NULL")
  }
  if (!is.null(parameters) && !is.numeric(parameters)) {
    stop("parameters must be numeric or NULL")
  }

  list(
    alpha = alpha,
    value = value,
    gradient = gradient,
    slope = slope,
    parameters = parameters
  )
}

as_line_point <- function(point) {
  if (!is.list(point)) {
    stop("line point must be a list")
  }
  if (
    identical(
      names(point),
      c("alpha", "value", "gradient", "slope", "parameters")
    )
  ) {
    return(point)
  }
  alpha <- point$alpha
  value <- point$f
  gradient <- point$df
  slope <- point$d
  parameters <- point$par
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      !is.numeric(value) ||
      length(value) != 1L ||
      !is.numeric(slope) ||
      length(slope) != 1L ||
      (!is.null(gradient) && !is.numeric(gradient)) ||
      (!is.null(parameters) && !is.numeric(parameters))
  ) {
    stop("invalid line point")
  }
  list(
    alpha = alpha,
    value = value,
    gradient = gradient,
    slope = slope,
    parameters = parameters
  )
}

as_legacy_line_step <- function(point) {
  step <- list(alpha = point$alpha, f = point$value)
  if (!is.null(point$gradient)) {
    step$df <- point$gradient
  }
  step$d <- point$slope
  if (!is.null(point$parameters)) {
    step$par <- point$parameters
  }
  step
}

line_point_has_finite_values <- function(point) {
  isTRUE(is.finite(point$value)) && isTRUE(is.finite(point$slope))
}

line_point_is_usable <- function(point) {
  isTRUE(is.finite(point$alpha)) &&
    line_point_has_finite_values(point) &&
    (is.null(point$gradient) || all(is.finite(point$gradient))) &&
    (is.null(point$parameters) || all(is.finite(point$parameters)))
}

new_line_evaluator <- function(
  phi,
  initial_point = NULL,
  max_evaluations = Inf,
  initial_step = NULL
) {
  if (!is.function(phi)) {
    stop("phi must be a function")
  }
  if (
    !is.numeric(max_evaluations) ||
      length(max_evaluations) != 1L ||
      is.na(max_evaluations) ||
      max_evaluations < 0
  ) {
    stop("max_evaluations must be a nonnegative numeric scalar")
  }
  if (!is.null(initial_point) && is.null(initial_point$value)) {
    initial_point <- as_line_point(initial_point)
  }

  evaluation_count <- 0L
  best_decrease <- NULL
  recovered_nonfinite <- FALSE
  evaluator <- function(alpha, calc_gradient = TRUE) {
    if (evaluation_count >= max_evaluations) {
      stop("line evaluator budget exhausted")
    }
    step <- phi(alpha, calc_gradient = calc_gradient)
    evaluation_count <<- evaluation_count + 1L

    has_finite_values <- is.finite(step$f) && is.finite(step$d)
    if (!has_finite_values) {
      recovered_nonfinite <<- TRUE
    }
    if (
      !is.null(initial_point) &&
        is.finite(step$alpha) &&
        has_finite_values &&
        (is.null(step$df) || all(is.finite(step$df))) &&
        (is.null(step$par) || all(is.finite(step$par))) &&
        step$f < initial_point$value &&
        (is.null(best_decrease) || step$f < best_decrease$f)
    ) {
      best_decrease <<- step
    }
    step
  }
  attr(evaluator, "line_evaluator") <- evaluator
  evaluator
}

remaining_line_evaluations <- function(evaluator) {
  state <- environment(evaluator)
  max(0, state$max_evaluations - state$evaluation_count)
}

line_evaluator_evaluations <- function(evaluator) {
  environment(evaluator)$evaluation_count
}

line_evaluator_best_decrease <- function(evaluator) {
  environment(evaluator)$best_decrease
}

evaluate_line_point <- function(evaluator, alpha, calc_gradient = TRUE) {
  as_line_point(evaluate_legacy_line_step(evaluator, alpha, calc_gradient))
}

evaluate_legacy_line_step <- function(evaluator, alpha, calc_gradient = TRUE) {
  step <- evaluator(alpha, calc_gradient)
  if (
    !is.list(step) ||
      !is.numeric(step$alpha) ||
      length(step$alpha) != 1L ||
      !is.numeric(step$f) ||
      length(step$f) != 1L ||
      !is.numeric(step$d) ||
      length(step$d) != 1L ||
      (!is.null(step$df) && !is.numeric(step$df)) ||
      (!is.null(step$par) && !is.numeric(step$par))
  ) {
    stop("invalid line point")
  }
  step
}

recover_finite_line_point <- function(
  evaluator,
  alpha,
  minimum_alpha = 0,
  maximum_evaluations = Inf
) {
  recovered <- recover_finite_legacy_step(
    evaluator,
    alpha,
    minimum_alpha,
    maximum_evaluations
  )
  recovered$point <- if (is.null(recovered$step)) {
    NULL
  } else {
    as_line_point(recovered$step)
  }
  recovered$step <- NULL
  recovered
}

recover_finite_legacy_step <- function(
  evaluator,
  alpha,
  minimum_alpha = 0,
  maximum_evaluations = Inf
) {
  evaluator_state <- environment(evaluator)
  evaluations_before <- evaluator_state$evaluation_count
  step <- NULL
  ok <- FALSE
  while (
    evaluator_state$evaluation_count - evaluations_before <
      maximum_evaluations &&
      evaluator_state$evaluation_count < evaluator_state$max_evaluations &&
      alpha >= minimum_alpha
  ) {
    step <- evaluate_legacy_line_step(evaluator, alpha)
    if (isTRUE(is.finite(step$f)) && isTRUE(is.finite(step$d))) {
      ok <- TRUE
      break
    }
    alpha <- (minimum_alpha + alpha) / 2
  }
  list(
    step = step,
    evaluations = evaluator_state$evaluation_count - evaluations_before,
    ok = ok
  )
}

new_line_condition_policy <- function(
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
    trial_point$f <=
      initial_point$f +
        armijo_constant * trial_point$alpha * initial_point$d
  }
  selected_armijo <- if (approximate_armijo) {
    function(initial_point, trial_point, ...) {
      exact_armijo(initial_point, trial_point) ||
        (trial_point$f <=
          initial_point$f +
            approximation_tolerance * abs(initial_point$f) &&
          (2 * armijo_constant - 1) * initial_point$d >= trial_point$d)
    }
  } else {
    exact_armijo
  }

  selected_curvature <- if (strong_curvature) {
    function(initial_point, trial_point, ...) {
      abs(trial_point$d) <= -curvature_constant * initial_point$d
    }
  } else {
    function(initial_point, trial_point, ...) {
      trial_point$d > curvature_constant * initial_point$d
    }
  }

  classify <- function(initial_point, trial_point) {
    if (is.null(initial_point$value)) {
      initial_value <- initial_point$f
      initial_slope <- initial_point$d
    } else {
      initial_value <- initial_point$value
      initial_slope <- initial_point$slope
    }
    if (is.null(trial_point$value)) {
      trial_alpha <- trial_point$alpha
      trial_value <- trial_point$f
      trial_slope <- trial_point$d
    } else {
      trial_alpha <- trial_point$alpha
      trial_value <- trial_point$value
      trial_slope <- trial_point$slope
    }

    exact_armijo <- isTRUE(
      trial_value <=
        initial_value + armijo_constant * trial_alpha * initial_slope
    )
    approximate_armijo_result <- isTRUE(
      trial_value <=
        initial_value + approximation_tolerance * abs(initial_value)
    ) &&
      isTRUE(
        (2 * armijo_constant - 1) * initial_slope >= trial_slope
      )
    armijo <- exact_armijo ||
      (approximate_armijo && approximate_armijo_result)
    weak_curvature <- isTRUE(
      trial_slope > curvature_constant * initial_slope
    )
    strong_curvature_result <- isTRUE(
      abs(trial_slope) <= -curvature_constant * initial_slope
    )
    curvature <- if (strong_curvature) {
      strong_curvature_result
    } else {
      weak_curvature
    }

    list(
      exact_armijo = exact_armijo,
      approximate_armijo = approximate_armijo_result,
      armijo = armijo,
      weak_curvature = weak_curvature,
      strong_curvature = strong_curvature_result,
      curvature = curvature,
      wolfe = armijo && curvature
    )
  }

  list(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    approximation_tolerance = approximation_tolerance,
    classify = classify,
    armijo = selected_armijo,
    curvature = selected_curvature,
    wolfe = function(initial_point, trial_point, ...) {
      selected_armijo(initial_point, trial_point) &&
        selected_curvature(initial_point, trial_point)
    }
  )
}

finalize_line_search_result <- function(
  candidate,
  initial_point,
  evaluator,
  condition_policy,
  termination_reason = NULL,
  gradient_is_current = NULL
) {
  evaluator_state <- environment(evaluator)
  candidate_step <- candidate
  if (!is.null(candidate) && is.null(candidate$value)) {
    candidate <- list(
      alpha = candidate$alpha,
      value = candidate$f,
      gradient = candidate$df,
      slope = candidate$d,
      parameters = candidate$par
    )
  }

  accepted <- !is.null(candidate) &&
    line_point_is_usable(candidate) &&
    identical(termination_reason, "wolfe")

  if (accepted) {
    point <- candidate
    selected_step <- if (
      identical(
        names(candidate_step),
        c("alpha", "value", "gradient", "slope", "parameters")
      )
    ) {
      as_legacy_line_step(candidate_step)
    } else {
      candidate_step
    }
    selection <- "wolfe"
    termination_reason <- "wolfe"
  } else if (!is.null(best_decrease <- evaluator_state$best_decrease)) {
    point <- list(
      alpha = best_decrease$alpha,
      value = best_decrease$f,
      gradient = best_decrease$df,
      slope = best_decrease$d,
      parameters = best_decrease$par
    )
    selected_step <- best_decrease
    selection <- "strict_decrease"
  } else {
    point <- initial_point
    selected_step <- evaluator_state$initial_step
    if (is.null(selected_step)) {
      selected_step <- as_legacy_line_step(initial_point)
    }
    selection <- "unchanged_start"
  }

  if (is.null(termination_reason)) {
    if (evaluator_state$evaluation_count >= evaluator_state$max_evaluations) {
      termination_reason <- "budget_exhausted"
    } else if (evaluator_state$recovered_nonfinite) {
      termination_reason <- "nonfinite_recovery"
    } else if (selection == "unchanged_start") {
      termination_reason <- "unchanged_start"
    } else {
      termination_reason <- "progress_failure"
    }
  }

  if (accepted) {
    selected_armijo <- TRUE
    selected_curvature <- TRUE
  } else {
    classification <- condition_policy$classify(initial_point, point)
    selected_armijo <- classification$armijo
    selected_curvature <- classification$curvature
  }

  evaluations <- evaluator_state$evaluation_count
  legacy_result <- list(
    step = selected_step,
    nfn = evaluations,
    ngr = evaluations
  )
  if (!is.null(gradient_is_current)) {
    legacy_result$is_gr_curr <- gradient_is_current
  }

  list(
    point = point,
    accepted = accepted,
    selection = selection,
    termination_reason = termination_reason,
    conditions = list(
      armijo = selected_armijo,
      curvature = selected_curvature,
      wolfe = selected_armijo && selected_curvature
    ),
    evaluations = evaluations,
    recovered_nonfinite = evaluator_state$recovered_nonfinite,
    gradient_is_current = gradient_is_current,
    legacy_result = legacy_result
  )
}

finalize_legacy_line_search_result <- function(
  candidate,
  evaluator,
  termination_reason,
  gradient_is_current = NULL
) {
  evaluator_state <- environment(evaluator)
  candidate_step <- candidate
  if (!is.null(candidate) && !is.null(candidate$value)) {
    candidate_step <- as_legacy_line_step(candidate)
  }

  candidate_is_usable <- !is.null(candidate_step) &&
    is.finite(candidate_step$alpha) &&
    is.finite(candidate_step$f) &&
    is.finite(candidate_step$d) &&
    (is.null(candidate_step$df) || all(is.finite(candidate_step$df))) &&
    (is.null(candidate_step$par) || all(is.finite(candidate_step$par)))

  if (candidate_is_usable && identical(termination_reason, "wolfe")) {
    selected_step <- candidate_step
  } else if (!is.null(evaluator_state$best_decrease)) {
    selected_step <- evaluator_state$best_decrease
  } else {
    selected_step <- evaluator_state$initial_step
  }

  evaluations <- evaluator_state$evaluation_count
  result <- list(step = selected_step, nfn = evaluations, ngr = evaluations)
  if (!is.null(gradient_is_current)) {
    result$is_gr_curr <- gradient_is_current
  }
  result
}

new_wolfe_line_search <- function(
  core,
  armijo_constant,
  curvature_constant,
  max_evaluations = Inf,
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE,
  options = list()
) {
  evaluator <- new_line_evaluator(
    function(alpha, calc_gradient = TRUE) NULL,
    max_evaluations = 0
  )
  evaluator_state <- environment(evaluator)
  condition_policy <- new_line_condition_policy(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    approximation_tolerance = approximation_tolerance
  )

  function(
    phi,
    step0,
    alpha,
    total_max_fn = Inf,
    total_max_gr = Inf,
    total_max_fg = Inf,
    pm = NULL
  ) {
    evaluation_limit <- min(
      max_evaluations,
      total_max_fn,
      total_max_gr,
      floor(total_max_fg / 2)
    )
    initial_point <- list(
      alpha = step0$alpha,
      value = step0$f,
      gradient = step0$df,
      slope = step0$d,
      parameters = step0$par
    )
    evaluator_state$phi <- phi
    evaluator_state$initial_point <- initial_point
    evaluator_state$max_evaluations <- evaluation_limit
    evaluator_state$initial_step <- step0
    evaluator_state$evaluation_count <- 0L
    evaluator_state$best_decrease <- NULL
    evaluator_state$recovered_nonfinite <- FALSE

    if (evaluation_limit <= 0) {
      core_result <- list(
        candidate = NULL,
        termination_reason = "budget_exhausted"
      )
    } else {
      core_result <- core(
        evaluator = evaluator,
        initial_point = initial_point,
        initial_alpha = alpha,
        condition_policy = condition_policy,
        direction = pm,
        options = options
      )
    }

    finalize_legacy_line_search_result(
      candidate = core_result$candidate,
      evaluator = evaluator,
      termination_reason = core_result$termination_reason,
      gradient_is_current = core_result$gradient_is_current
    )
  }
}
