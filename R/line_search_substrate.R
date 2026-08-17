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

validate_line_step <- function(step) {
  if (!is.list(step)) {
    stop("line step must be a list")
  }
  for (name in c("alpha", "f", "d")) {
    if (!is.numeric(step[[name]]) || length(step[[name]]) != 1L) {
      stop("line step ", name, " must be a numeric scalar")
    }
  }
  if (!is.null(step$df) && !is.numeric(step$df)) {
    stop("line step df must be numeric or NULL")
  }
  if (!is.null(step$par) && !is.numeric(step$par)) {
    stop("line step par must be numeric or NULL")
  }
  step
}

line_step_is_usable <- function(step) {
  isTRUE(is.finite(step$alpha)) &&
    isTRUE(is.finite(step$f)) &&
    isTRUE(is.finite(step$d)) &&
    (is.null(step$df) || all(is.finite(step$df))) &&
    (is.null(step$par) || all(is.finite(step$par)))
}

new_line_evaluator <- function(
  phi,
  initial_step = NULL,
  max_evaluations = Inf
) {
  if (!is.function(phi)) {
    stop("phi must be a function")
  }
  validate_line_evaluation_limit(max_evaluations, "max_evaluations")
  if (!is.null(initial_step)) {
    initial_step <- validate_line_step(initial_step)
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
    step <- validate_line_step(step)

    has_finite_values <- isTRUE(is.finite(step$f)) &&
      isTRUE(is.finite(step$d))
    if (!has_finite_values) {
      recovered_nonfinite <<- TRUE
    }
    if (
      !is.null(initial_step) &&
        isTRUE(is.finite(step$alpha)) &&
        has_finite_values &&
        (is.null(step$df) || all(is.finite(step$df))) &&
        (is.null(step$par) || all(is.finite(step$par))) &&
        step$f < initial_step$f &&
        (is.null(best_decrease) || step$f < best_decrease$f)
    ) {
      best_decrease <<- step
    }
    step
  }
  attr(evaluator, "line_evaluator") <- evaluator
  evaluator
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
  exact_armijo <- function(initial_step, trial_step, ...) {
    trial_step$f <=
      initial_step$f +
        armijo_constant * trial_step$alpha * initial_step$d
  }
  selected_armijo <- if (approximate_armijo) {
    function(initial_step, trial_step, ...) {
      exact_armijo(initial_step, trial_step) ||
        (trial_step$f <=
          initial_step$f +
            approximation_tolerance * abs(initial_step$f) &&
          (2 * armijo_constant - 1) * initial_step$d >= trial_step$d)
    }
  } else {
    exact_armijo
  }

  selected_curvature <- if (strong_curvature) {
    function(initial_step, trial_step, ...) {
      abs(trial_step$d) <= -curvature_constant * initial_step$d
    }
  } else {
    function(initial_step, trial_step, ...) {
      trial_step$d > curvature_constant * initial_step$d
    }
  }

  list(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    armijo = selected_armijo,
    curvature = selected_curvature,
    wolfe = function(initial_step, trial_step, ...) {
      selected_armijo(initial_step, trial_step) &&
        selected_curvature(initial_step, trial_step)
    }
  )
}

finalize_wolfe_line_search_result <- function(
  candidate,
  evaluator,
  termination_reason,
  gradient_is_current = NULL
) {
  evaluator_state <- environment(evaluator)
  candidate_is_usable <- !is.null(candidate) &&
    line_step_is_usable(candidate)

  if (candidate_is_usable && identical(termination_reason, "wolfe")) {
    selected_step <- candidate
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
  method_policy = list()
) {
  validate_line_evaluation_limit(max_evaluations, "max_evaluations")
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
    step0 <- validate_line_step(step0)
    evaluation_limit <- min(
      max_evaluations,
      total_max_fn,
      total_max_gr,
      floor(total_max_fg / 2)
    )
    evaluator_state$phi <- phi
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
        initial_step = step0,
        initial_alpha = alpha,
        condition_policy = condition_policy,
        search_direction = pm,
        method_policy = method_policy
      )
    }

    finalize_wolfe_line_search_result(
      candidate = core_result$candidate,
      evaluator = evaluator,
      termination_reason = core_result$termination_reason,
      gradient_is_current = core_result$gradient_is_current
    )
  }
}
