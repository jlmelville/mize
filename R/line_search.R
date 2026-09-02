# Shared private helpers for Wolfe line searches.

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
    !is.null(line_point$parameters)
}

project_line_parameters <- function(parameters, alpha, search_direction) {
  parameters + (alpha * search_direction)
}

line_parameters_are_identical <- function(first, second) {
  line_numeric_vector_is_finite(first) &&
    line_numeric_vector_is_finite(second) &&
    line_parameters_have_same_values(first, second)
}

line_numeric_vector_is_finite <- function(parameters) {
  is.numeric(parameters) &&
    is.null(dim(parameters)) &&
    all(is.finite(parameters))
}

line_parameters_have_same_values <- function(first, second) {
  if (
    !is.numeric(first) ||
      !is.numeric(second) ||
      length(first) != length(second)
  ) {
    return(FALSE)
  }
  if (
    identical(typeof(first), typeof(second)) &&
      is.null(attributes(first)) &&
      is.null(attributes(second))
  ) {
    return(identical(first, second))
  }
  isTRUE(all(first == second))
}

line_parameter_mapping_is_available <- function(
  initial_parameters,
  search_direction,
  first_endpoint,
  second_endpoint
) {
  vectors <- list(
    initial_parameters,
    search_direction,
    first_endpoint$parameters,
    second_endpoint$parameters
  )
  vector_lengths <- vapply(vectors, length, integer(1L))
  vectors_are_conformable <- all(vapply(
    vectors,
    line_numeric_vector_is_finite,
    logical(1L)
  )) &&
    length(unique(vector_lengths)) == 1L
  endpoint_alphas_are_finite <- is_finite_numeric(first_endpoint$alpha) &&
    is_finite_numeric(second_endpoint$alpha)
  if (!vectors_are_conformable || !endpoint_alphas_are_finite) {
    return(FALSE)
  }

  first_projection <- project_line_parameters(
    initial_parameters,
    first_endpoint$alpha,
    search_direction
  )
  second_projection <- project_line_parameters(
    initial_parameters,
    second_endpoint$alpha,
    search_direction
  )
  line_parameters_have_same_values(
    first_projection,
    first_endpoint$parameters
  ) &&
    line_parameters_have_same_values(
      second_projection,
      second_endpoint$parameters
    )
}

line_search_midpoint <- function(lower_alpha, upper_alpha) {
  midpoint <- lower_alpha + (upper_alpha - lower_alpha) / 2
  if (!isTRUE(is.finite(midpoint))) {
    midpoint <- lower_alpha / 2 + upper_alpha / 2
  }
  if (
    !isTRUE(is.finite(midpoint)) ||
      !isTRUE(midpoint > lower_alpha) ||
      !isTRUE(midpoint < upper_alpha)
  ) {
    return(NULL)
  }
  midpoint
}

find_novel_bracketed_line_parameters <- function(
  first_endpoint,
  second_endpoint,
  initial_parameters,
  search_direction
) {
  if (first_endpoint$alpha <= second_endpoint$alpha) {
    lower_endpoint <- first_endpoint
    upper_endpoint <- second_endpoint
  } else {
    lower_endpoint <- second_endpoint
    upper_endpoint <- first_endpoint
  }
  lower_alpha <- lower_endpoint$alpha
  upper_alpha <- upper_endpoint$alpha

  repeat {
    midpoint <- line_search_midpoint(lower_alpha, upper_alpha)
    if (is.null(midpoint)) {
      return(NULL)
    }
    midpoint_parameters <- project_line_parameters(
      initial_parameters,
      midpoint,
      search_direction
    )
    if (
      line_parameters_have_same_values(
        midpoint_parameters,
        lower_endpoint$parameters
      )
    ) {
      lower_alpha <- midpoint
    } else if (
      line_parameters_have_same_values(
        midpoint_parameters,
        upper_endpoint$parameters
      )
    ) {
      upper_alpha <- midpoint
    } else {
      return(list(alpha = midpoint, parameters = midpoint_parameters))
    }
  }
}

bracketed_line_parameters_are_exhausted <- function(
  first_endpoint,
  second_endpoint,
  initial_parameters,
  search_direction
) {
  line_parameter_mapping_is_available(
    initial_parameters,
    search_direction,
    first_endpoint,
    second_endpoint
  ) &&
    is.null(find_novel_bracketed_line_parameters(
      first_endpoint,
      second_endpoint,
      initial_parameters,
      search_direction
    ))
}

make_condition_only_line_point <- function(alpha, parameters, endpoint) {
  list(
    alpha = alpha,
    value = endpoint$value,
    gradient = endpoint$gradient,
    slope = endpoint$slope,
    parameters = parameters
  )
}

admit_bracketed_line_alpha <- function(
  alpha,
  first_endpoint,
  second_endpoint,
  initial_parameters,
  search_direction
) {
  unchecked <- list(
    status = "unchecked",
    alpha = alpha,
    parameters = NULL,
    matching_endpoint = NULL,
    condition_point = NULL,
    replacement = NULL,
    exhausted = FALSE
  )
  if (
    !is_finite_numeric(alpha) ||
      !wolfe_trial_point_is_usable(first_endpoint) ||
      !wolfe_trial_point_is_usable(second_endpoint) ||
      !line_parameter_mapping_is_available(
        initial_parameters,
        search_direction,
        first_endpoint,
        second_endpoint
      )
  ) {
    return(unchecked)
  }

  parameters <- project_line_parameters(
    initial_parameters,
    alpha,
    search_direction
  )
  matches <- c(
    line_parameters_have_same_values(parameters, first_endpoint$parameters),
    line_parameters_have_same_values(parameters, second_endpoint$parameters)
  )
  if (!any(matches)) {
    return(list(
      status = "novel",
      alpha = alpha,
      parameters = parameters,
      matching_endpoint = NULL,
      condition_point = NULL,
      replacement = NULL,
      exhausted = FALSE
    ))
  }

  matching_index <- which(matches)
  if (length(matching_index) > 1L) {
    alpha_distances <- abs(c(
      first_endpoint$alpha - alpha,
      second_endpoint$alpha - alpha
    ))
    matching_index <- matching_index[[which.min(alpha_distances[matches])]]
  } else {
    matching_index <- matching_index[[1L]]
  }
  matching_endpoint <- if (matching_index == 1L) {
    first_endpoint
  } else {
    second_endpoint
  }
  list(
    status = "collision",
    alpha = alpha,
    parameters = parameters,
    matching_endpoint = if (matching_index == 1L) "first" else "second",
    condition_point = make_condition_only_line_point(
      alpha,
      parameters,
      matching_endpoint
    ),
    replacement = NULL,
    exhausted = NULL
  )
}

make_line_evaluator <- function(
  evaluate_line,
  initial_point = NULL,
  max_evaluations = Inf
) {
  if (!is.function(evaluate_line)) {
    stop("evaluate_line must be a function")
  }
  evaluation_count <- 0L
  best_decreasing_point <- NULL
  parameter_projection_witness_initialized <- FALSE
  parameter_projection_witness <- NULL
  evaluator <- function(alpha, calc_gradient = TRUE) {
    if (evaluation_count >= max_evaluations) {
      stop("line evaluator budget exhausted")
    }
    line_point <- evaluate_line(alpha, calc_gradient = calc_gradient)
    evaluation_count <<- evaluation_count + 1L
    validate_wolfe_trial_point(line_point)
  }
  evaluator
}

initialize_line_evaluator_parameter_witness <- function(
  evaluator,
  initial_point,
  search_direction
) {
  evaluator_state <- environment(evaluator)
  initial_parameters <- initial_point$parameters
  witness_is_available <-
    line_numeric_vector_is_finite(initial_parameters) &&
    line_numeric_vector_is_finite(search_direction) &&
    length(initial_parameters) == length(search_direction) &&
    length(search_direction) > 0L &&
    any(search_direction != 0)
  witness <- NULL
  if (witness_is_available) {
    witness_index <- which.max(abs(search_direction))
    witness <- list(
      index = witness_index,
      initial_parameter = initial_parameters[[witness_index]],
      direction = search_direction[[witness_index]]
    )
  }
  evaluator_state$parameter_projection_witness_initialized <- TRUE
  evaluator_state$parameter_projection_witness <- witness
  invisible(witness)
}

bracketed_line_proposal_is_definitely_novel <- function(
  evaluator,
  alpha,
  first_endpoint,
  second_endpoint
) {
  evaluator_state <- environment(evaluator)
  if (!isTRUE(evaluator_state$parameter_projection_witness_initialized)) {
    return(FALSE)
  }
  witness <- evaluator_state$parameter_projection_witness
  if (is.null(witness)) {
    return(FALSE)
  }
  witness_index <- witness$index
  if (
    witness_index > length(first_endpoint$parameters) ||
      witness_index > length(second_endpoint$parameters)
  ) {
    return(FALSE)
  }
  projected_witness <- witness$initial_parameter + alpha * witness$direction
  isTRUE(is.finite(projected_witness)) &&
    projected_witness != first_endpoint$parameters[[witness_index]] &&
    projected_witness != second_endpoint$parameters[[witness_index]]
}

continue_bracketed_line_recovery <- function(
  evaluator,
  recovery,
  failed_alpha,
  min_alpha,
  first_endpoint,
  second_endpoint,
  initial_point,
  condition_policy,
  search_direction,
  max_evaluations
) {
  current_parameters <- project_line_parameters(
    initial_point$parameters,
    failed_alpha,
    search_direction
  )
  callback_mapping_is_consistent <-
    !is.null(recovery$line_point) &&
    line_parameters_have_same_values(
      current_parameters,
      recovery$line_point$parameters
    )
  next_alpha <- line_search_midpoint(min_alpha, failed_alpha)
  if (!callback_mapping_is_consistent || is.null(next_alpha)) {
    recovery$accepted <- FALSE
    recovery$termination_reason <- NULL
    return(recovery)
  }
  remaining_evaluations <- if (is.finite(max_evaluations)) {
    max(0, max_evaluations - 1)
  } else {
    max_evaluations
  }
  recover_bracketed_line_point_slow(
    evaluator = evaluator,
    alpha = next_alpha,
    min_alpha = min_alpha,
    first_endpoint = first_endpoint,
    second_endpoint = second_endpoint,
    initial_point = initial_point,
    condition_policy = condition_policy,
    search_direction = search_direction,
    max_evaluations = remaining_evaluations,
    initial_nonfinite_recovery = recovery,
    initial_nonfinite_alpha = failed_alpha,
    initial_nonfinite_parameters = current_parameters
  )
}

recover_bracketed_line_point <- function(
  evaluator,
  alpha,
  min_alpha,
  first_endpoint,
  second_endpoint,
  initial_point,
  condition_policy,
  search_direction,
  max_evaluations = 20
) {
  evaluator_state <- environment(evaluator)
  final_evaluation_count <- min(
    evaluator_state$max_evaluations,
    evaluator_state$evaluation_count + max_evaluations
  )
  witness_initialized <-
    evaluator_state$parameter_projection_witness_initialized
  if (isTRUE(witness_initialized)) {
    witness <- evaluator_state$parameter_projection_witness
    if (is.null(witness)) {
      recovery <- recover_finite_line_point(
        evaluator,
        alpha,
        min_alpha = min_alpha,
        max_evaluations = max_evaluations
      )
      recovery$accepted <- FALSE
      recovery$termination_reason <- NULL
      return(recovery)
    }
    if (
      bracketed_line_proposal_is_definitely_novel(
        evaluator,
        alpha,
        first_endpoint,
        second_endpoint
      )
    ) {
      recovery <- recover_finite_line_point(
        evaluator,
        alpha,
        min_alpha = alpha,
        max_evaluations = 1
      )
      if (recovery$succeeded) {
        recovery$accepted <- FALSE
        recovery$termination_reason <- NULL
        return(recovery)
      }
      return(continue_bracketed_line_recovery(
        evaluator = evaluator,
        recovery = recovery,
        failed_alpha = alpha,
        min_alpha = min_alpha,
        first_endpoint = first_endpoint,
        second_endpoint = second_endpoint,
        initial_point = initial_point,
        condition_policy = condition_policy,
        search_direction = search_direction,
        max_evaluations = final_evaluation_count -
          evaluator_state$evaluation_count +
          1
      ))
    }
  }

  recover_bracketed_line_point_slow(
    evaluator = evaluator,
    alpha = alpha,
    min_alpha = min_alpha,
    first_endpoint = first_endpoint,
    second_endpoint = second_endpoint,
    initial_point = initial_point,
    condition_policy = condition_policy,
    search_direction = search_direction,
    max_evaluations = max_evaluations
  )
}

recover_bracketed_line_point_slow <- function(
  evaluator,
  alpha,
  min_alpha,
  first_endpoint,
  second_endpoint,
  initial_point,
  condition_policy,
  search_direction,
  max_evaluations,
  initial_nonfinite_recovery = NULL,
  initial_nonfinite_alpha = NULL,
  initial_nonfinite_parameters = NULL
) {
  evaluator_state <- environment(evaluator)
  final_evaluation_count <- min(
    evaluator_state$max_evaluations,
    evaluator_state$evaluation_count + max_evaluations
  )
  if (first_endpoint$alpha <= second_endpoint$alpha) {
    resolution_lower_endpoint <- first_endpoint
    resolution_upper_endpoint <- second_endpoint
  } else {
    resolution_lower_endpoint <- second_endpoint
    resolution_upper_endpoint <- first_endpoint
  }
  recovering_nonfinite_point <- !is.null(initial_nonfinite_recovery)
  last_recovery <- initial_nonfinite_recovery
  if (recovering_nonfinite_point) {
    resolution_upper_endpoint <- list(
      alpha = initial_nonfinite_alpha,
      parameters = initial_nonfinite_parameters
    )
  }
  current_alpha <- alpha

  repeat {
    admission <- admit_bracketed_line_alpha(
      alpha = current_alpha,
      first_endpoint = first_endpoint,
      second_endpoint = second_endpoint,
      initial_parameters = initial_point$parameters,
      search_direction = search_direction
    )
    if (identical(admission$status, "unchecked")) {
      recovery <- recover_finite_line_point(
        evaluator,
        current_alpha,
        min_alpha = min_alpha,
        max_evaluations = final_evaluation_count -
          evaluator_state$evaluation_count
      )
      recovery$accepted <- FALSE
      recovery$termination_reason <- NULL
      return(recovery)
    }
    if (identical(admission$status, "collision")) {
      if (condition_policy$wolfe(initial_point, admission$condition_point)) {
        return(list(
          line_point = admission$condition_point,
          succeeded = TRUE,
          accepted = TRUE,
          termination_reason = "wolfe"
        ))
      }
      replacement <- find_novel_bracketed_line_parameters(
        resolution_lower_endpoint,
        resolution_upper_endpoint,
        initial_point$parameters,
        search_direction
      )
      if (is.null(replacement)) {
        if (recovering_nonfinite_point) {
          last_recovery$accepted <- FALSE
          last_recovery$termination_reason <- NULL
          return(last_recovery)
        }
        return(list(
          line_point = NULL,
          succeeded = FALSE,
          accepted = FALSE,
          termination_reason = "rounding_stagnation"
        ))
      }
      current_alpha <- replacement$alpha
      current_parameters <- replacement$parameters
    } else {
      current_alpha <- admission$alpha
      current_parameters <- admission$parameters
    }

    if (evaluator_state$evaluation_count >= final_evaluation_count) {
      return(list(
        line_point = NULL,
        succeeded = FALSE,
        accepted = FALSE,
        termination_reason = NULL
      ))
    }
    recovery <- recover_finite_line_point(
      evaluator,
      current_alpha,
      min_alpha = current_alpha,
      max_evaluations = 1
    )
    if (recovery$succeeded) {
      recovery$accepted <- FALSE
      recovery$termination_reason <- NULL
      return(recovery)
    }
    recovering_nonfinite_point <- TRUE
    last_recovery <- recovery

    resolution_upper_endpoint <- list(
      alpha = current_alpha,
      parameters = current_parameters
    )
    next_alpha <- line_search_midpoint(
      resolution_lower_endpoint$alpha,
      resolution_upper_endpoint$alpha
    )
    if (is.null(next_alpha)) {
      recovery$accepted <- FALSE
      recovery$termination_reason <- NULL
      return(recovery)
    }
    current_alpha <- next_alpha
  }
}

make_line_condition_policy <- function(
  armijo_constant,
  curvature_constant,
  approximate_armijo = FALSE,
  strong_curvature = TRUE,
  approximation_tolerance = 1e-6
) {
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
    exact_armijo = exact_armijo,
    selected_armijo = selected_armijo,
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

  if (!is.null(accepted_point) && identical(termination_reason, "wolfe")) {
    selected_point <- accepted_point
    outcome <- "wolfe"
  } else if (!is.null(evaluator_state$best_decreasing_point)) {
    selected_point <- evaluator_state$best_decreasing_point
    outcome <- "improving_fallback"
  } else {
    selected_point <- evaluator_state$initial_point
    outcome <- "no_step"
  }

  evaluations <- evaluator_state$evaluation_count
  list(
    line_point = selected_point,
    function_evaluations = evaluations,
    gradient_evaluations = evaluations,
    outcome = outcome
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
    initialize_line_evaluator_parameter_witness(
      evaluator,
      initial_point,
      search_direction
    )

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

    accepted_point_is_usable <-
      !is.null(core_result$accepted_point) &&
      wolfe_trial_point_is_usable(core_result$accepted_point)
    if (
      identical(core_result$termination_reason, "wolfe") &&
        !accepted_point_is_usable
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
