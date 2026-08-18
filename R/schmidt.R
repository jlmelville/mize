# Schmidt line searches implementing methods used by Mark Schmidt's minFunc.
# https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.

# Wolfe search ------------------------------------------------------------

make_schmidt_wolfe_search <- function(
  armijo_constant,
  curvature_constant,
  max_evaluations = Inf,
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE
) {
  method_policy <- make_schmidt_wolfe_policy()
  make_wolfe_line_search(
    core = run_bracket_zoom,
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    max_evaluations = max_evaluations,
    approximation_tolerance = approximation_tolerance,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    method_policy = method_policy
  )
}

make_schmidt_wolfe_policy <- function(
  expansion_factor = 10,
  minimum_expansion_fraction = 0.01,
  interior_fraction = 0.1,
  parameter_tolerance = 1e-9
) {
  validate_bracket_zoom_control(
    expansion_factor,
    "expansion_factor",
    1,
    minimum_open = TRUE
  )
  validate_bracket_zoom_control(
    minimum_expansion_fraction,
    "minimum_expansion_fraction",
    0
  )
  validate_bracket_zoom_control(
    interior_fraction,
    "interior_fraction",
    0,
    0.5,
    maximum_open = TRUE
  )
  validate_bracket_zoom_control(
    parameter_tolerance,
    "parameter_tolerance",
    0
  )

  list(
    expansion_recovery_lower_bound = function(expansion_state) {
      if (expansion_state$iteration == 0L) {
        0
      } else {
        expansion_state$previous_point$alpha
      }
    },
    classify_expansion = classify_schmidt_expansion,
    propose_expansion = function(expansion_state, trial_point) {
      minimum_alpha <- trial_point$alpha +
        minimum_expansion_fraction *
          (trial_point$alpha - expansion_state$previous_point$alpha)
      maximum_alpha <- trial_point$alpha * expansion_factor
      propose_schmidt_cubic_alpha(
        expansion_state$previous_point,
        trial_point,
        lower_alpha = minimum_alpha,
        upper_alpha = maximum_alpha
      )
    },
    initialize_zoom = initialize_schmidt_zoom_state,
    propose_zoom = function(zoom_state, initial_point) {
      propose_schmidt_zoom_alpha(zoom_state, interior_fraction)
    },
    zoom_recovery_lower_bound = function(zoom_state) {
      min(
        zoom_state$endpoints[[1L]]$alpha,
        zoom_state$endpoints[[2L]]$alpha
      )
    },
    zoom_bounds = function(zoom_state) {
      sort(c(
        zoom_state$endpoints[[1L]]$alpha,
        zoom_state$endpoints[[2L]]$alpha
      ))
    },
    process_zoom_trial = function(
      zoom_state,
      trial_point,
      initial_point,
      condition_policy,
      direction_scale
    ) {
      process_schmidt_zoom_trial(
        zoom_state,
        trial_point,
        initial_point,
        condition_policy,
        direction_scale,
        parameter_tolerance
      )
    }
  )
}

classify_schmidt_expansion <- function(
  expansion_state,
  trial_point,
  condition_policy
) {
  initial_point <- expansion_state$initial_point
  armijo_failed <- !condition_policy$armijo(initial_point, trial_point)
  objective_stopped_decreasing <- expansion_state$iteration >= 1L &&
    trial_point$value >= expansion_state$previous_point$value
  trial_satisfies_wolfe <- !armijo_failed &&
    condition_policy$curvature(initial_point, trial_point)
  trial_accepted <- trial_satisfies_wolfe &&
    !objective_stopped_decreasing
  bracket_found <- armijo_failed ||
    objective_stopped_decreasing ||
    (!trial_accepted && trial_point$slope >= 0)

  list(
    accepted = trial_accepted,
    bracket = if (bracket_found) {
      list(expansion_state$previous_point, trial_point)
    } else {
      NULL
    }
  )
}

initialize_schmidt_zoom_state <- function(bracket) {
  # Endpoint position breaks objective-value ties during later updates.
  list(
    endpoints = bracket,
    previous_proposal_not_safely_interior = FALSE
  )
}

propose_schmidt_zoom_alpha <- function(zoom_state, interior_fraction) {
  bracket_alphas <- c(
    zoom_state$endpoints[[1L]]$alpha,
    zoom_state$endpoints[[2L]]$alpha
  )
  lower_alpha <- min(bracket_alphas)
  upper_alpha <- max(bracket_alphas)
  interval_width <- upper_alpha - lower_alpha
  proposed_alpha <- propose_schmidt_cubic_alpha(
    zoom_state$endpoints[[1L]],
    zoom_state$endpoints[[2L]],
    lower_alpha = lower_alpha,
    upper_alpha = upper_alpha
  )
  if (!is.finite(proposed_alpha)) {
    proposed_alpha <- mean(bracket_alphas)
  }

  proposal_not_safely_interior <-
    interval_width > 0 &&
    min(
      upper_alpha - proposed_alpha,
      proposed_alpha - lower_alpha
    ) /
      interval_width <
      interior_fraction
  if (proposal_not_safely_interior) {
    proposal_not_strictly_interior <- proposed_alpha >= upper_alpha ||
      proposed_alpha <= lower_alpha
    if (
      zoom_state$previous_proposal_not_safely_interior ||
        proposal_not_strictly_interior
    ) {
      if (
        abs(proposed_alpha - upper_alpha) < abs(proposed_alpha - lower_alpha)
      ) {
        proposed_alpha <- upper_alpha - interior_fraction * interval_width
      } else {
        proposed_alpha <- lower_alpha + interior_fraction * interval_width
      }
      zoom_state$previous_proposal_not_safely_interior <- FALSE
    } else {
      zoom_state$previous_proposal_not_safely_interior <- TRUE
    }
  } else {
    zoom_state$previous_proposal_not_safely_interior <- FALSE
  }

  list(alpha = proposed_alpha, state = zoom_state)
}

update_schmidt_zoom <- function(
  zoom_state,
  trial_point,
  initial_point,
  condition_policy
) {
  best_index <- which.min(c(
    zoom_state$endpoints[[1L]]$value,
    zoom_state$endpoints[[2L]]$value
  ))
  other_index <- 3L - best_index
  best_point <- zoom_state$endpoints[[best_index]]
  other_point <- zoom_state$endpoints[[other_index]]

  if (
    !condition_policy$armijo(initial_point, trial_point) ||
      trial_point$value >= best_point$value
  ) {
    zoom_state$endpoints[[other_index]] <- trial_point
  } else {
    if (
      trial_point$slope *
        (other_point$alpha - best_point$alpha) >=
        0
    ) {
      zoom_state$endpoints[[other_index]] <- best_point
    }
    zoom_state$endpoints[[best_index]] <- trial_point
  }
  zoom_state
}

process_schmidt_zoom_trial <- function(
  zoom_state,
  trial_point,
  initial_point,
  condition_policy,
  direction_scale,
  parameter_tolerance
) {
  zoom_state <- update_schmidt_zoom(
    zoom_state,
    trial_point,
    initial_point,
    condition_policy
  )
  alpha_interval_width <- abs(
    zoom_state$endpoints[[2L]]$alpha -
      zoom_state$endpoints[[1L]]$alpha
  )
  list(
    state = zoom_state,
    termination_reason = if (
      alpha_interval_width * direction_scale < parameter_tolerance
    ) {
      "parameter_tolerance"
    } else {
      NULL
    }
  )
}

# Armijo search ------------------------------------------------------------

make_schmidt_armijo_search <- function(
  armijo_constant = 0.05,
  step_down = NULL,
  max_evaluations = Inf,
  parameter_tolerance = 1e-9
) {
  validate_line_scalar(armijo_constant, "armijo_constant")
  validate_line_evaluation_limit(max_evaluations, "max_evaluations")
  validate_line_scalar(parameter_tolerance, "parameter_tolerance")
  if (
    is.na(armijo_constant) ||
      !is.finite(armijo_constant) ||
      armijo_constant < 0 ||
      armijo_constant > 1
  ) {
    stop("armijo_constant must be between zero and one")
  }
  if (
    is.na(parameter_tolerance) ||
      !is.finite(parameter_tolerance) ||
      parameter_tolerance < 0
  ) {
    stop("parameter_tolerance must be nonnegative and finite")
  }
  if (!is.null(step_down)) {
    validate_line_scalar(step_down, "step_down")
    if (
      is.na(step_down) ||
        !is.finite(step_down) ||
        step_down < 0 ||
        step_down > 1
    ) {
      stop("step_down must be between zero and one")
    }
  }

  fixed_reduction_factor <- step_down
  evaluates_gradient <- is.null(fixed_reduction_factor)

  function(
    evaluate_line,
    initial_point,
    initial_alpha,
    remaining_function_evaluations = Inf,
    remaining_gradient_evaluations = Inf,
    remaining_combined_evaluations = Inf,
    search_direction
  ) {
    evaluation_limit <- if (evaluates_gradient) {
      min(
        max_evaluations,
        remaining_function_evaluations,
        remaining_gradient_evaluations,
        floor(remaining_combined_evaluations / 2)
      )
    } else {
      min(
        max_evaluations,
        remaining_function_evaluations,
        remaining_combined_evaluations
      )
    }
    if (evaluation_limit <= 0) {
      return(list(
        line_point = initial_point,
        function_evaluations = 0L,
        gradient_evaluations = 0L,
        gradient_is_current = !is.null(initial_point$gradient),
        termination_reason = "budget_exhausted"
      ))
    }

    run_schmidt_armijo_search(
      evaluate_line = evaluate_line,
      initial_point = initial_point,
      initial_alpha = initial_alpha,
      armijo_constant = armijo_constant,
      fixed_reduction_factor = fixed_reduction_factor,
      max_evaluations = evaluation_limit,
      direction_scale = norm_inf(search_direction),
      parameter_tolerance = parameter_tolerance
    )
  }
}

run_schmidt_armijo_search <- function(
  evaluate_line,
  initial_point,
  initial_alpha,
  armijo_constant,
  fixed_reduction_factor,
  max_evaluations,
  direction_scale,
  parameter_tolerance
) {
  evaluates_gradient <- is.null(fixed_reduction_factor)
  function_evaluations <- 0L
  gradient_evaluations <- 0L
  best_decreasing_point <- NULL
  trial_alpha <- initial_alpha
  is_initial_trial <- TRUE

  repeat {
    trial_point <- evaluate_line(
      trial_alpha,
      calc_gradient = evaluates_gradient
    )
    function_evaluations <- function_evaluations + 1L
    if (evaluates_gradient) {
      gradient_evaluations <- gradient_evaluations + 1L
    }

    trial_point_is_usable <- schmidt_armijo_point_is_usable(trial_point)
    has_strict_decrease <- trial_point_is_usable &&
      isTRUE(trial_point$value < initial_point$value)
    if (
      has_strict_decrease &&
        (is.null(best_decreasing_point) ||
          trial_point$value < best_decreasing_point$value)
    ) {
      best_decreasing_point <- trial_point
    }

    if (
      has_strict_decrease &&
        isTRUE(line_point_satisfies_armijo(
          initial_point,
          trial_point,
          armijo_constant
        ))
    ) {
      return(finalize_schmidt_armijo_result(
        line_point = trial_point,
        function_evaluations = function_evaluations,
        gradient_evaluations = gradient_evaluations,
        termination_reason = "armijo"
      ))
    }

    if (
      !is_initial_trial &&
        direction_scale * trial_alpha <= parameter_tolerance
    ) {
      termination_reason <- "parameter_tolerance"
      break
    }
    if (function_evaluations >= max_evaluations) {
      termination_reason <- "budget_exhausted"
      break
    }

    proposed_alpha <- propose_schmidt_armijo_alpha(
      initial_point = initial_point,
      trial_point = trial_point,
      fixed_reduction_factor = fixed_reduction_factor
    )
    trial_alpha <- safeguard_schmidt_armijo_alpha(
      proposed_alpha,
      previous_alpha = trial_alpha
    )
    is_initial_trial <- FALSE
  }

  selected_point <- if (is.null(best_decreasing_point)) {
    initial_point
  } else {
    best_decreasing_point
  }
  finalize_schmidt_armijo_result(
    line_point = selected_point,
    function_evaluations = function_evaluations,
    gradient_evaluations = gradient_evaluations,
    termination_reason = termination_reason
  )
}

finalize_schmidt_armijo_result <- function(
  line_point,
  function_evaluations,
  gradient_evaluations,
  termination_reason
) {
  list(
    line_point = line_point,
    function_evaluations = function_evaluations,
    gradient_evaluations = gradient_evaluations,
    gradient_is_current = !is.null(line_point$gradient),
    termination_reason = termination_reason
  )
}

schmidt_armijo_point_is_usable <- function(line_point) {
  isTRUE(is.finite(line_point$alpha)) &&
    isTRUE(is.finite(line_point$value)) &&
    (is.null(line_point$gradient) || all(is.finite(line_point$gradient))) &&
    (is.null(line_point$slope) || isTRUE(is.finite(line_point$slope))) &&
    (is.null(line_point$parameters) || all(is.finite(line_point$parameters)))
}

propose_schmidt_armijo_alpha <- function(
  initial_point,
  trial_point,
  fixed_reduction_factor
) {
  if (!isTRUE(is.finite(trial_point$value))) {
    reduction_factor <- if (is.null(fixed_reduction_factor)) {
      0.5
    } else {
      fixed_reduction_factor
    }
    return(reduction_factor * trial_point$alpha)
  }
  if (!is.null(fixed_reduction_factor)) {
    return(fixed_reduction_factor * trial_point$alpha)
  }
  if (
    is.null(trial_point$gradient) ||
      any(!is.finite(trial_point$gradient)) ||
      !isTRUE(is.finite(trial_point$slope))
  ) {
    proposed_alpha <- propose_quadratic_alpha(initial_point, trial_point)
    return(min(max(proposed_alpha, 0), trial_point$alpha))
  }

  propose_schmidt_cubic_alpha(
    initial_point,
    trial_point,
    lower_alpha = 0,
    upper_alpha = trial_point$alpha
  )
}

safeguard_schmidt_armijo_alpha <- function(
  proposed_alpha,
  previous_alpha,
  minimum_fraction = 1e-3,
  maximum_fraction = 0.6
) {
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(previous_alpha * maximum_fraction)
  }
  min(
    max(proposed_alpha, previous_alpha * minimum_fraction),
    previous_alpha * maximum_fraction
  )
}

# Schmidt cubic proposal ---------------------------------------------------

propose_schmidt_cubic_alpha <- function(
  first_point,
  second_point,
  lower_alpha = min(first_point$alpha, second_point$alpha),
  upper_alpha = max(first_point$alpha, second_point$alpha)
) {
  if (first_point$alpha <= second_point$alpha) {
    lower_point <- first_point
    upper_point <- second_point
  } else {
    lower_point <- second_point
    upper_point <- first_point
  }
  if (lower_point$alpha == upper_point$alpha) {
    return(lower_point$alpha)
  }
  midpoint_alpha <- lower_alpha + (upper_alpha - lower_alpha) / 2

  cubic_shape <- lower_point$slope +
    upper_point$slope -
    3 *
      (lower_point$value - upper_point$value) /
      (lower_point$alpha - upper_point$alpha)
  discriminant <- cubic_shape^2 - lower_point$slope * upper_point$slope
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return(midpoint_alpha)
  }

  discriminant_root <- sqrt(discriminant)
  denominator <- upper_point$slope - lower_point$slope + 2 * discriminant_root
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(midpoint_alpha)
  }
  proposed_alpha <- upper_point$alpha -
    (upper_point$alpha - lower_point$alpha) *
      ((upper_point$slope + discriminant_root - cubic_shape) / denominator)
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(midpoint_alpha)
  }

  min(max(proposed_alpha, lower_alpha), upper_alpha)
}
