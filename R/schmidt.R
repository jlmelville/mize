# Schmidt line searches implementing methods used by Mark Schmidt's minFunc.
# https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.

# Wolfe search ------------------------------------------------------------

new_schmidt_wolfe_search <- function(
  armijo_constant,
  curvature_constant,
  max_evaluations = Inf,
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE
) {
  method_policy <- new_schmidt_wolfe_policy()
  new_wolfe_line_search(
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

new_schmidt_wolfe_policy <- function(
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
        expansion_state$previous_step$alpha
      }
    },
    classify_expansion = classify_schmidt_expansion,
    propose_expansion = function(expansion_state, trial_step) {
      minimum_alpha <- trial_step$alpha +
        minimum_expansion_fraction *
          (trial_step$alpha - expansion_state$previous_step$alpha)
      maximum_alpha <- trial_step$alpha * expansion_factor
      propose_schmidt_cubic_alpha(
        expansion_state$previous_step,
        trial_step,
        lower_alpha = minimum_alpha,
        upper_alpha = maximum_alpha
      )
    },
    initialize_zoom = initialize_schmidt_zoom_state,
    propose_zoom = function(zoom_state, initial_step) {
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
      trial_step,
      initial_step,
      condition_policy,
      direction_scale
    ) {
      process_schmidt_zoom_trial(
        zoom_state,
        trial_step,
        initial_step,
        condition_policy,
        direction_scale,
        parameter_tolerance
      )
    }
  )
}

classify_schmidt_expansion <- function(
  expansion_state,
  trial_step,
  condition_policy
) {
  initial_step <- expansion_state$initial_step
  armijo_failed <- !condition_policy$armijo(initial_step, trial_step)
  objective_stopped_decreasing <- expansion_state$iteration >= 1L &&
    trial_step$f >= expansion_state$previous_step$f
  trial_satisfies_wolfe <- !armijo_failed &&
    condition_policy$curvature(initial_step, trial_step)
  trial_accepted <- trial_satisfies_wolfe &&
    !objective_stopped_decreasing
  bracket_found <- armijo_failed ||
    objective_stopped_decreasing ||
    (!trial_accepted && trial_step$d >= 0)

  list(
    accepted = trial_accepted,
    bracket = if (bracket_found) {
      list(expansion_state$previous_step, trial_step)
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
  trial_step,
  initial_step,
  condition_policy
) {
  best_index <- which.min(c(
    zoom_state$endpoints[[1L]]$f,
    zoom_state$endpoints[[2L]]$f
  ))
  other_index <- 3L - best_index
  best_step <- zoom_state$endpoints[[best_index]]
  other_step <- zoom_state$endpoints[[other_index]]

  if (
    !condition_policy$armijo(initial_step, trial_step) ||
      trial_step$f >= best_step$f
  ) {
    zoom_state$endpoints[[other_index]] <- trial_step
  } else {
    if (
      trial_step$d *
        (other_step$alpha - best_step$alpha) >=
        0
    ) {
      zoom_state$endpoints[[other_index]] <- best_step
    }
    zoom_state$endpoints[[best_index]] <- trial_step
  }
  zoom_state
}

process_schmidt_zoom_trial <- function(
  zoom_state,
  trial_step,
  initial_step,
  condition_policy,
  direction_scale,
  parameter_tolerance
) {
  zoom_state <- update_schmidt_zoom(
    zoom_state,
    trial_step,
    initial_step,
    condition_policy
  )
  alpha_interval_width <- abs(
    zoom_state$endpoints[[2L]]$alpha -
      zoom_state$endpoints[[1L]]$alpha
  )
  list(
    state = zoom_state,
    progress_stalled = alpha_interval_width * direction_scale <
      parameter_tolerance
  )
}

# Armijo search ------------------------------------------------------------

new_schmidt_armijo_search <- function(
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
    phi,
    step0,
    alpha,
    total_max_fn = Inf,
    total_max_gr = Inf,
    total_max_fg = Inf,
    pm
  ) {
    evaluation_limit <- if (evaluates_gradient) {
      min(
        max_evaluations,
        total_max_fn,
        total_max_gr,
        floor(total_max_fg / 2)
      )
    } else {
      min(max_evaluations, total_max_fn, total_max_fg)
    }
    if (evaluation_limit <= 0) {
      return(list(
        step = step0,
        nfn = 0L,
        ngr = 0L,
        is_gr_curr = !is.null(step0$df)
      ))
    }

    run_schmidt_armijo_search(
      phi = phi,
      initial_step = step0,
      initial_alpha = alpha,
      armijo_constant = armijo_constant,
      fixed_reduction_factor = fixed_reduction_factor,
      max_evaluations = evaluation_limit,
      direction_scale = norm_inf(pm),
      parameter_tolerance = parameter_tolerance
    )
  }
}

run_schmidt_armijo_search <- function(
  phi,
  initial_step,
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
  best_decreasing_step <- NULL
  trial_alpha <- initial_alpha
  is_initial_trial <- TRUE

  repeat {
    trial_step <- phi(trial_alpha, calc_gradient = evaluates_gradient)
    function_evaluations <- function_evaluations + 1L
    if (evaluates_gradient) {
      gradient_evaluations <- gradient_evaluations + 1L
    }

    trial_step_is_usable <- schmidt_armijo_step_is_usable(trial_step)
    has_strict_decrease <- trial_step_is_usable &&
      isTRUE(trial_step$f < initial_step$f)
    if (
      has_strict_decrease &&
        (is.null(best_decreasing_step) ||
          trial_step$f < best_decreasing_step$f)
    ) {
      best_decreasing_step <- trial_step
    }

    if (
      has_strict_decrease &&
        isTRUE(armijo_ok_step(initial_step, trial_step, armijo_constant))
    ) {
      return(finalize_schmidt_armijo_result(
        step = trial_step,
        function_evaluations = function_evaluations,
        gradient_evaluations = gradient_evaluations
      ))
    }

    if (
      !is_initial_trial &&
        direction_scale * trial_alpha <= parameter_tolerance
    ) {
      break
    }
    if (function_evaluations >= max_evaluations) {
      break
    }

    proposed_alpha <- propose_schmidt_armijo_alpha(
      initial_step = initial_step,
      trial_step = trial_step,
      fixed_reduction_factor = fixed_reduction_factor
    )
    trial_alpha <- safeguard_schmidt_armijo_alpha(
      proposed_alpha,
      previous_alpha = trial_alpha
    )
    is_initial_trial <- FALSE
  }

  selected_step <- if (is.null(best_decreasing_step)) {
    initial_step
  } else {
    best_decreasing_step
  }
  finalize_schmidt_armijo_result(
    step = selected_step,
    function_evaluations = function_evaluations,
    gradient_evaluations = gradient_evaluations
  )
}

finalize_schmidt_armijo_result <- function(
  step,
  function_evaluations,
  gradient_evaluations
) {
  list(
    step = step,
    nfn = function_evaluations,
    ngr = gradient_evaluations,
    is_gr_curr = !is.null(step$df)
  )
}

schmidt_armijo_step_is_usable <- function(step) {
  isTRUE(is.finite(step$alpha)) &&
    isTRUE(is.finite(step$f)) &&
    (is.null(step$df) || all(is.finite(step$df))) &&
    (is.null(step$d) || isTRUE(is.finite(step$d))) &&
    (is.null(step$par) || all(is.finite(step$par)))
}

propose_schmidt_armijo_alpha <- function(
  initial_step,
  trial_step,
  fixed_reduction_factor
) {
  if (!isTRUE(is.finite(trial_step$f))) {
    reduction_factor <- if (is.null(fixed_reduction_factor)) {
      0.5
    } else {
      fixed_reduction_factor
    }
    return(reduction_factor * trial_step$alpha)
  }
  if (!is.null(fixed_reduction_factor)) {
    return(fixed_reduction_factor * trial_step$alpha)
  }
  if (
    is.null(trial_step$df) ||
      any(!is.finite(trial_step$df)) ||
      !isTRUE(is.finite(trial_step$d))
  ) {
    proposed_alpha <- propose_quadratic_alpha(initial_step, trial_step)
    return(min(max(proposed_alpha, 0), trial_step$alpha))
  }

  propose_schmidt_cubic_alpha(
    initial_step,
    trial_step,
    lower_alpha = 0,
    upper_alpha = trial_step$alpha
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
  first_step,
  second_step,
  lower_alpha = min(first_step$alpha, second_step$alpha),
  upper_alpha = max(first_step$alpha, second_step$alpha)
) {
  if (first_step$alpha <= second_step$alpha) {
    lower_step <- first_step
    upper_step <- second_step
  } else {
    lower_step <- second_step
    upper_step <- first_step
  }
  if (lower_step$alpha == upper_step$alpha) {
    return(lower_step$alpha)
  }
  midpoint_alpha <- lower_alpha + (upper_alpha - lower_alpha) / 2

  cubic_shape <- lower_step$d +
    upper_step$d -
    3 *
      (lower_step$f - upper_step$f) /
      (lower_step$alpha - upper_step$alpha)
  discriminant <- cubic_shape^2 - lower_step$d * upper_step$d
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return(midpoint_alpha)
  }

  discriminant_root <- sqrt(discriminant)
  denominator <- upper_step$d - lower_step$d + 2 * discriminant_root
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(midpoint_alpha)
  }
  proposed_alpha <- upper_step$alpha -
    (upper_step$alpha - lower_step$alpha) *
      ((upper_step$d + discriminant_root - cubic_shape) / denominator)
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(midpoint_alpha)
  }

  min(max(proposed_alpha, lower_alpha), upper_alpha)
}
