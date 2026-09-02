# Rasmussen Wolfe line search.
#
# @seealso This implementation was informed by Carl Edward Rasmussen's
#  `minimize.m` routine in the Matlab
#  [GPML](https://www.gaussianprocess.org/gpml/code/matlab/doc/) package.

make_rasmussen_wolfe_search <- function(
  armijo_constant,
  curvature_constant,
  max_evaluations = Inf,
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE,
  interior_fraction = 0.1,
  expansion_factor = 3,
  relative_interval_tolerance = 1e-6
) {
  method_policy <- make_rasmussen_wolfe_policy(
    interior_fraction = interior_fraction,
    expansion_factor = expansion_factor,
    relative_interval_tolerance = relative_interval_tolerance
  )
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

make_rasmussen_wolfe_policy <- function(
  interior_fraction = 0.1,
  expansion_factor = 3,
  relative_interval_tolerance = 1e-6
) {
  list(
    expansion_recovery_lower_bound = function(expansion_state) 0,
    classify_expansion = classify_rasmussen_expansion,
    propose_expansion = function(expansion_state, trial_point) {
      propose_rasmussen_expansion_alpha(
        expansion_state$initial_point,
        trial_point,
        expansion_factor,
        interior_fraction
      )
    },
    initialize_zoom = initialize_rasmussen_zoom_state,
    propose_zoom = function(zoom_state, initial_point) {
      list(
        alpha = propose_rasmussen_zoom_alpha(
          zoom_state,
          initial_point,
          interior_fraction
        ),
        state = zoom_state
      )
    },
    zoom_recovery_lower_bound = function(zoom_state) {
      zoom_state$lower_point$alpha
    },
    zoom_bounds = function(zoom_state) {
      c(zoom_state$lower_point$alpha, zoom_state$upper_point$alpha)
    },
    zoom_endpoints = function(zoom_state) {
      list(zoom_state$lower_point, zoom_state$upper_point)
    },
    process_zoom_trial = function(
      zoom_state,
      trial_point,
      initial_point,
      condition_policy,
      direction_scale
    ) {
      process_rasmussen_zoom_trial(
        zoom_state,
        trial_point,
        initial_point,
        condition_policy,
        relative_interval_tolerance
      )
    }
  )
}

classify_rasmussen_expansion <- function(
  expansion_state,
  trial_point,
  condition_policy
) {
  initial_point <- expansion_state$initial_point
  armijo_failed <- !condition_policy$selected_armijo(initial_point, trial_point)
  trial_accepted <- !armijo_failed &&
    condition_policy$curvature(initial_point, trial_point)
  bracket_found <- trial_point$slope >
    condition_policy$curvature_constant * initial_point$slope ||
    armijo_failed

  list(
    accepted = trial_accepted,
    bracket = if (bracket_found) {
      list(initial_point, trial_point)
    } else {
      NULL
    }
  )
}

propose_rasmussen_expansion_alpha <- function(
  initial_point,
  trial_point,
  expansion_factor,
  interior_fraction
) {
  proposed_alpha <- propose_rasmussen_expansion_cubic_alpha(
    initial_point,
    trial_point
  )
  minimum_alpha <- trial_point$alpha +
    interior_fraction * (trial_point$alpha - initial_point$alpha)
  maximum_alpha <- trial_point$alpha * expansion_factor

  if (!isTRUE(is.finite(proposed_alpha)) || proposed_alpha < 0) {
    return(maximum_alpha)
  }
  min(max(proposed_alpha, minimum_alpha), maximum_alpha)
}

propose_rasmussen_expansion_cubic_alpha <- function(
  first_point,
  second_point
) {
  alpha_difference <- second_point$alpha - first_point$alpha
  if (!isTRUE(is.finite(alpha_difference)) || alpha_difference == 0) {
    return(NA_real_)
  }
  scaled_cubic_term <- 6 *
    (first_point$value - second_point$value) +
    3 * (second_point$slope + first_point$slope) * alpha_difference
  scaled_quadratic_term <- 3 *
    (second_point$value - first_point$value) -
    (2 * first_point$slope + second_point$slope) * alpha_difference
  discriminant <- scaled_quadratic_term^2 -
    scaled_cubic_term * first_point$slope * alpha_difference
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return(NA_real_)
  }

  denominator <- scaled_quadratic_term + sqrt(discriminant)
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(NA_real_)
  }
  proposed_alpha <- first_point$alpha -
    first_point$slope * alpha_difference^2 / denominator
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(NA_real_)
  }
  proposed_alpha
}

initialize_rasmussen_zoom_state <- function(bracket) {
  if (bracket[[1L]]$alpha <= bracket[[2L]]$alpha) {
    lower_point <- bracket[[1L]]
    upper_point <- bracket[[2L]]
  } else {
    lower_point <- bracket[[2L]]
    upper_point <- bracket[[1L]]
  }
  list(lower_point = lower_point, upper_point = upper_point)
}

propose_rasmussen_zoom_alpha <- function(
  zoom_state,
  initial_point,
  interior_fraction
) {
  proposed_alpha <- if (zoom_state$upper_point$value > initial_point$value) {
    propose_quadratic_alpha(zoom_state$lower_point, zoom_state$upper_point)
  } else {
    cubic_interpolate(
      zoom_state$lower_point$alpha,
      zoom_state$lower_point$value,
      zoom_state$lower_point$slope,
      zoom_state$upper_point$alpha,
      zoom_state$upper_point$value,
      zoom_state$upper_point$slope
    )
  }
  safeguard_rasmussen_zoom_alpha(
    proposed_alpha,
    zoom_state$lower_point$alpha,
    zoom_state$upper_point$alpha,
    interior_fraction
  )
}

safeguard_rasmussen_zoom_alpha <- function(
  proposed_alpha,
  lower_alpha,
  upper_alpha,
  interior_fraction
) {
  if (!isTRUE(is.finite(proposed_alpha))) {
    proposed_alpha <- lower_alpha + (upper_alpha - lower_alpha) / 2
  }
  interval_width <- upper_alpha - lower_alpha
  min(
    max(proposed_alpha, lower_alpha + interior_fraction * interval_width),
    upper_alpha - interior_fraction * interval_width
  )
}

update_rasmussen_zoom <- function(
  zoom_state,
  trial_point,
  initial_point,
  condition_policy
) {
  if (
    trial_point$slope > 0 ||
      !condition_policy$selected_armijo(initial_point, trial_point)
  ) {
    zoom_state$upper_point <- trial_point
  } else {
    zoom_state$lower_point <- trial_point
  }
  zoom_state
}

process_rasmussen_zoom_trial <- function(
  zoom_state,
  trial_point,
  initial_point,
  condition_policy,
  relative_interval_tolerance
) {
  relative_interval_tolerance_reached <- rasmussen_zoom_interval_is_small(
    zoom_state,
    trial_point,
    relative_interval_tolerance
  )
  list(
    state = update_rasmussen_zoom(
      zoom_state,
      trial_point,
      initial_point,
      condition_policy
    ),
    termination_reason = if (relative_interval_tolerance_reached) {
      "relative_interval_tolerance"
    } else {
      NULL
    }
  )
}

rasmussen_zoom_interval_is_small <- function(
  zoom_state,
  trial_point,
  relative_interval_tolerance
) {
  abs(zoom_state$upper_point$alpha - zoom_state$lower_point$alpha) <
    relative_interval_tolerance * trial_point$alpha
}
