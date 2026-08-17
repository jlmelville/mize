# Rasmussen Wolfe line search.
#
# @seealso This implementation was informed by Carl Edward Rasmussen's
#  `minimize.m` routine in the Matlab
#  [GPML](https://www.gaussianprocess.org/gpml/code/matlab/doc/) package.

new_rasmussen_wolfe_search <- function(
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
  method_policy <- new_rasmussen_wolfe_policy(
    interior_fraction = interior_fraction,
    expansion_factor = expansion_factor,
    relative_interval_tolerance = relative_interval_tolerance
  )
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

new_rasmussen_wolfe_policy <- function(
  interior_fraction = 0.1,
  expansion_factor = 3,
  relative_interval_tolerance = 1e-6
) {
  validate_bracket_zoom_control(
    interior_fraction,
    "interior_fraction",
    0,
    0.5,
    maximum_open = TRUE
  )
  validate_bracket_zoom_control(
    expansion_factor,
    "expansion_factor",
    1,
    minimum_open = TRUE
  )
  validate_bracket_zoom_control(
    relative_interval_tolerance,
    "relative_interval_tolerance",
    0
  )

  list(
    expansion_recovery_lower_bound = function(expansion_state) 0,
    classify_expansion = classify_rasmussen_expansion,
    propose_expansion = function(expansion_state, trial_step) {
      propose_rasmussen_expansion_alpha(
        expansion_state$initial_step,
        trial_step,
        expansion_factor,
        interior_fraction
      )
    },
    initialize_zoom = new_rasmussen_zoom,
    propose_zoom = function(zoom_state, initial_step) {
      list(
        alpha = propose_rasmussen_zoom_alpha(
          zoom_state,
          initial_step,
          interior_fraction
        ),
        state = zoom_state
      )
    },
    zoom_recovery_lower_bound = function(zoom_state) {
      zoom_state$lower_step$alpha
    },
    zoom_bounds = function(zoom_state) {
      c(zoom_state$lower_step$alpha, zoom_state$upper_step$alpha)
    },
    process_zoom_trial = function(
      zoom_state,
      trial_step,
      initial_step,
      condition_policy,
      direction_scale
    ) {
      process_rasmussen_zoom_trial(
        zoom_state,
        trial_step,
        initial_step,
        condition_policy,
        relative_interval_tolerance
      )
    }
  )
}

classify_rasmussen_expansion <- function(
  expansion_state,
  trial_step,
  condition_policy
) {
  initial_step <- expansion_state$initial_step
  has_bracket <- trial_step$d >
    condition_policy$curvature_constant * initial_step$d ||
    !condition_policy$armijo(initial_step, trial_step)

  list(
    accepted = has_bracket &&
      condition_policy$wolfe(initial_step, trial_step),
    bracket = if (has_bracket) {
      list(initial_step, trial_step)
    } else {
      NULL
    }
  )
}

propose_rasmussen_expansion_alpha <- function(
  initial_step,
  trial_step,
  expansion_factor,
  interior_fraction
) {
  proposed_alpha <- propose_rasmussen_expansion_cubic_alpha(
    initial_step,
    trial_step
  )
  minimum_alpha <- trial_step$alpha +
    interior_fraction * (trial_step$alpha - initial_step$alpha)
  maximum_alpha <- trial_step$alpha * expansion_factor

  if (!isTRUE(is.finite(proposed_alpha)) || proposed_alpha < 0) {
    return(maximum_alpha)
  }
  min(max(proposed_alpha, minimum_alpha), maximum_alpha)
}

propose_rasmussen_expansion_cubic_alpha <- function(
  first_step,
  second_step
) {
  alpha_difference <- second_step$alpha - first_step$alpha
  if (!isTRUE(is.finite(alpha_difference)) || alpha_difference == 0) {
    return(NA_real_)
  }
  linear_coefficient <- 6 *
    (first_step$f - second_step$f) +
    3 * (second_step$d + first_step$d) * alpha_difference
  quadratic_coefficient <- 3 *
    (second_step$f - first_step$f) -
    (2 * first_step$d + second_step$d) * alpha_difference
  discriminant <- quadratic_coefficient^2 -
    linear_coefficient * first_step$d * alpha_difference
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return(NA_real_)
  }

  denominator <- quadratic_coefficient + sqrt(discriminant)
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(NA_real_)
  }
  proposed_alpha <- first_step$alpha -
    first_step$d * alpha_difference^2 / denominator
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(NA_real_)
  }
  proposed_alpha
}

propose_rasmussen_zoom_cubic_alpha <- function(first_step, second_step) {
  cubic_interpolate(
    first_step$alpha,
    first_step$f,
    first_step$d,
    second_step$alpha,
    second_step$f,
    second_step$d
  )
}

new_rasmussen_zoom <- function(bracket) {
  if (bracket[[1L]]$alpha <= bracket[[2L]]$alpha) {
    lower_step <- bracket[[1L]]
    upper_step <- bracket[[2L]]
  } else {
    lower_step <- bracket[[2L]]
    upper_step <- bracket[[1L]]
  }
  list(lower_step = lower_step, upper_step = upper_step)
}

propose_rasmussen_zoom_alpha <- function(
  zoom_state,
  initial_step,
  interior_fraction
) {
  proposed_alpha <- if (zoom_state$upper_step$f > initial_step$f) {
    propose_quadratic_alpha(zoom_state$lower_step, zoom_state$upper_step)
  } else {
    propose_rasmussen_zoom_cubic_alpha(
      zoom_state$lower_step,
      zoom_state$upper_step
    )
  }
  safeguard_rasmussen_zoom_alpha(
    proposed_alpha,
    zoom_state$lower_step$alpha,
    zoom_state$upper_step$alpha,
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
  trial_step,
  initial_step,
  condition_policy
) {
  if (
    trial_step$d > 0 ||
      !condition_policy$armijo(initial_step, trial_step)
  ) {
    zoom_state$upper_step <- trial_step
  } else {
    zoom_state$lower_step <- trial_step
  }
  zoom_state
}

process_rasmussen_zoom_trial <- function(
  zoom_state,
  trial_step,
  initial_step,
  condition_policy,
  relative_interval_tolerance
) {
  progress_stalled <- rasmussen_zoom_has_stalled(
    zoom_state,
    trial_step,
    relative_interval_tolerance
  )
  list(
    state = update_rasmussen_zoom(
      zoom_state,
      trial_step,
      initial_step,
      condition_policy
    ),
    progress_stalled = progress_stalled
  )
}

rasmussen_zoom_has_stalled <- function(
  zoom_state,
  trial_step,
  relative_interval_tolerance
) {
  abs(zoom_state$upper_step$alpha - zoom_state$lower_step$alpha) <
    relative_interval_tolerance * trial_step$alpha
}
