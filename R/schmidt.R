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
  progress_tolerance = 1e-9
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
    progress_tolerance,
    "progress_tolerance",
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
    initialize_zoom = function(bracket) {
      list(
        first_step = bracket[[1L]],
        second_step = bracket[[2L]],
        insufficient_progress = FALSE
      )
    },
    propose_zoom = function(zoom_state, initial_step) {
      propose_schmidt_zoom_alpha(zoom_state, interior_fraction)
    },
    zoom_recovery_lower_bound = function(zoom_state) {
      min(zoom_state$first_step$alpha, zoom_state$second_step$alpha)
    },
    zoom_bounds = function(zoom_state) {
      sort(c(zoom_state$first_step$alpha, zoom_state$second_step$alpha))
    },
    process_zoom_trial = function(
      zoom_state,
      trial_step,
      initial_step,
      condition_policy,
      direction_scale
    ) {
      zoom_state <- update_schmidt_zoom(
        zoom_state,
        trial_step,
        initial_step,
        condition_policy
      )
      list(
        state = zoom_state,
        progress_stalled = abs(
          zoom_state$second_step$alpha - zoom_state$first_step$alpha
        ) *
          direction_scale <
          progress_tolerance
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
  accepted <- !armijo_failed &&
    !objective_stopped_decreasing &&
    condition_policy$curvature(initial_step, trial_step)
  has_bracket <- armijo_failed ||
    objective_stopped_decreasing ||
    (!accepted && trial_step$d >= 0)

  list(
    accepted = accepted,
    bracket = if (has_bracket) {
      list(expansion_state$previous_step, trial_step)
    } else {
      NULL
    }
  )
}

propose_schmidt_zoom_alpha <- function(zoom_state, interior_fraction) {
  bracket_alphas <- c(
    zoom_state$first_step$alpha,
    zoom_state$second_step$alpha
  )
  lower_alpha <- min(bracket_alphas)
  upper_alpha <- max(bracket_alphas)
  interval_width <- upper_alpha - lower_alpha
  proposed_alpha <- propose_schmidt_cubic_alpha(
    zoom_state$first_step,
    zoom_state$second_step,
    lower_alpha = lower_alpha,
    upper_alpha = upper_alpha
  )
  if (!is.finite(proposed_alpha)) {
    proposed_alpha <- mean(bracket_alphas)
  }

  if (
    interval_width > 0 &&
      min(
        upper_alpha - proposed_alpha,
        proposed_alpha - lower_alpha
      ) /
        interval_width <
        interior_fraction
  ) {
    outside_interval <- proposed_alpha >= upper_alpha ||
      proposed_alpha <= lower_alpha
    if (zoom_state$insufficient_progress || outside_interval) {
      if (
        abs(proposed_alpha - upper_alpha) < abs(proposed_alpha - lower_alpha)
      ) {
        proposed_alpha <- upper_alpha - interior_fraction * interval_width
      } else {
        proposed_alpha <- lower_alpha + interior_fraction * interval_width
      }
      zoom_state$insufficient_progress <- FALSE
    } else {
      zoom_state$insufficient_progress <- TRUE
    }
  } else {
    zoom_state$insufficient_progress <- FALSE
  }

  list(alpha = proposed_alpha, state = zoom_state)
}

update_schmidt_zoom <- function(
  zoom_state,
  trial_step,
  initial_step,
  condition_policy
) {
  step_names <- c("first_step", "second_step")
  lower_value_position <- which.min(c(
    zoom_state$first_step$f,
    zoom_state$second_step$f
  ))
  other_position <- 3L - lower_value_position
  lower_value_name <- step_names[[lower_value_position]]
  other_name <- step_names[[other_position]]
  lower_value_step <- zoom_state[[lower_value_name]]
  other_step <- zoom_state[[other_name]]

  if (
    !condition_policy$armijo(initial_step, trial_step) ||
      trial_step$f >= lower_value_step$f
  ) {
    zoom_state[[other_name]] <- trial_step
  } else {
    if (
      trial_step$d *
        (other_step$alpha - lower_value_step$alpha) >=
        0
    ) {
      zoom_state[[other_name]] <- lower_value_step
    }
    zoom_state[[lower_value_name]] <- trial_step
  }
  zoom_state
}

# Armijo search ------------------------------------------------------------

new_schmidt_armijo_search <- function(
  armijo_constant = 0.05,
  step_down = NULL,
  max_fn = Inf,
  progress_tolerance = 1e-9
) {
  validate_line_scalar(armijo_constant, "armijo_constant")
  validate_line_evaluation_limit(max_fn, "max_fn")
  validate_line_scalar(progress_tolerance, "progress_tolerance")
  if (
    is.na(armijo_constant) ||
      !is.finite(armijo_constant) ||
      armijo_constant < 0 ||
      armijo_constant > 1
  ) {
    stop("armijo_constant must be between zero and one")
  }
  if (
    is.na(progress_tolerance) ||
      !is.finite(progress_tolerance) ||
      progress_tolerance < 0
  ) {
    stop("progress_tolerance must be nonnegative and finite")
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
    max_evaluations <- if (evaluates_gradient) {
      min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    } else {
      min(max_fn, total_max_fn, total_max_fg)
    }
    if (max_evaluations <= 0) {
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
      max_evaluations = max_evaluations,
      direction_scale = norm_inf(pm),
      progress_tolerance = progress_tolerance
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
  progress_tolerance
) {
  evaluates_gradient <- is.null(fixed_reduction_factor)
  function_evaluations <- 0L
  gradient_evaluations <- 0L
  best_decrease <- NULL
  trial_alpha <- initial_alpha
  is_initial_trial <- TRUE

  repeat {
    trial_step <- phi(trial_alpha, calc_gradient = evaluates_gradient)
    function_evaluations <- function_evaluations + 1L
    if (evaluates_gradient) {
      gradient_evaluations <- gradient_evaluations + 1L
    }

    trial_is_usable <- schmidt_armijo_trial_is_usable(trial_step)
    has_strict_decrease <- trial_is_usable &&
      isTRUE(trial_step$f < initial_step$f)
    if (
      has_strict_decrease &&
        (is.null(best_decrease) || trial_step$f < best_decrease$f)
    ) {
      best_decrease <- trial_step
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
        direction_scale * trial_alpha <= progress_tolerance
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

  selected_step <- if (is.null(best_decrease)) initial_step else best_decrease
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

schmidt_armijo_trial_is_usable <- function(step) {
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
    return(propose_schmidt_quadratic_alpha(
      initial_step,
      trial_step,
      lower_alpha = 0,
      upper_alpha = trial_step$alpha
    ))
  }

  propose_schmidt_cubic_alpha(
    initial_step,
    trial_step,
    lower_alpha = 0,
    upper_alpha = trial_step$alpha
  )
}

# A finite value without a usable trial gradient can still define a quadratic
# proposal with the initial slope.
propose_schmidt_quadratic_alpha <- function(
  initial_point,
  trial_point,
  lower_alpha,
  upper_alpha
) {
  proposed_alpha <- propose_quadratic_alpha(initial_point, trial_point)
  min(max(proposed_alpha, lower_alpha), upper_alpha)
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

  cubic_shape <- lower_point$d +
    upper_point$d -
    3 *
      (lower_point$f - upper_point$f) /
      (lower_point$alpha - upper_point$alpha)
  discriminant <- cubic_shape^2 - lower_point$d * upper_point$d
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return(midpoint_alpha)
  }

  cubic_root <- sqrt(discriminant)
  denominator <- upper_point$d - lower_point$d + 2 * cubic_root
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(midpoint_alpha)
  }
  proposed_alpha <- upper_point$alpha -
    (upper_point$alpha - lower_point$alpha) *
      ((upper_point$d + cubic_root - cubic_shape) / denominator)
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(midpoint_alpha)
  }

  min(max(proposed_alpha, lower_alpha), upper_alpha)
}
