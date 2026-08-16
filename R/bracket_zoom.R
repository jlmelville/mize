# Shared bracket-and-zoom control for the supported Rasmussen and Schmidt
# Wolfe profiles.

validate_bracket_zoom_control <- function(
  value,
  name,
  minimum,
  maximum = Inf,
  minimum_open = FALSE,
  maximum_open = FALSE
) {
  validate_line_scalar(value, name)
  below_minimum <- if (minimum_open) value <= minimum else value < minimum
  above_maximum <- if (maximum_open) value >= maximum else value > maximum
  if (is.na(value) || !is.finite(value) || below_minimum || above_maximum) {
    stop(name, " is outside its supported range")
  }
}

run_bracket_zoom <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  direction,
  policy
) {
  evaluator_state <- environment(evaluator)
  initial_step <- evaluator_state$initial_step
  previous_point <- initial_step
  trial_alpha <- initial_alpha
  expansion_iteration <- 0L
  bracket <- NULL

  while (evaluator_state$evaluation_count < evaluator_state$max_evaluations) {
    recovery <- recover_finite_legacy_step(
      evaluator,
      trial_alpha,
      minimum_alpha = policy$expansion_recovery_lower_bound(
        previous_point,
        expansion_iteration
      )
    )
    trial_point <- recovery$step

    if (!recovery$ok) {
      return(list(
        candidate = trial_point,
        termination_reason = "nonfinite_recovery"
      ))
    }

    budget_exhausted <- evaluator_state$evaluation_count >=
      evaluator_state$max_evaluations
    if (
      budget_exhausted &&
        !policy$classify_final_expansion_trial
    ) {
      return(list(
        candidate = trial_point,
        termination_reason = "budget_exhausted"
      ))
    }

    expansion <- policy$classify_expansion(
      initial_step,
      previous_point,
      trial_point,
      expansion_iteration,
      condition_policy
    )
    if (expansion$accepted) {
      return(list(candidate = trial_point, termination_reason = "wolfe"))
    }
    if (budget_exhausted) {
      if (condition_policy$wolfe(initial_step, trial_point)) {
        return(list(candidate = trial_point, termination_reason = "wolfe"))
      }
      return(list(
        candidate = trial_point,
        termination_reason = "budget_exhausted"
      ))
    }
    if (expansion$is_bracketed) {
      bracket <- expansion$bracket
      break
    }

    trial_alpha <- policy$propose_expansion(
      initial_step,
      previous_point,
      trial_point
    )
    previous_point <- trial_point
    expansion_iteration <- expansion_iteration + 1L
  }

  if (is.null(bracket)) {
    return(list(
      candidate = trial_point,
      termination_reason = "budget_exhausted"
    ))
  }

  zoom_state <- policy$initialize_zoom(
    initial_step,
    bracket,
    trial_point
  )
  direction_scale <- max(abs(direction))

  while (evaluator_state$evaluation_count < evaluator_state$max_evaluations) {
    zoom_state <- policy$prepare_zoom(
      zoom_state,
      trial_point,
      initial_step,
      condition_policy
    )
    proposal <- policy$propose_zoom(zoom_state, initial_step)
    zoom_state <- proposal$state

    recovery <- recover_finite_legacy_step(
      evaluator,
      proposal$alpha,
      minimum_alpha = policy$zoom_recovery_lower_bound(zoom_state)
    )
    trial_point <- recovery$step
    if (!recovery$ok) {
      return(list(
        candidate = trial_point,
        termination_reason = "nonfinite_recovery"
      ))
    }

    if (condition_policy$wolfe(initial_step, trial_point)) {
      return(list(candidate = trial_point, termination_reason = "wolfe"))
    }

    if (
      policy$check_progress_before_update &&
        policy$progress_stalled(
          zoom_state,
          trial_point,
          direction_scale
        )
    ) {
      return(list(
        candidate = trial_point,
        termination_reason = "progress_failure"
      ))
    }

    zoom_state <- policy$update_zoom(
      zoom_state,
      trial_point,
      initial_step,
      condition_policy
    )
    if (
      !policy$check_progress_before_update &&
        policy$progress_stalled(
          zoom_state,
          trial_point,
          direction_scale
        )
    ) {
      return(list(
        candidate = trial_point,
        termination_reason = "progress_failure"
      ))
    }
  }

  list(candidate = trial_point, termination_reason = "budget_exhausted")
}

new_rasmussen_bracket_zoom_policy <- function(
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
    interior_fraction = interior_fraction,
    expansion_factor = expansion_factor,
    relative_interval_tolerance = relative_interval_tolerance,
    classify_final_expansion_trial = TRUE,
    check_progress_before_update = TRUE,
    expansion_recovery_lower_bound = function(previous_point, iteration) 0,
    classify_expansion = function(
      initial_point,
      previous_point,
      trial_point,
      iteration,
      condition_policy
    ) {
      is_bracketed <- trial_point$d >
        condition_policy$curvature_constant * initial_point$d ||
        !condition_policy$armijo(initial_point, trial_point)
      list(
        accepted = is_bracketed &&
          condition_policy$wolfe(initial_point, trial_point),
        is_bracketed = is_bracketed,
        bracket = if (is_bracketed) {
          list(initial_point, trial_point)
        } else {
          NULL
        }
      )
    },
    propose_expansion = function(
      initial_point,
      previous_point,
      trial_point
    ) {
      proposed_alpha <- cubic_extrapolate_step(
        initial_point,
        trial_point
      )
      tweak_extrapolation(
        proposed_alpha,
        initial_point$alpha,
        trial_point$alpha,
        expansion_factor,
        interior_fraction
      )
    },
    initialize_zoom = function(initial_point, bracket, trial_point) {
      list(
        lower_point = initial_point,
        upper_point = bracket[[2L]]
      )
    },
    prepare_zoom = function(
      state,
      trial_point,
      initial_point,
      condition_policy
    ) {
      if (
        trial_point$d > 0 ||
          !condition_policy$armijo(initial_point, trial_point)
      ) {
        state$upper_point <- trial_point
      } else {
        state$lower_point <- trial_point
      }
      state
    },
    propose_zoom = function(state, initial_point) {
      proposed_alpha <- if (state$upper_point$f > initial_point$f) {
        quadratic_interpolate_step(state$lower_point, state$upper_point)
      } else {
        cubic_interpolate_step(state$lower_point, state$upper_point)
      }
      proposed_alpha <- tweak_interpolation(
        proposed_alpha,
        state$lower_point$alpha,
        state$upper_point$alpha,
        interior_fraction
      )
      list(alpha = proposed_alpha, state = state)
    },
    zoom_recovery_lower_bound = function(state) {
      min(state$lower_point$alpha, state$upper_point$alpha)
    },
    update_zoom = function(
      state,
      trial_point,
      initial_point,
      condition_policy
    ) {
      state
    },
    progress_stalled = function(state, trial_point, direction_scale) {
      abs(state$upper_point$alpha - state$lower_point$alpha) <
        relative_interval_tolerance * trial_point$alpha
    }
  )
}

new_schmidt_bracket_zoom_policy <- function(
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
    expansion_factor = expansion_factor,
    minimum_expansion_fraction = minimum_expansion_fraction,
    interior_fraction = interior_fraction,
    progress_tolerance = progress_tolerance,
    classify_final_expansion_trial = FALSE,
    check_progress_before_update = FALSE,
    expansion_recovery_lower_bound = function(previous_point, iteration) {
      if (iteration == 0L) 0 else previous_point$alpha
    },
    classify_expansion = function(
      initial_point,
      previous_point,
      trial_point,
      iteration,
      condition_policy
    ) {
      armijo_failed <- !condition_policy$armijo(
        initial_point,
        trial_point
      )
      insufficient_decrease <- iteration > 1L &&
        trial_point$f >= previous_point$f
      accepted <- !armijo_failed &&
        !insufficient_decrease &&
        condition_policy$curvature(initial_point, trial_point)
      is_bracketed <- armijo_failed ||
        insufficient_decrease ||
        (!accepted && trial_point$d >= 0)
      list(
        accepted = accepted,
        is_bracketed = is_bracketed,
        bracket = if (is_bracketed) {
          list(previous_point, trial_point)
        } else if (accepted) {
          list(trial_point)
        } else {
          NULL
        }
      )
    },
    propose_expansion = function(
      initial_point,
      previous_point,
      trial_point
    ) {
      minimum_alpha <- trial_point$alpha +
        minimum_expansion_fraction *
          (trial_point$alpha - previous_point$alpha)
      maximum_alpha <- trial_point$alpha * expansion_factor
      propose_schmidt_cubic_alpha(
        previous_point,
        trial_point,
        lower_alpha = minimum_alpha,
        upper_alpha = maximum_alpha
      )
    },
    initialize_zoom = function(initial_point, bracket, trial_point) {
      list(
        first_point = bracket[[1L]],
        second_point = bracket[[2L]],
        lowest_value_position = NULL,
        other_position = NULL,
        insufficient_progress = FALSE
      )
    },
    prepare_zoom = function(
      state,
      trial_point,
      initial_point,
      condition_policy
    ) {
      state$lowest_value_position <- which.min(c(
        state$first_point$f,
        state$second_point$f
      ))
      state$other_position <- 3L - state$lowest_value_position
      state
    },
    propose_zoom = function(state, initial_point) {
      bracket_alphas <- c(
        state$first_point$alpha,
        state$second_point$alpha
      )
      alpha_min <- min(bracket_alphas)
      alpha_max <- max(bracket_alphas)
      alpha_range <- alpha_max - alpha_min
      proposed_alpha <- propose_schmidt_cubic_alpha(
        state$first_point,
        state$second_point,
        lower_alpha = alpha_min,
        upper_alpha = alpha_max
      )
      if (!is.finite(proposed_alpha)) {
        proposed_alpha <- mean(bracket_alphas)
      }

      if (
        alpha_range > 0 &&
          min(
            alpha_max - proposed_alpha,
            proposed_alpha - alpha_min
          ) /
            alpha_range <
            interior_fraction
      ) {
        outside_interval <- proposed_alpha >= alpha_max ||
          proposed_alpha <= alpha_min
        if (state$insufficient_progress || outside_interval) {
          if (
            abs(proposed_alpha - alpha_max) < abs(proposed_alpha - alpha_min)
          ) {
            proposed_alpha <- alpha_max - interior_fraction * alpha_range
          } else {
            proposed_alpha <- alpha_min + interior_fraction * alpha_range
          }
          state$insufficient_progress <- FALSE
        } else {
          state$insufficient_progress <- TRUE
        }
      } else {
        state$insufficient_progress <- FALSE
      }
      list(alpha = proposed_alpha, state = state)
    },
    zoom_recovery_lower_bound = function(state) {
      min(state$first_point$alpha, state$second_point$alpha)
    },
    update_zoom = function(
      state,
      trial_point,
      initial_point,
      condition_policy
    ) {
      point_names <- c("first_point", "second_point")
      lowest_value_name <- point_names[[state$lowest_value_position]]
      other_name <- point_names[[state$other_position]]
      lowest_value_point <- state[[lowest_value_name]]
      other_point <- state[[other_name]]
      if (
        !condition_policy$armijo(initial_point, trial_point) ||
          trial_point$f >= lowest_value_point$f
      ) {
        state[[other_name]] <- trial_point
      } else {
        if (
          trial_point$d *
            (other_point$alpha - lowest_value_point$alpha) >=
            0
        ) {
          state[[other_name]] <- lowest_value_point
        }
        state[[lowest_value_name]] <- trial_point
      }
      state
    },
    progress_stalled = function(state, trial_point, direction_scale) {
      abs(state$second_point$alpha - state$first_point$alpha) *
        direction_scale <
        progress_tolerance
    }
  )
}
