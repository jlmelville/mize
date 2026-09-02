# More'-Thuente Line Search
#
# Adapter from the shared Wolfe substrate to the More-Thuente transition.
# @references
# More, J. J., & Thuente, D. J. (1994).
# Line search algorithms with guaranteed sufficient decrease.
# *ACM Transactions on Mathematical Software (TOMS)*, *20*(3),
# 286-307.
# @seealso This implementation was informed by
#  [Dianne O'Leary's Matlab implementation](https://www.cs.umd.edu/users/oleary/software/)
#  of the More-Thuente method.
more_thuente_core <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  search_direction,
  method_policy
) {
  run_more_thuente_search(
    evaluator = evaluator,
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    condition_policy = condition_policy,
    policy = method_policy,
    search_direction = search_direction
  )
}

make_more_thuente_policy <- function(
  relative_interval_tolerance = .Machine$double.eps,
  contraction_factor = 0.66,
  expansion_factor = 4,
  alpha_min = 0,
  alpha_max = Inf,
  safeguard_cubic = FALSE,
  cubic_interior_fraction = 0.001
) {
  list(
    relative_interval_tolerance = relative_interval_tolerance,
    contraction_factor = contraction_factor,
    expansion_factor = expansion_factor,
    alpha_min = alpha_min,
    alpha_max = alpha_max,
    safeguard_cubic = safeguard_cubic,
    cubic_interior_fraction = cubic_interior_fraction
  )
}

initialize_more_thuente_search_state <- function(
  initial_point,
  initial_alpha,
  policy
) {
  interval_width <- policy$alpha_max - policy$alpha_min
  list(
    reference_endpoint = initial_point,
    other_endpoint = initial_point,
    next_trial_alpha = initial_alpha,
    is_bracketed = FALSE,
    stage = "auxiliary",
    interval_width = interval_width,
    two_trial_reference_width = 2 * interval_width
  )
}

derive_more_thuente_trial_bounds <- function(state, policy) {
  if (state$is_bracketed) {
    return(list(
      lower = min(
        state$reference_endpoint$alpha,
        state$other_endpoint$alpha
      ),
      upper = max(
        state$reference_endpoint$alpha,
        state$other_endpoint$alpha
      )
    ))
  }

  list(
    lower = state$reference_endpoint$alpha,
    upper = state$next_trial_alpha +
      policy$expansion_factor *
        (state$next_trial_alpha - state$reference_endpoint$alpha)
  )
}

run_more_thuente_search <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  policy,
  search_direction = NULL
) {
  state <- initialize_more_thuente_search_state(
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    policy = policy
  )
  if (initial_point$slope >= 0) {
    return(make_line_search_core_result("non_descent_direction"))
  }

  armijo_slope <- condition_policy$armijo_constant * initial_point$slope
  evaluator_state <- environment(evaluator)
  parameter_projection_witness <- if (
    isTRUE(evaluator_state$parameter_projection_witness_initialized)
  ) {
    evaluator_state$parameter_projection_witness
  } else {
    NULL
  }
  previous_transition_is_valid <- TRUE

  repeat {
    trial_bounds <- derive_more_thuente_trial_bounds(state, policy)
    state$next_trial_alpha <- max(state$next_trial_alpha, policy$alpha_min)
    state$next_trial_alpha <- min(state$next_trial_alpha, policy$alpha_max)

    if (state$is_bracketed) {
      proposal_is_definitely_novel <- FALSE
      if (
        !is.null(parameter_projection_witness) &&
          parameter_projection_witness$index <=
            length(state$reference_endpoint$parameters) &&
          parameter_projection_witness$index <=
            length(state$other_endpoint$parameters)
      ) {
        witness_index <- parameter_projection_witness$index
        projected_witness <- parameter_projection_witness$initial_parameter +
          state$next_trial_alpha * parameter_projection_witness$direction
        proposal_is_definitely_novel <- is.finite(projected_witness) &&
          projected_witness !=
            state$reference_endpoint$parameters[[witness_index]] &&
          projected_witness != state$other_endpoint$parameters[[witness_index]]
      }
      if (proposal_is_definitely_novel) {
        recovered <- recover_finite_line_point(
          evaluator,
          state$next_trial_alpha,
          min_alpha = state$next_trial_alpha,
          max_evaluations = 1
        )
        if (recovered$succeeded) {
          recovered$accepted <- FALSE
          recovered$termination_reason <- NULL
        } else {
          recovered <- continue_bracketed_line_recovery(
            evaluator = evaluator,
            recovery = recovered,
            failed_alpha = state$next_trial_alpha,
            min_alpha = trial_bounds$lower,
            first_endpoint = state$reference_endpoint,
            second_endpoint = state$other_endpoint,
            initial_point = initial_point,
            condition_policy = condition_policy,
            search_direction = search_direction,
            max_evaluations = Inf
          )
        }
      } else {
        recovered <- recover_bracketed_line_point(
          evaluator = evaluator,
          alpha = state$next_trial_alpha,
          min_alpha = trial_bounds$lower,
          first_endpoint = state$reference_endpoint,
          second_endpoint = state$other_endpoint,
          initial_point = initial_point,
          condition_policy = condition_policy,
          search_direction = search_direction,
          max_evaluations = Inf
        )
      }
    } else {
      recovered <- recover_finite_line_point(
        evaluator = evaluator,
        alpha = state$next_trial_alpha,
        min_alpha = trial_bounds$lower,
        max_evaluations = Inf
      )
      recovered$accepted <- FALSE
      recovered$termination_reason <- NULL
    }
    if (isTRUE(recovered$accepted)) {
      return(make_line_search_core_result("wolfe", recovered$line_point))
    }
    if (!is.null(recovered$termination_reason)) {
      return(make_line_search_core_result(recovered$termination_reason))
    }
    if (!recovered$succeeded) {
      return(make_line_search_core_result("nonfinite_recovery"))
    }
    trial_point <- recovered$line_point
    state$next_trial_alpha <- trial_point$alpha

    termination_reason <- classify_more_thuente_termination(
      state = state,
      trial_point = trial_point,
      trial_bounds = trial_bounds,
      previous_transition_is_valid = previous_transition_is_valid,
      initial_point = initial_point,
      condition_policy = condition_policy,
      policy = policy,
      evaluation_count = evaluator_state$evaluation_count,
      max_evaluations = evaluator_state$max_evaluations
    )
    if (!is.null(termination_reason)) {
      return(make_line_search_core_result(
        termination_reason,
        accepted_point = if (identical(termination_reason, "wolfe")) {
          trial_point
        } else {
          NULL
        }
      ))
    }

    if (
      identical(state$stage, "auxiliary") &&
        condition_policy$exact_armijo(initial_point, trial_point) &&
        line_point_satisfies_weak_curvature(
          initial_point,
          trial_point,
          min(
            condition_policy$armijo_constant,
            condition_policy$curvature_constant
          )
        )
    ) {
      state$stage <- "objective"
    }

    use_auxiliary_measure <- identical(state$stage, "auxiliary") &&
      trial_point$value <= state$reference_endpoint$value &&
      !condition_policy$selected_armijo(initial_point, trial_point)
    advance <- advance_more_thuente_interval(
      state = state,
      trial_point = trial_point,
      trial_bounds = trial_bounds,
      armijo_slope = if (use_auxiliary_measure) armijo_slope else 0,
      policy = policy
    )
    state <- advance$state
    previous_transition_is_valid <- advance$transition_is_valid
    if (!advance$proposal_is_finite) {
      return(make_line_search_core_result("nonfinite_recovery"))
    }

    if (state$is_bracketed) {
      new_interval_width <- abs(
        state$other_endpoint$alpha - state$reference_endpoint$alpha
      )
      if (
        new_interval_width >=
          policy$contraction_factor * state$two_trial_reference_width
      ) {
        state$next_trial_alpha <- state$reference_endpoint$alpha +
          0.5 *
            (state$other_endpoint$alpha - state$reference_endpoint$alpha)
      }
      state$two_trial_reference_width <- state$interval_width
      state$interval_width <- new_interval_width
    }
  }
}

classify_more_thuente_termination <- function(
  state,
  trial_point,
  trial_bounds,
  previous_transition_is_valid,
  initial_point,
  condition_policy,
  policy,
  evaluation_count,
  max_evaluations
) {
  reason <- NULL

  if (
    (state$is_bracketed &&
      (trial_point$alpha <= trial_bounds$lower ||
        trial_point$alpha >= trial_bounds$upper)) ||
      !previous_transition_is_valid
  ) {
    reason <- "rounding_stagnation"
  }
  if (
    trial_point$alpha == policy$alpha_max &&
      condition_policy$selected_armijo(initial_point, trial_point) &&
      !line_point_satisfies_weak_curvature(
        initial_point,
        trial_point,
        condition_policy$armijo_constant
      )
  ) {
    reason <- "alpha_max"
  }
  if (
    trial_point$alpha == policy$alpha_min &&
      (!condition_policy$selected_armijo(initial_point, trial_point) ||
        line_point_satisfies_weak_curvature(
          initial_point,
          trial_point,
          condition_policy$armijo_constant
        ))
  ) {
    reason <- "alpha_min"
  }
  if (evaluation_count >= max_evaluations) {
    reason <- "budget_exhausted"
  }
  if (
    state$is_bracketed &&
      trial_bounds$upper - trial_bounds$lower <=
        policy$relative_interval_tolerance * trial_bounds$upper
  ) {
    reason <- "relative_interval_tolerance"
  }
  if (condition_policy$wolfe(initial_point, trial_point)) {
    reason <- "wolfe"
  }

  reason
}

measure_more_thuente_point <- function(point, armijo_slope = 0) {
  list(
    value = point$value - point$alpha * armijo_slope,
    slope = point$slope - armijo_slope
  )
}

classify_more_thuente_case <- function(
  state,
  trial_point,
  trial_bounds,
  reference_measure,
  trial_measure
) {
  reference <- state$reference_endpoint
  other <- state$other_endpoint
  if (
    (state$is_bracketed &&
      (trial_point$alpha <= min(reference$alpha, other$alpha) ||
        trial_point$alpha >= max(reference$alpha, other$alpha))) ||
      reference_measure$slope *
        (trial_point$alpha - reference$alpha) >=
        0 ||
      trial_bounds$upper < trial_bounds$lower
  ) {
    return(list(
      classification = "invalid",
      is_bounded = FALSE
    ))
  }

  opposite_slope <- trial_measure$slope *
    (reference_measure$slope / abs(reference_measure$slope)) <
    0
  if (trial_measure$value > reference_measure$value) {
    classification <- "higher_trial_value"
    is_bounded <- TRUE
  } else if (opposite_slope) {
    classification <- "lower_value_opposite_slope"
    is_bounded <- FALSE
  } else if (abs(trial_measure$slope) < abs(reference_measure$slope)) {
    classification <- "lower_value_reduced_slope_magnitude"
    is_bounded <- TRUE
  } else {
    classification <- "lower_value_unreduced_slope_magnitude"
    is_bounded <- FALSE
  }

  list(
    classification = classification,
    is_bounded = is_bounded
  )
}

propose_more_thuente_trial <- function(
  state,
  trial_point,
  trial_bounds,
  reference_measure,
  other_measure,
  trial_measure,
  case,
  policy
) {
  reference <- state$reference_endpoint
  other <- state$other_endpoint
  classification <- case$classification

  if (identical(classification, "higher_trial_value")) {
    cubic_proposal <- cubic_interpolate(
      reference$alpha,
      reference_measure$value,
      reference_measure$slope,
      trial_point$alpha,
      trial_measure$value,
      trial_measure$slope
    )
    quadratic_proposal <- quadratic_interpolate(
      reference$alpha,
      reference_measure$value,
      reference_measure$slope,
      trial_point$alpha,
      trial_measure$value
    )
    cubic_is_valid <- isTRUE(is.finite(cubic_proposal))

    if (!cubic_is_valid) {
      next_alpha <- quadratic_proposal
    } else if (
      abs(cubic_proposal - reference$alpha) <
        abs(quadratic_proposal - reference$alpha)
    ) {
      next_alpha <- safeguard_more_thuente_cubic(
        cubic_proposal,
        reference$alpha,
        other$alpha,
        policy
      )
    } else {
      next_alpha <- cubic_proposal +
        (quadratic_proposal - cubic_proposal) / 2
    }
  } else if (identical(classification, "lower_value_opposite_slope")) {
    cubic_proposal <- cubic_interpolate(
      reference$alpha,
      reference_measure$value,
      reference_measure$slope,
      trial_point$alpha,
      trial_measure$value,
      trial_measure$slope
    )
    secant_proposal <- propose_slope_secant_alpha(
      trial_point$alpha,
      trial_measure$slope,
      reference$alpha,
      reference_measure$slope
    )
    cubic_is_valid <- isTRUE(is.finite(cubic_proposal))

    if (!cubic_is_valid) {
      next_alpha <- secant_proposal
    } else if (
      abs(cubic_proposal - trial_point$alpha) >
        abs(secant_proposal - trial_point$alpha)
    ) {
      next_alpha <- safeguard_more_thuente_cubic(
        cubic_proposal,
        reference$alpha,
        other$alpha,
        policy
      )
    } else {
      next_alpha <- secant_proposal
    }
  } else if (identical(classification, "lower_value_reduced_slope_magnitude")) {
    theta <- 3 *
      (reference_measure$value - trial_measure$value) /
      (trial_point$alpha - reference$alpha) +
      reference_measure$slope +
      trial_measure$slope
    scale <- norm(
      rbind(theta, reference_measure$slope, trial_measure$slope),
      "i"
    )
    gamma <- scale *
      sqrt(
        max(
          0,
          (theta / scale)^2 -
            (reference_measure$slope / scale) *
              (trial_measure$slope / scale)
        )
      )
    if (trial_point$alpha > reference$alpha) {
      gamma <- -gamma
    }
    p <- (gamma - trial_measure$slope) + theta
    q <- (gamma +
      (reference_measure$slope - trial_measure$slope)) +
      gamma
    ratio <- p / q

    if (ratio < 0 && gamma != 0) {
      cubic_proposal <- trial_point$alpha +
        ratio * (reference$alpha - trial_point$alpha)
    } else if (trial_point$alpha > reference$alpha) {
      cubic_proposal <- trial_bounds$upper
    } else {
      cubic_proposal <- trial_bounds$lower
    }
    secant_proposal <- propose_slope_secant_alpha(
      trial_point$alpha,
      trial_measure$slope,
      reference$alpha,
      reference_measure$slope
    )

    if (state$is_bracketed) {
      if (
        abs(trial_point$alpha - cubic_proposal) <
          abs(trial_point$alpha - secant_proposal)
      ) {
        next_alpha <- cubic_proposal
      } else {
        next_alpha <- secant_proposal
      }
    } else if (
      abs(trial_point$alpha - cubic_proposal) >
        abs(trial_point$alpha - secant_proposal)
    ) {
      next_alpha <- cubic_proposal
    } else {
      next_alpha <- secant_proposal
    }
  } else if (state$is_bracketed) {
    cubic_proposal <- cubic_interpolate(
      other$alpha,
      other_measure$value,
      other_measure$slope,
      trial_point$alpha,
      trial_measure$value,
      trial_measure$slope
    )
    cubic_is_valid <- isTRUE(is.finite(cubic_proposal))
    if (!cubic_is_valid) {
      cubic_proposal <- (other$alpha + trial_point$alpha) / 2
    }
    next_alpha <- safeguard_more_thuente_cubic(
      cubic_proposal,
      reference$alpha,
      other$alpha,
      policy
    )
  } else if (trial_point$alpha > reference$alpha) {
    next_alpha <- trial_bounds$upper
  } else {
    next_alpha <- trial_bounds$lower
  }

  list(alpha = next_alpha, is_bounded = case$is_bounded)
}

update_more_thuente_endpoints <- function(
  reference_endpoint,
  other_endpoint,
  trial_point,
  classification
) {
  if (identical(classification, "higher_trial_value")) {
    other_endpoint <- trial_point
  } else if (identical(classification, "lower_value_opposite_slope")) {
    other_endpoint <- reference_endpoint
    reference_endpoint <- trial_point
  } else {
    reference_endpoint <- trial_point
  }

  list(
    reference_endpoint = reference_endpoint,
    other_endpoint = other_endpoint
  )
}

safeguard_more_thuente_trial <- function(
  proposal,
  endpoints,
  is_bracketed,
  trial_bounds,
  policy
) {
  next_alpha <- min(trial_bounds$upper, proposal$alpha)
  next_alpha <- max(trial_bounds$lower, next_alpha)
  if (is_bracketed && proposal$is_bounded) {
    weighted_midpoint <- endpoints$reference_endpoint$alpha +
      policy$contraction_factor *
        (endpoints$other_endpoint$alpha -
          endpoints$reference_endpoint$alpha)
    if (endpoints$other_endpoint$alpha > endpoints$reference_endpoint$alpha) {
      next_alpha <- min(weighted_midpoint, next_alpha)
    } else {
      next_alpha <- max(weighted_midpoint, next_alpha)
    }
  }

  proposal_is_finite <- isTRUE(is.finite(next_alpha))

  list(
    alpha = next_alpha,
    is_finite = proposal_is_finite
  )
}

advance_more_thuente_interval <- function(
  state,
  trial_point,
  trial_bounds,
  armijo_slope,
  policy
) {
  if (armijo_slope == 0) {
    reference_measure <- state$reference_endpoint
    other_measure <- state$other_endpoint
    trial_measure <- trial_point
  } else {
    reference_measure <- measure_more_thuente_point(
      state$reference_endpoint,
      armijo_slope
    )
    other_measure <- measure_more_thuente_point(
      state$other_endpoint,
      armijo_slope
    )
    trial_measure <- measure_more_thuente_point(trial_point, armijo_slope)
  }
  case <- classify_more_thuente_case(
    state,
    trial_point,
    trial_bounds,
    reference_measure,
    trial_measure
  )
  if (identical(case$classification, "invalid")) {
    return(list(
      state = state,
      transition_is_valid = FALSE,
      proposal_is_finite = TRUE
    ))
  }

  proposal <- propose_more_thuente_trial(
    state,
    trial_point,
    trial_bounds,
    reference_measure,
    other_measure,
    trial_measure,
    case,
    policy
  )
  endpoints <- update_more_thuente_endpoints(
    state$reference_endpoint,
    state$other_endpoint,
    trial_point,
    case$classification
  )
  is_bracketed <- state$is_bracketed ||
    case$classification %in%
      c("higher_trial_value", "lower_value_opposite_slope")
  safeguarded <- safeguard_more_thuente_trial(
    proposal,
    endpoints,
    is_bracketed,
    trial_bounds,
    policy
  )

  state$reference_endpoint <- endpoints$reference_endpoint
  state$other_endpoint <- endpoints$other_endpoint
  state$is_bracketed <- is_bracketed
  if (safeguarded$is_finite) {
    state$next_trial_alpha <- safeguarded$alpha
  }

  list(
    state = state,
    transition_is_valid = TRUE,
    proposal_is_finite = safeguarded$is_finite
  )
}

safeguard_more_thuente_cubic <- function(
  alpha,
  first_endpoint_alpha,
  second_endpoint_alpha,
  policy
) {
  if (!policy$safeguard_cubic) {
    return(alpha)
  }
  lower_alpha <- min(first_endpoint_alpha, second_endpoint_alpha)
  interval_width <- abs(first_endpoint_alpha - second_endpoint_alpha)
  max(
    lower_alpha + policy$cubic_interior_fraction * interval_width,
    alpha
  )
}
