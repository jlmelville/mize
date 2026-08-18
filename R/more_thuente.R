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
  result <- run_more_thuente_search(
    evaluator = evaluator,
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    condition_policy = condition_policy,
    policy = method_policy
  )

  make_line_search_core_result(
    termination_reason = result$state$termination_reason,
    accepted_point = if (identical(result$state$termination_reason, "wolfe")) {
      result$candidate
    } else {
      NULL
    }
  )
}

make_more_thuente_policy <- function(
  relative_interval_tolerance = .Machine$double.eps,
  contraction_factor = 0.66,
  expansion_factor = 4,
  alpha_min = 0,
  alpha_max = Inf,
  safeguard_cubic = FALSE,
  cubic_interior_fraction = 0.001,
  verbose = FALSE
) {
  validate_line_scalar(
    relative_interval_tolerance,
    "relative_interval_tolerance"
  )
  validate_line_scalar(contraction_factor, "contraction_factor")
  validate_line_scalar(expansion_factor, "expansion_factor")
  validate_line_scalar(alpha_min, "alpha_min")
  validate_line_scalar(alpha_max, "alpha_max")
  validate_line_scalar(cubic_interior_fraction, "cubic_interior_fraction")
  if (is.na(relative_interval_tolerance) || relative_interval_tolerance < 0) {
    stop("relative_interval_tolerance must be nonnegative")
  }
  if (
    is.na(contraction_factor) ||
      contraction_factor <= 0 ||
      contraction_factor >= 1
  ) {
    stop("contraction_factor must be between zero and one")
  }
  if (is.na(expansion_factor) || expansion_factor <= 0) {
    stop("expansion_factor must be positive")
  }
  if (is.na(alpha_min) || alpha_min < 0) {
    stop("alpha_min must be nonnegative")
  }
  if (is.na(alpha_max) || alpha_max < alpha_min) {
    stop("alpha_max must be at least alpha_min")
  }
  if (
    is.na(cubic_interior_fraction) ||
      cubic_interior_fraction < 0 ||
      cubic_interior_fraction >= 1
  ) {
    stop("cubic_interior_fraction must be between zero and one")
  }
  if (
    !is.logical(safeguard_cubic) ||
      length(safeguard_cubic) != 1L ||
      is.na(safeguard_cubic)
  ) {
    stop("safeguard_cubic must be TRUE or FALSE")
  }
  if (
    !is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)
  ) {
    stop("verbose must be TRUE or FALSE")
  }

  structure(
    list(
      relative_interval_tolerance = relative_interval_tolerance,
      contraction_factor = contraction_factor,
      expansion_factor = expansion_factor,
      alpha_min = alpha_min,
      alpha_max = alpha_max,
      safeguard_cubic = safeguard_cubic,
      cubic_interior_fraction = cubic_interior_fraction,
      verbose = verbose
    ),
    class = "more_thuente_policy"
  )
}

initialize_more_thuente_search_state <- function(
  initial_point,
  initial_alpha,
  policy
) {
  interval_width <- policy$alpha_max - policy$alpha_min
  list(
    best_endpoint = initial_point,
    other_endpoint = initial_point,
    trial_point = list(alpha = initial_alpha),
    is_bracketed = FALSE,
    modified_function_stage = TRUE,
    current_interval_width = interval_width,
    previous_interval_width = 2 * interval_width,
    trial_lower_bound = policy$alpha_min,
    trial_upper_bound = policy$alpha_max,
    interval_update_case = "initial",
    termination_reason = NULL
  )
}

run_more_thuente_search <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  policy
) {
  state <- initialize_more_thuente_search_state(
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    policy = policy
  )
  if (initial_point$slope >= 0) {
    state$termination_reason <- "non_descent_direction"
    return(list(candidate = initial_point, state = state))
  }

  armijo_slope <- condition_policy$armijo_constant * initial_point$slope
  evaluator_state <- environment(evaluator)

  repeat {
    if (state$is_bracketed) {
      state$trial_lower_bound <- min(
        state$best_endpoint$alpha,
        state$other_endpoint$alpha
      )
      state$trial_upper_bound <- max(
        state$best_endpoint$alpha,
        state$other_endpoint$alpha
      )
    } else {
      state$trial_lower_bound <- state$best_endpoint$alpha
      state$trial_upper_bound <- state$trial_point$alpha +
        policy$expansion_factor *
          (state$trial_point$alpha - state$best_endpoint$alpha)
    }

    state$trial_point$alpha <- max(
      state$trial_point$alpha,
      policy$alpha_min
    )
    state$trial_point$alpha <- min(
      state$trial_point$alpha,
      policy$alpha_max
    )

    if (policy$verbose) {
      message(
        "Bracket: [",
        formatC(state$trial_lower_bound),
        ", ",
        formatC(state$trial_upper_bound),
        "] alpha = ",
        formatC(state$trial_point$alpha)
      )
    }

    recovered <- recover_finite_line_point(
      evaluate_line = evaluator,
      alpha = state$trial_point$alpha,
      min_alpha = state$trial_lower_bound,
      max_evaluations = Inf
    )
    if (!recovered$succeeded) {
      state$termination_reason <- "nonfinite_recovery"
      return(list(candidate = initial_point, state = state))
    }
    state$trial_point <- recovered$line_point

    state$termination_reason <- classify_more_thuente_termination(
      state = state,
      initial_point = initial_point,
      condition_policy = condition_policy,
      policy = policy,
      evaluation_count = evaluator_state$evaluation_count,
      max_evaluations = evaluator_state$max_evaluations
    )
    if (!is.null(state$termination_reason)) {
      candidate <- select_more_thuente_candidate(state)
      if (policy$verbose) {
        message("alpha = ", formatC(candidate$alpha))
      }
      return(list(candidate = candidate, state = state))
    }

    if (
      state$modified_function_stage &&
        line_point_satisfies_weak_wolfe(
          initial_point,
          state$trial_point,
          condition_policy$armijo_constant,
          min(
            condition_policy$armijo_constant,
            condition_policy$curvature_constant
          )
        )
    ) {
      state$modified_function_stage <- FALSE
    }

    use_modified_function <- state$modified_function_stage &&
      state$trial_point$value <= state$best_endpoint$value &&
      !condition_policy$armijo(initial_point, state$trial_point)

    if (use_modified_function) {
      modified_state <- state
      modified_state$best_endpoint <- modify_more_thuente_point(
        state$best_endpoint,
        armijo_slope
      )
      modified_state$other_endpoint <- modify_more_thuente_point(
        state$other_endpoint,
        armijo_slope
      )
      modified_state$trial_point <- modify_more_thuente_point(
        state$trial_point,
        armijo_slope
      )

      update <- update_more_thuente_interval(modified_state, policy)
      state$best_endpoint <- unmodify_more_thuente_point(
        update$state$best_endpoint,
        armijo_slope
      )
      state$other_endpoint <- unmodify_more_thuente_point(
        update$state$other_endpoint,
        armijo_slope
      )
      state$trial_point$alpha <- update$state$trial_point$alpha
      state$is_bracketed <- update$state$is_bracketed
      state$interval_update_case <- update$classification
    } else {
      update <- update_more_thuente_interval(state, policy)
      state <- update$state
      state$interval_update_case <- update$classification
    }

    if (state$is_bracketed) {
      new_interval_width <- abs(
        state$other_endpoint$alpha - state$best_endpoint$alpha
      )
      if (
        new_interval_width >=
          policy$contraction_factor * state$previous_interval_width
      ) {
        if (policy$verbose) {
          message("Interval did not decrease sufficiently: bisecting")
        }
        state$trial_point$alpha <- state$best_endpoint$alpha +
          0.5 *
            (state$other_endpoint$alpha - state$best_endpoint$alpha)
      }
      state$previous_interval_width <- state$current_interval_width
      state$current_interval_width <- new_interval_width
    }
  }
}

classify_more_thuente_termination <- function(
  state,
  initial_point,
  condition_policy,
  policy,
  evaluation_count,
  max_evaluations
) {
  reason <- NULL
  trial <- state$trial_point

  if (
    (state$is_bracketed &&
      (trial$alpha <= state$trial_lower_bound ||
        trial$alpha >= state$trial_upper_bound)) ||
      identical(state$interval_update_case, "invalid")
  ) {
    reason <- "rounding_stagnation"
  }
  if (
    trial$alpha == policy$alpha_max &&
      condition_policy$armijo(initial_point, trial) &&
      !line_point_satisfies_weak_curvature(
        initial_point,
        trial,
        condition_policy$armijo_constant
      )
  ) {
    reason <- "alpha_max"
  }
  if (
    trial$alpha == policy$alpha_min &&
      (!condition_policy$armijo(initial_point, trial) ||
        line_point_satisfies_weak_curvature(
          initial_point,
          trial,
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
      state$trial_upper_bound - state$trial_lower_bound <=
        policy$relative_interval_tolerance * state$trial_upper_bound
  ) {
    reason <- "relative_interval_tolerance"
  }
  if (condition_policy$wolfe(initial_point, trial)) {
    reason <- "wolfe"
  }

  reason
}

select_more_thuente_candidate <- function(state) {
  if (
    state$termination_reason %in%
      c(
        "budget_exhausted",
        "relative_interval_tolerance",
        "rounding_stagnation"
      )
  ) {
    state$best_endpoint
  } else {
    state$trial_point
  }
}

modify_more_thuente_point <- function(point, armijo_slope) {
  modified_point <- point
  modified_point$value <- point$value - point$alpha * armijo_slope
  modified_point$slope <- point$slope - armijo_slope
  modified_point
}

unmodify_more_thuente_point <- function(point, armijo_slope) {
  point$value <- point$value + point$alpha * armijo_slope
  point$slope <- point$slope + armijo_slope
  point
}

update_more_thuente_interval <- function(state, policy) {
  best <- state$best_endpoint
  other <- state$other_endpoint
  trial <- state$trial_point

  if (
    (state$is_bracketed &&
      (trial$alpha <= min(best$alpha, other$alpha) ||
        trial$alpha >= max(best$alpha, other$alpha))) ||
      best$slope * (trial$alpha - best$alpha) >= 0 ||
      state$trial_upper_bound < state$trial_lower_bound
  ) {
    return(list(state = state, classification = "invalid"))
  }

  opposite_slope <- trial$slope * (best$slope / abs(best$slope)) < 0
  is_bounded <- FALSE

  if (trial$value > best$value) {
    classification <- "higher_trial_value"
    is_bounded <- TRUE
    cubic_proposal <- cubic_interpolate(
      best$alpha,
      best$value,
      best$slope,
      trial$alpha,
      trial$value,
      trial$slope
    )
    quadratic_proposal <- quadratic_interpolate(
      best$alpha,
      best$value,
      best$slope,
      trial$alpha,
      trial$value
    )

    if (!isTRUE(is.finite(cubic_proposal))) {
      next_alpha <- quadratic_proposal
    } else if (
      abs(cubic_proposal - best$alpha) < abs(quadratic_proposal - best$alpha)
    ) {
      next_alpha <- safeguard_more_thuente_cubic(
        cubic_proposal,
        best$alpha,
        other$alpha,
        policy
      )
    } else {
      next_alpha <- cubic_proposal +
        (quadratic_proposal - cubic_proposal) / 2
    }
    state$is_bracketed <- TRUE
  } else if (opposite_slope) {
    classification <- "lower_value_opposite_slope"
    cubic_proposal <- cubic_interpolate(
      best$alpha,
      best$value,
      best$slope,
      trial$alpha,
      trial$value,
      trial$slope
    )
    secant_proposal <- propose_slope_secant_alpha(
      trial$alpha,
      trial$slope,
      best$alpha,
      best$slope
    )

    if (!isTRUE(is.finite(cubic_proposal))) {
      next_alpha <- secant_proposal
    } else if (
      abs(cubic_proposal - trial$alpha) > abs(secant_proposal - trial$alpha)
    ) {
      next_alpha <- safeguard_more_thuente_cubic(
        cubic_proposal,
        best$alpha,
        other$alpha,
        policy
      )
    } else {
      next_alpha <- secant_proposal
    }
    state$is_bracketed <- TRUE
  } else if (abs(trial$slope) < abs(best$slope)) {
    classification <- "lower_value_reduced_slope_magnitude"
    is_bounded <- TRUE

    theta <- 3 *
      (best$value - trial$value) /
      (trial$alpha - best$alpha) +
      best$slope +
      trial$slope
    scale <- norm(rbind(theta, best$slope, trial$slope), "i")
    gamma <- scale *
      sqrt(
        max(
          0,
          (theta / scale)^2 -
            (best$slope / scale) * (trial$slope / scale)
        )
      )
    if (trial$alpha > best$alpha) {
      gamma <- -gamma
    }
    p <- (gamma - trial$slope) + theta
    q <- (gamma + (best$slope - trial$slope)) + gamma
    ratio <- p / q

    if (ratio < 0 && gamma != 0) {
      cubic_proposal <- trial$alpha + ratio * (best$alpha - trial$alpha)
    } else if (trial$alpha > best$alpha) {
      cubic_proposal <- state$trial_upper_bound
    } else {
      cubic_proposal <- state$trial_lower_bound
    }
    secant_proposal <- propose_slope_secant_alpha(
      trial$alpha,
      trial$slope,
      best$alpha,
      best$slope
    )

    if (state$is_bracketed) {
      if (
        abs(trial$alpha - cubic_proposal) < abs(trial$alpha - secant_proposal)
      ) {
        next_alpha <- cubic_proposal
      } else {
        next_alpha <- secant_proposal
      }
    } else if (
      abs(trial$alpha - cubic_proposal) > abs(trial$alpha - secant_proposal)
    ) {
      next_alpha <- cubic_proposal
    } else {
      next_alpha <- secant_proposal
    }
  } else {
    classification <- "lower_value_unreduced_slope_magnitude"
    if (state$is_bracketed) {
      cubic_proposal <- cubic_interpolate(
        other$alpha,
        other$value,
        other$slope,
        trial$alpha,
        trial$value,
        trial$slope
      )
      if (!isTRUE(is.finite(cubic_proposal))) {
        cubic_proposal <- (other$alpha + trial$alpha) / 2
      }
      next_alpha <- safeguard_more_thuente_cubic(
        cubic_proposal,
        best$alpha,
        other$alpha,
        policy
      )
    } else if (trial$alpha > best$alpha) {
      next_alpha <- state$trial_upper_bound
    } else {
      next_alpha <- state$trial_lower_bound
    }
  }

  if (trial$value > best$value) {
    state$other_endpoint <- trial
  } else {
    if (opposite_slope) {
      state$other_endpoint <- best
    }
    state$best_endpoint <- trial
  }

  next_alpha <- min(state$trial_upper_bound, next_alpha)
  next_alpha <- max(state$trial_lower_bound, next_alpha)
  if (state$is_bracketed && is_bounded) {
    weighted_midpoint <- state$best_endpoint$alpha +
      policy$contraction_factor *
        (state$other_endpoint$alpha - state$best_endpoint$alpha)
    if (state$other_endpoint$alpha > state$best_endpoint$alpha) {
      next_alpha <- min(weighted_midpoint, next_alpha)
    } else {
      next_alpha <- max(weighted_midpoint, next_alpha)
    }
  }

  state$trial_point$alpha <- next_alpha
  state$interval_update_case <- classification
  list(state = state, classification = classification)
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
