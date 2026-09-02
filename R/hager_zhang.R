# Hager-Zhang line search
#
# Implements the line search described in:
#
# Hager, W. W., & Zhang, H. (2005).
# A new conjugate gradient method with guaranteed descent and an efficient line
# search. SIAM Journal on Optimization, 16(1), 170-192.
#
# Hager, W. W., & Zhang, H. (2006).
# Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed
# descent. ACM Transactions on Mathematical Software, 32(1), 113-137.

make_hager_zhang_search <- function(
  armijo_constant = curvature_constant / 2,
  curvature_constant = 0.1,
  max_evaluations = Inf,
  strong_curvature = FALSE,
  approximate_armijo = TRUE,
  method_policy = make_hager_zhang_policy()
) {
  make_wolfe_line_search(
    core = run_hager_zhang_search,
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    max_evaluations = max_evaluations,
    approximation_tolerance = method_policy$approximation_tolerance,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    method_policy = method_policy
  )
}

make_hager_zhang_policy <- function(
  approximation_tolerance = 1e-6,
  bisection_fraction = 0.5,
  expansion_factor = 5,
  interval_contraction_factor = 0.66,
  relative_interval_tolerance = 1e-6
) {
  list(
    approximation_tolerance = approximation_tolerance,
    bisection_fraction = bisection_fraction,
    expansion_factor = expansion_factor,
    interval_contraction_factor = interval_contraction_factor,
    relative_interval_tolerance = relative_interval_tolerance
  )
}

run_hager_zhang_search <- function(
  initial_alpha,
  initial_point,
  evaluator,
  condition_policy,
  method_policy,
  search_direction = NULL
) {
  if (initial_point$slope >= 0) {
    return(make_line_search_core_result("non_descent_direction"))
  }

  approximate_decrease_tolerance <-
    method_policy$approximation_tolerance * abs(initial_point$value)
  recovery <- recover_finite_line_point(
    evaluator = evaluator,
    alpha = initial_alpha,
    min_alpha = 0,
    max_evaluations = Inf
  )
  if (!isTRUE(recovery$succeeded)) {
    return(make_line_search_core_result("nonfinite_recovery"))
  }
  trial_point <- recovery$line_point

  if (condition_policy$wolfe(initial_point, trial_point)) {
    return(make_line_search_core_result("wolfe", trial_point))
  }
  if (!hager_zhang_budget_available(evaluator)) {
    return(make_line_search_core_result("budget_exhausted"))
  }

  bracket_result <- find_hager_zhang_bracket(
    trial_point = trial_point,
    initial_point = initial_point,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy,
    search_direction = search_direction
  )
  bracket <- bracket_result$bracket
  if (!is.null(bracket_result$accepted_point)) {
    return(make_line_search_core_result(
      "wolfe",
      bracket_result$accepted_point
    ))
  }
  if (!bracket_result$succeeded) {
    return(make_line_search_core_result(bracket_result$termination_reason))
  }

  if (!hager_zhang_budget_available(evaluator)) {
    return(make_line_search_core_result("budget_exhausted"))
  }

  previous_bracket <- bracket
  repeat {
    proposed_alpha <- propose_hager_zhang_secant_alpha(
      previous_bracket$lower_endpoint,
      previous_bracket$upper_endpoint
    )
    proposed_alpha <- safeguard_zoom_alpha(
      proposed_alpha,
      hager_zhang_bracket_lower_alpha(previous_bracket),
      hager_zhang_bracket_upper_alpha(previous_bracket)
    )
    if (is.null(proposed_alpha)) {
      return(make_line_search_core_result("progress_failure"))
    }

    recovery <- recover_hager_zhang_bracketed_point(
      evaluator = evaluator,
      alpha = proposed_alpha,
      bracket = previous_bracket,
      initial_point = initial_point,
      condition_policy = condition_policy,
      search_direction = search_direction,
      max_evaluations = Inf
    )
    if (isTRUE(recovery$accepted)) {
      return(make_line_search_core_result("wolfe", recovery$line_point))
    }
    if (!is.null(recovery$termination_reason)) {
      return(make_line_search_core_result(recovery$termination_reason))
    }
    if (!isTRUE(recovery$succeeded)) {
      return(make_line_search_core_result("nonfinite_recovery"))
    }
    trial_point <- recovery$line_point

    if (condition_policy$wolfe(initial_point, trial_point)) {
      return(make_line_search_core_result("wolfe", trial_point))
    }
    if (!hager_zhang_budget_available(evaluator)) {
      return(make_line_search_core_result("budget_exhausted"))
    }

    secant_result <- refine_hager_zhang_bracket_with_secants(
      bracket = previous_bracket,
      trial_point = trial_point,
      initial_point = initial_point,
      evaluator = evaluator,
      approximate_decrease_tolerance = approximate_decrease_tolerance,
      condition_policy = condition_policy,
      method_policy = method_policy,
      search_direction = search_direction
    )
    bracket <- secant_result$bracket
    if (!is.null(secant_result$accepted_point)) {
      return(make_line_search_core_result(
        "wolfe",
        secant_result$accepted_point
      ))
    }
    if (!secant_result$succeeded) {
      return(make_line_search_core_result(secant_result$termination_reason))
    }

    if (!hager_zhang_budget_available(evaluator)) {
      return(make_line_search_core_result("budget_exhausted"))
    }

    previous_width <- hager_zhang_bracket_width(previous_bracket)
    current_width <- hager_zhang_bracket_width(bracket)
    if (hager_zhang_bracket_is_small(bracket, method_policy)) {
      return(make_line_search_core_result(
        hager_zhang_small_bracket_termination_reason(
          bracket,
          initial_point,
          search_direction
        )
      ))
    }

    if (
      current_width > method_policy$interval_contraction_factor * previous_width
    ) {
      midpoint_alpha <- safeguard_zoom_alpha(
        NA_real_,
        hager_zhang_bracket_lower_alpha(bracket),
        hager_zhang_bracket_upper_alpha(bracket)
      )
      if (is.null(midpoint_alpha)) {
        return(make_line_search_core_result("progress_failure"))
      }

      recovery <- recover_hager_zhang_bracketed_point(
        evaluator = evaluator,
        alpha = midpoint_alpha,
        bracket = bracket,
        initial_point = initial_point,
        condition_policy = condition_policy,
        search_direction = search_direction,
        max_evaluations = Inf
      )
      if (isTRUE(recovery$accepted)) {
        return(make_line_search_core_result("wolfe", recovery$line_point))
      }
      if (!is.null(recovery$termination_reason)) {
        return(make_line_search_core_result(recovery$termination_reason))
      }
      if (!isTRUE(recovery$succeeded)) {
        return(make_line_search_core_result("nonfinite_recovery"))
      }
      midpoint_point <- recovery$line_point

      update_result <- update_hager_zhang_bracket(
        bracket = bracket,
        trial_point = midpoint_point,
        initial_point = initial_point,
        evaluator = evaluator,
        approximate_decrease_tolerance = approximate_decrease_tolerance,
        condition_policy = condition_policy,
        method_policy = method_policy,
        search_direction = search_direction
      )
      bracket <- update_result$bracket
      if (!is.null(update_result$accepted_point)) {
        return(make_line_search_core_result(
          "wolfe",
          update_result$accepted_point
        ))
      }
      if (!update_result$succeeded) {
        return(make_line_search_core_result(update_result$termination_reason))
      }

      if (!hager_zhang_budget_available(evaluator)) {
        return(make_line_search_core_result("budget_exhausted"))
      }
    }

    previous_bracket <- bracket
  }
}

find_hager_zhang_bracket <- function(
  trial_point,
  initial_point,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy,
  search_direction = NULL
) {
  previous_point <- initial_point
  penultimate_point <- initial_point

  repeat {
    if (isTRUE(condition_policy$wolfe(initial_point, trial_point))) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(previous_point, trial_point),
        TRUE,
        accepted_point = trial_point
      ))
    }
    if (trial_point$slope >= 0) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(previous_point, trial_point),
        TRUE
      ))
    }
    if (
      trial_point$value > initial_point$value + approximate_decrease_tolerance
    ) {
      if (!hager_zhang_budget_available(evaluator)) {
        return(make_hager_zhang_bracket_result(
          make_hager_zhang_bracket(penultimate_point, previous_point),
          FALSE,
          "budget_exhausted"
        ))
      }
      return(bisect_hager_zhang_bracket(
        bracket = make_hager_zhang_bracket(initial_point, trial_point),
        initial_point = initial_point,
        evaluator = evaluator,
        approximate_decrease_tolerance = approximate_decrease_tolerance,
        condition_policy = condition_policy,
        method_policy = method_policy,
        search_direction = search_direction
      ))
    }

    if (!hager_zhang_budget_available(evaluator)) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(previous_point, trial_point),
        FALSE,
        "budget_exhausted"
      ))
    }

    penultimate_point <- previous_point
    previous_point <- trial_point
    proposed_alpha <- safeguard_expansion_alpha(
      trial_point$alpha * method_policy$expansion_factor,
      trial_point$alpha
    )
    if (is.null(proposed_alpha)) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(penultimate_point, previous_point),
        FALSE,
        "progress_failure"
      ))
    }

    recovery <- recover_finite_line_point(
      evaluator = evaluator,
      alpha = proposed_alpha,
      min_alpha = trial_point$alpha,
      max_evaluations = Inf
    )
    if (!isTRUE(recovery$succeeded)) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(penultimate_point, previous_point),
        FALSE,
        "nonfinite_recovery"
      ))
    }
    trial_point <- recovery$line_point
  }
}

update_hager_zhang_bracket <- function(
  bracket,
  trial_point,
  initial_point,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy,
  search_direction = NULL
) {
  if (
    wolfe_trial_point_is_usable(trial_point) &&
      isTRUE(condition_policy$wolfe(initial_point, trial_point))
  ) {
    return(make_hager_zhang_bracket_result(
      bracket,
      TRUE,
      accepted_point = trial_point
    ))
  }
  classification <- classify_hager_zhang_trial(
    bracket,
    trial_point,
    initial_point,
    approximate_decrease_tolerance
  )

  if (identical(classification, "nonfinite_trial")) {
    return(make_hager_zhang_bracket_result(
      bracket,
      FALSE,
      "nonfinite_recovery"
    ))
  }
  if (identical(classification, "outside_bracket")) {
    return(make_hager_zhang_bracket_result(bracket, TRUE))
  }
  if (identical(classification, "upper_endpoint")) {
    return(make_hager_zhang_bracket_result(
      make_hager_zhang_bracket(bracket$lower_endpoint, trial_point),
      TRUE
    ))
  }
  if (identical(classification, "lower_endpoint")) {
    return(make_hager_zhang_bracket_result(
      make_hager_zhang_bracket(trial_point, bracket$upper_endpoint),
      TRUE
    ))
  }

  bisect_hager_zhang_bracket(
    bracket = make_hager_zhang_bracket(bracket$lower_endpoint, trial_point),
    initial_point = initial_point,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy,
    search_direction = search_direction
  )
}

classify_hager_zhang_trial <- function(
  bracket,
  trial_point,
  initial_point,
  approximate_decrease_tolerance
) {
  if (!wolfe_trial_point_is_usable(trial_point)) {
    return("nonfinite_trial")
  }
  if (!hager_zhang_alpha_is_in_bracket(bracket, trial_point$alpha)) {
    return("outside_bracket")
  }
  if (trial_point$slope >= 0) {
    return("upper_endpoint")
  }
  if (
    trial_point$value <= initial_point$value + approximate_decrease_tolerance
  ) {
    return("lower_endpoint")
  }
  "needs_bisection"
}

bisect_hager_zhang_bracket <- function(
  bracket,
  initial_point,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy,
  search_direction = NULL
) {
  current_bracket <- bracket

  repeat {
    if (!hager_zhang_budget_available(evaluator)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        "budget_exhausted"
      ))
    }
    if (hager_zhang_bracket_is_small(current_bracket, method_policy)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        hager_zhang_small_bracket_termination_reason(
          current_bracket,
          initial_point,
          search_direction
        )
      ))
    }

    lower_alpha <- hager_zhang_bracket_lower_alpha(current_bracket)
    upper_alpha <- hager_zhang_bracket_upper_alpha(current_bracket)
    proposed_alpha <-
      (1 - method_policy$bisection_fraction) *
      lower_alpha +
      method_policy$bisection_fraction * upper_alpha
    if (
      !isTRUE(is.finite(proposed_alpha)) ||
        !isTRUE(proposed_alpha > lower_alpha) ||
        !isTRUE(proposed_alpha < upper_alpha)
    ) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        "progress_failure"
      ))
    }

    recovery <- recover_hager_zhang_bracketed_point(
      evaluator = evaluator,
      alpha = proposed_alpha,
      bracket = current_bracket,
      initial_point = initial_point,
      condition_policy = condition_policy,
      search_direction = search_direction,
      max_evaluations = Inf
    )
    if (isTRUE(recovery$accepted)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        TRUE,
        accepted_point = recovery$line_point
      ))
    }
    if (!is.null(recovery$termination_reason)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        recovery$termination_reason
      ))
    }
    if (!isTRUE(recovery$succeeded)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        "nonfinite_recovery"
      ))
    }
    trial_point <- recovery$line_point

    if (isTRUE(condition_policy$wolfe(initial_point, trial_point))) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        TRUE,
        accepted_point = trial_point
      ))
    }
    if (trial_point$slope >= 0) {
      current_bracket$upper_endpoint <- trial_point
      return(make_hager_zhang_bracket_result(current_bracket, TRUE))
    }
    if (
      trial_point$value <= initial_point$value + approximate_decrease_tolerance
    ) {
      current_bracket$lower_endpoint <- trial_point
    } else {
      current_bracket$upper_endpoint <- trial_point
    }
  }
}

refine_hager_zhang_bracket_with_secants <- function(
  bracket,
  trial_point,
  initial_point,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy,
  search_direction = NULL
) {
  first_update <- update_hager_zhang_bracket(
    bracket = bracket,
    trial_point = trial_point,
    initial_point = initial_point,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy,
    search_direction = search_direction
  )
  if (!is.null(first_update$accepted_point) || !first_update$succeeded) {
    return(first_update)
  }
  updated_bracket <- first_update$bracket

  trial_is_upper_endpoint <-
    trial_point$alpha == updated_bracket$upper_endpoint$alpha
  trial_is_lower_endpoint <-
    trial_point$alpha == updated_bracket$lower_endpoint$alpha
  second_secant_alpha <- NULL
  if (trial_is_upper_endpoint) {
    second_secant_alpha <- propose_hager_zhang_secant_alpha(
      bracket$upper_endpoint,
      updated_bracket$upper_endpoint
    )
  } else if (trial_is_lower_endpoint) {
    second_secant_alpha <- propose_hager_zhang_secant_alpha(
      bracket$lower_endpoint,
      updated_bracket$lower_endpoint
    )
  }

  second_secant_is_usable <- !is.null(second_secant_alpha) &&
    hager_zhang_alpha_is_strictly_in_bracket(
      updated_bracket,
      second_secant_alpha
    )
  if (!second_secant_is_usable) {
    return(make_hager_zhang_bracket_result(updated_bracket, TRUE))
  }
  if (!hager_zhang_budget_available(evaluator)) {
    return(make_hager_zhang_bracket_result(
      updated_bracket,
      FALSE,
      "budget_exhausted"
    ))
  }

  recovery <- recover_hager_zhang_bracketed_point(
    evaluator = evaluator,
    alpha = second_secant_alpha,
    bracket = updated_bracket,
    initial_point = initial_point,
    condition_policy = condition_policy,
    search_direction = search_direction,
    max_evaluations = Inf
  )
  if (isTRUE(recovery$accepted)) {
    return(make_hager_zhang_bracket_result(
      updated_bracket,
      TRUE,
      accepted_point = recovery$line_point
    ))
  }
  if (!is.null(recovery$termination_reason)) {
    return(make_hager_zhang_bracket_result(
      updated_bracket,
      FALSE,
      recovery$termination_reason
    ))
  }
  if (!isTRUE(recovery$succeeded)) {
    return(make_hager_zhang_bracket_result(
      updated_bracket,
      FALSE,
      "nonfinite_recovery"
    ))
  }

  second_update <- update_hager_zhang_bracket(
    bracket = updated_bracket,
    trial_point = recovery$line_point,
    initial_point = initial_point,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy,
    search_direction = search_direction
  )
  second_update
}

recover_hager_zhang_bracketed_point <- function(
  evaluator,
  alpha,
  bracket,
  initial_point,
  condition_policy,
  search_direction,
  max_evaluations = Inf
) {
  recovery_lower_bound <- hager_zhang_bracket_lower_alpha(bracket)
  evaluator_state <- environment(evaluator)
  parameter_projection_witness <- if (
    isTRUE(evaluator_state$parameter_projection_witness_initialized)
  ) {
    evaluator_state$parameter_projection_witness
  } else {
    NULL
  }
  proposal_is_definitely_novel <- FALSE
  if (
    !is.null(parameter_projection_witness) &&
      parameter_projection_witness$index <=
        length(bracket$lower_endpoint$parameters) &&
      parameter_projection_witness$index <=
        length(bracket$upper_endpoint$parameters)
  ) {
    witness_index <- parameter_projection_witness$index
    projected_witness <- parameter_projection_witness$initial_parameter +
      alpha * parameter_projection_witness$direction
    proposal_is_definitely_novel <- is.finite(projected_witness) &&
      projected_witness != bracket$lower_endpoint$parameters[[witness_index]] &&
      projected_witness != bracket$upper_endpoint$parameters[[witness_index]]
  }
  if (proposal_is_definitely_novel) {
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
      min_alpha = recovery_lower_bound,
      first_endpoint = bracket$lower_endpoint,
      second_endpoint = bracket$upper_endpoint,
      initial_point = initial_point,
      condition_policy = condition_policy,
      search_direction = search_direction,
      max_evaluations = max_evaluations
    ))
  }

  recover_bracketed_line_point(
    evaluator = evaluator,
    alpha = alpha,
    min_alpha = recovery_lower_bound,
    first_endpoint = bracket$lower_endpoint,
    second_endpoint = bracket$upper_endpoint,
    initial_point = initial_point,
    condition_policy = condition_policy,
    search_direction = search_direction,
    max_evaluations = max_evaluations
  )
}

propose_hager_zhang_secant_alpha <- function(first_point, second_point) {
  denominator <- second_point$slope - first_point$slope
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(NA_real_)
  }
  proposed_alpha <-
    (first_point$alpha *
      second_point$slope -
      second_point$alpha * first_point$slope) /
    denominator
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(NA_real_)
  }
  proposed_alpha
}

hager_zhang_budget_available <- function(evaluator) {
  evaluator_state <- environment(evaluator)
  evaluator_state$evaluation_count < evaluator_state$max_evaluations
}

make_hager_zhang_bracket <- function(lower_endpoint, upper_endpoint) {
  # The lower endpoint has acceptable value and negative slope; the upper
  # endpoint has nonnegative slope.
  list(
    lower_endpoint = lower_endpoint,
    upper_endpoint = upper_endpoint
  )
}

hager_zhang_bracket_lower_alpha <- function(bracket) {
  bracket$lower_endpoint$alpha
}

hager_zhang_bracket_upper_alpha <- function(bracket) {
  bracket$upper_endpoint$alpha
}

hager_zhang_bracket_width <- function(bracket) {
  abs(
    hager_zhang_bracket_upper_alpha(bracket) -
      hager_zhang_bracket_lower_alpha(bracket)
  )
}

hager_zhang_bracket_is_small <- function(bracket, method_policy) {
  hager_zhang_bracket_width(bracket) <=
    method_policy$relative_interval_tolerance *
      hager_zhang_bracket_upper_alpha(bracket)
}

hager_zhang_small_bracket_termination_reason <- function(
  bracket,
  initial_point,
  search_direction
) {
  if (
    bracketed_line_parameters_are_exhausted(
      bracket$lower_endpoint,
      bracket$upper_endpoint,
      initial_point$parameters,
      search_direction
    )
  ) {
    "rounding_stagnation"
  } else {
    "relative_interval_tolerance"
  }
}

hager_zhang_alpha_is_in_bracket <- function(bracket, alpha) {
  isTRUE(is.finite(alpha)) &&
    alpha >= hager_zhang_bracket_lower_alpha(bracket) &&
    alpha <= hager_zhang_bracket_upper_alpha(bracket)
}

hager_zhang_alpha_is_strictly_in_bracket <- function(bracket, alpha) {
  isTRUE(is.finite(alpha)) &&
    alpha > hager_zhang_bracket_lower_alpha(bracket) &&
    alpha < hager_zhang_bracket_upper_alpha(bracket)
}

make_hager_zhang_bracket_result <- function(
  bracket,
  succeeded,
  termination_reason = NULL,
  accepted_point = NULL
) {
  list(
    bracket = bracket,
    succeeded = succeeded,
    termination_reason = termination_reason,
    accepted_point = accepted_point
  )
}
