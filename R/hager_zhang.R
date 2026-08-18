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
  strong_curvature = TRUE,
  approximate_armijo = TRUE,
  method_policy = make_hager_zhang_policy()
) {
  validate_line_evaluation_limit(max_evaluations, "max_evaluations")
  if (!inherits(method_policy, "hager_zhang_policy")) {
    stop("method_policy must be a Hager-Zhang policy")
  }
  condition_policy <- new_line_condition_policy(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    approximation_tolerance = method_policy$approximation_tolerance
  )
  evaluator <- new_line_evaluator(
    function(alpha, calc_gradient = TRUE) NULL,
    max_evaluations = 0
  )
  evaluator_state <- environment(evaluator)

  function(
    phi,
    step0,
    alpha,
    total_max_fn = Inf,
    total_max_gr = Inf,
    total_max_fg = Inf,
    pm = NULL
  ) {
    step0 <- validate_line_step(step0)
    evaluation_limit <- min(
      max_evaluations,
      total_max_fn,
      total_max_gr,
      floor(total_max_fg / 2)
    )
    evaluator_state$phi <- phi
    evaluator_state$max_evaluations <- evaluation_limit
    evaluator_state$initial_step <- step0
    evaluator_state$evaluation_count <- 0L
    evaluator_state$best_decrease <- NULL
    evaluator_state$recovered_nonfinite <- FALSE

    search_result <- if (evaluation_limit <= 0) {
      make_hager_zhang_search_result(step0, "budget_exhausted")
    } else {
      run_hager_zhang_search(
        initial_alpha = alpha,
        initial_step = step0,
        evaluator = evaluator,
        condition_policy = condition_policy,
        method_policy = method_policy
      )
    }
    finalized_result <- finalize_wolfe_line_search_result(
      candidate = search_result$candidate,
      evaluator = evaluator,
      termination_reason = search_result$termination_reason
    )
    finalized_result$termination_reason <- search_result$termination_reason
    finalized_result
  }
}

make_hager_zhang_policy <- function(
  approximation_tolerance = 1e-6,
  bisection_fraction = 0.5,
  expansion_factor = 5,
  interval_contraction_factor = 0.66,
  relative_interval_tolerance = 1e-6
) {
  validate_bracket_zoom_control(
    approximation_tolerance,
    "approximation_tolerance",
    minimum = 0
  )
  validate_bracket_zoom_control(
    bisection_fraction,
    "bisection_fraction",
    minimum = 0,
    maximum = 1,
    minimum_open = TRUE,
    maximum_open = TRUE
  )
  validate_bracket_zoom_control(
    expansion_factor,
    "expansion_factor",
    minimum = 1,
    minimum_open = TRUE
  )
  validate_bracket_zoom_control(
    interval_contraction_factor,
    "interval_contraction_factor",
    minimum = 0,
    maximum = 1,
    minimum_open = TRUE,
    maximum_open = TRUE
  )
  validate_bracket_zoom_control(
    relative_interval_tolerance,
    "relative_interval_tolerance",
    minimum = 0
  )

  structure(
    list(
      approximation_tolerance = approximation_tolerance,
      bisection_fraction = bisection_fraction,
      expansion_factor = expansion_factor,
      interval_contraction_factor = interval_contraction_factor,
      relative_interval_tolerance = relative_interval_tolerance
    ),
    class = "hager_zhang_policy"
  )
}

run_hager_zhang_search <- function(
  initial_alpha,
  initial_step,
  evaluator,
  condition_policy,
  method_policy
) {
  if (!line_step_is_usable(initial_step)) {
    return(make_hager_zhang_search_result(
      initial_step,
      "nonfinite_initial_step"
    ))
  }
  if (initial_step$d >= 0) {
    return(make_hager_zhang_search_result(
      initial_step,
      "non_descent_direction"
    ))
  }

  approximate_decrease_tolerance <-
    method_policy$approximation_tolerance * abs(initial_step$f)
  recovery <- find_finite(
    phi = evaluator,
    alpha = initial_alpha,
    min_alpha = 0,
    max_fn = Inf
  )
  if (!hager_zhang_recovery_succeeded(recovery)) {
    return(make_hager_zhang_search_result(
      initial_step,
      "nonfinite_recovery"
    ))
  }
  trial_step <- recovery$step

  if (condition_policy$wolfe(initial_step, trial_step)) {
    return(make_hager_zhang_search_result(trial_step, "wolfe"))
  }
  if (!hager_zhang_budget_available(evaluator)) {
    return(make_hager_zhang_search_result(
      trial_step,
      "budget_exhausted"
    ))
  }

  bracket_result <- find_hager_zhang_bracket(
    trial_step = trial_step,
    initial_step = initial_step,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy
  )
  bracket <- bracket_result$bracket
  if (!is.null(bracket_result$accepted_step)) {
    return(make_hager_zhang_search_result(
      bracket_result$accepted_step,
      "wolfe"
    ))
  }
  if (!bracket_result$succeeded) {
    return(make_hager_zhang_search_result(
      select_best_hager_zhang_step(bracket),
      bracket_result$termination_reason
    ))
  }

  if (!hager_zhang_budget_available(evaluator)) {
    return(make_hager_zhang_search_result(
      select_best_hager_zhang_step(bracket),
      "budget_exhausted"
    ))
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
      return(make_hager_zhang_search_result(
        select_best_hager_zhang_step(bracket),
        "progress_failure"
      ))
    }

    recovery <- find_finite(
      phi = evaluator,
      alpha = proposed_alpha,
      min_alpha = hager_zhang_bracket_lower_alpha(previous_bracket),
      max_fn = Inf
    )
    if (!hager_zhang_recovery_succeeded(recovery)) {
      return(make_hager_zhang_search_result(
        select_best_hager_zhang_step(bracket),
        "nonfinite_recovery"
      ))
    }
    trial_step <- recovery$step

    if (condition_policy$wolfe(initial_step, trial_step)) {
      return(make_hager_zhang_search_result(trial_step, "wolfe"))
    }
    if (!hager_zhang_budget_available(evaluator)) {
      return(make_hager_zhang_search_result(
        select_best_hager_zhang_step(c(unname(bracket), list(trial_step))),
        "budget_exhausted"
      ))
    }

    secant_result <- refine_hager_zhang_bracket_with_secants(
      bracket = previous_bracket,
      trial_step = trial_step,
      initial_step = initial_step,
      evaluator = evaluator,
      approximate_decrease_tolerance = approximate_decrease_tolerance,
      condition_policy = condition_policy,
      method_policy = method_policy
    )
    bracket <- secant_result$bracket
    if (!is.null(secant_result$accepted_step)) {
      return(make_hager_zhang_search_result(
        secant_result$accepted_step,
        "wolfe"
      ))
    }
    if (!secant_result$succeeded) {
      return(make_hager_zhang_search_result(
        select_best_hager_zhang_step(bracket),
        secant_result$termination_reason
      ))
    }

    if (!hager_zhang_budget_available(evaluator)) {
      return(make_hager_zhang_search_result(
        select_best_hager_zhang_step(bracket),
        "budget_exhausted"
      ))
    }

    previous_width <- hager_zhang_bracket_width(previous_bracket)
    current_width <- hager_zhang_bracket_width(bracket)
    if (hager_zhang_bracket_is_small(bracket, method_policy)) {
      return(make_hager_zhang_search_result(
        select_best_hager_zhang_step(bracket),
        "relative_interval_tolerance"
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
        return(make_hager_zhang_search_result(
          select_best_hager_zhang_step(bracket),
          "progress_failure"
        ))
      }

      recovery <- find_finite(
        phi = evaluator,
        alpha = midpoint_alpha,
        min_alpha = hager_zhang_bracket_lower_alpha(bracket),
        max_fn = Inf
      )
      if (!hager_zhang_recovery_succeeded(recovery)) {
        return(make_hager_zhang_search_result(
          select_best_hager_zhang_step(bracket),
          "nonfinite_recovery"
        ))
      }
      midpoint_step <- recovery$step

      update_result <- update_hager_zhang_bracket(
        bracket = bracket,
        trial_step = midpoint_step,
        initial_step = initial_step,
        evaluator = evaluator,
        approximate_decrease_tolerance = approximate_decrease_tolerance,
        condition_policy = condition_policy,
        method_policy = method_policy
      )
      bracket <- update_result$bracket
      if (!is.null(update_result$accepted_step)) {
        return(make_hager_zhang_search_result(
          update_result$accepted_step,
          "wolfe"
        ))
      }
      if (!update_result$succeeded) {
        return(make_hager_zhang_search_result(
          select_best_hager_zhang_step(bracket),
          update_result$termination_reason
        ))
      }

      if (!hager_zhang_budget_available(evaluator)) {
        return(make_hager_zhang_search_result(
          select_best_hager_zhang_step(bracket),
          "budget_exhausted"
        ))
      }
    }

    previous_bracket <- bracket
  }
}

find_hager_zhang_bracket <- function(
  trial_step,
  initial_step,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy
) {
  previous_trial <- initial_step
  penultimate_trial <- initial_step

  repeat {
    if (isTRUE(condition_policy$wolfe(initial_step, trial_step))) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(previous_trial, trial_step),
        TRUE,
        accepted_step = trial_step
      ))
    }
    if (trial_step$d >= 0) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(previous_trial, trial_step),
        TRUE
      ))
    }
    if (trial_step$f > initial_step$f + approximate_decrease_tolerance) {
      if (!hager_zhang_budget_available(evaluator)) {
        return(make_hager_zhang_bracket_result(
          make_hager_zhang_bracket(penultimate_trial, previous_trial),
          FALSE,
          "budget_exhausted"
        ))
      }
      return(bisect_hager_zhang_bracket(
        bracket = make_hager_zhang_bracket(initial_step, trial_step),
        initial_step = initial_step,
        evaluator = evaluator,
        approximate_decrease_tolerance = approximate_decrease_tolerance,
        condition_policy = condition_policy,
        method_policy = method_policy
      ))
    }

    if (!hager_zhang_budget_available(evaluator)) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(previous_trial, trial_step),
        FALSE,
        "budget_exhausted"
      ))
    }

    penultimate_trial <- previous_trial
    previous_trial <- trial_step
    proposed_alpha <- safeguard_expansion_alpha(
      trial_step$alpha * method_policy$expansion_factor,
      trial_step$alpha
    )
    if (is.null(proposed_alpha)) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(penultimate_trial, previous_trial),
        FALSE,
        "progress_failure"
      ))
    }

    recovery <- find_finite(
      phi = evaluator,
      alpha = proposed_alpha,
      min_alpha = trial_step$alpha,
      max_fn = Inf
    )
    if (!hager_zhang_recovery_succeeded(recovery)) {
      return(make_hager_zhang_bracket_result(
        make_hager_zhang_bracket(penultimate_trial, previous_trial),
        FALSE,
        "nonfinite_recovery"
      ))
    }
    trial_step <- recovery$step
  }
}

update_hager_zhang_bracket <- function(
  bracket,
  trial_step,
  initial_step,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy
) {
  if (
    line_step_is_usable(trial_step) &&
      isTRUE(condition_policy$wolfe(initial_step, trial_step))
  ) {
    return(make_hager_zhang_bracket_result(
      bracket,
      TRUE,
      accepted_step = trial_step
    ))
  }
  classification <- classify_hager_zhang_trial(
    bracket,
    trial_step,
    initial_step,
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
      make_hager_zhang_bracket(bracket$lower_endpoint, trial_step),
      TRUE
    ))
  }
  if (identical(classification, "lower_endpoint")) {
    return(make_hager_zhang_bracket_result(
      make_hager_zhang_bracket(trial_step, bracket$upper_endpoint),
      TRUE
    ))
  }

  bisect_hager_zhang_bracket(
    bracket = make_hager_zhang_bracket(bracket$lower_endpoint, trial_step),
    initial_step = initial_step,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy
  )
}

classify_hager_zhang_trial <- function(
  bracket,
  trial_step,
  initial_step,
  approximate_decrease_tolerance
) {
  if (!line_step_is_usable(trial_step)) {
    return("nonfinite_trial")
  }
  if (!hager_zhang_alpha_is_in_bracket(bracket, trial_step$alpha)) {
    return("outside_bracket")
  }
  if (trial_step$d >= 0) {
    return("upper_endpoint")
  }
  if (trial_step$f <= initial_step$f + approximate_decrease_tolerance) {
    return("lower_endpoint")
  }
  "needs_bisection"
}

bisect_hager_zhang_bracket <- function(
  bracket,
  initial_step,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy
) {
  current_bracket <- bracket

  repeat {
    if (hager_zhang_bracket_is_small(current_bracket, method_policy)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        "relative_interval_tolerance"
      ))
    }
    if (!hager_zhang_budget_available(evaluator)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        "budget_exhausted"
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

    recovery <- find_finite(
      phi = evaluator,
      alpha = proposed_alpha,
      min_alpha = lower_alpha,
      max_fn = Inf
    )
    if (!hager_zhang_recovery_succeeded(recovery)) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        FALSE,
        "nonfinite_recovery"
      ))
    }
    trial_step <- recovery$step

    if (isTRUE(condition_policy$wolfe(initial_step, trial_step))) {
      return(make_hager_zhang_bracket_result(
        current_bracket,
        TRUE,
        accepted_step = trial_step
      ))
    }
    if (trial_step$d >= 0) {
      current_bracket$upper_endpoint <- trial_step
      return(make_hager_zhang_bracket_result(current_bracket, TRUE))
    }
    if (trial_step$f <= initial_step$f + approximate_decrease_tolerance) {
      current_bracket$lower_endpoint <- trial_step
    } else {
      current_bracket$upper_endpoint <- trial_step
    }
  }
}

refine_hager_zhang_bracket_with_secants <- function(
  bracket,
  trial_step,
  initial_step,
  evaluator,
  approximate_decrease_tolerance,
  condition_policy,
  method_policy
) {
  first_update <- update_hager_zhang_bracket(
    bracket = bracket,
    trial_step = trial_step,
    initial_step = initial_step,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy
  )
  if (!is.null(first_update$accepted_step) || !first_update$succeeded) {
    return(first_update)
  }
  updated_bracket <- first_update$bracket

  trial_is_upper_endpoint <-
    trial_step$alpha == updated_bracket$upper_endpoint$alpha
  trial_is_lower_endpoint <-
    trial_step$alpha == updated_bracket$lower_endpoint$alpha
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

  recovery <- find_finite(
    phi = evaluator,
    alpha = second_secant_alpha,
    min_alpha = hager_zhang_bracket_lower_alpha(updated_bracket),
    max_fn = Inf
  )
  if (!hager_zhang_recovery_succeeded(recovery)) {
    return(make_hager_zhang_bracket_result(
      updated_bracket,
      FALSE,
      "nonfinite_recovery"
    ))
  }

  second_update <- update_hager_zhang_bracket(
    bracket = updated_bracket,
    trial_step = recovery$step,
    initial_step = initial_step,
    evaluator = evaluator,
    approximate_decrease_tolerance = approximate_decrease_tolerance,
    condition_policy = condition_policy,
    method_policy = method_policy
  )
  second_update
}

propose_hager_zhang_secant_alpha <- function(first_step, second_step) {
  denominator <- second_step$d - first_step$d
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(NA_real_)
  }
  proposed_alpha <-
    (first_step$alpha * second_step$d - second_step$alpha * first_step$d) /
    denominator
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(NA_real_)
  }
  proposed_alpha
}

select_best_hager_zhang_step <- function(steps) {
  values <- vapply(steps, `[[`, numeric(1), "f")
  steps[[which.min(values)]]
}

hager_zhang_budget_available <- function(evaluator) {
  evaluator_state <- environment(evaluator)
  evaluator_state$evaluation_count < evaluator_state$max_evaluations
}

hager_zhang_recovery_succeeded <- function(recovery) {
  isTRUE(recovery$ok) &&
    !is.null(recovery$step) &&
    line_step_is_usable(recovery$step)
}

make_hager_zhang_bracket <- function(lower_endpoint, upper_endpoint) {
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
  accepted_step = NULL
) {
  list(
    bracket = bracket,
    succeeded = succeeded,
    termination_reason = termination_reason,
    accepted_step = accepted_step
  )
}

make_hager_zhang_search_result <- function(candidate, termination_reason) {
  list(
    candidate = candidate,
    termination_reason = termination_reason
  )
}
