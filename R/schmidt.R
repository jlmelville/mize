# Schmidt line searches implementing methods used by Mark Schmidt's minFunc.
# https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.

# Wolfe search -------------------------------------------------------------

schmidt_core <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  direction,
  method_policy
) {
  run_bracket_zoom(
    evaluator = evaluator,
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    condition_policy = condition_policy,
    direction = direction,
    policy = method_policy
  )
}

# Armijo search ------------------------------------------------------------

new_schmidt_armijo_search <- function(
  armijo_constant = 0.05,
  step_down = NULL,
  max_fn = Inf,
  progress_tolerance = 1e-9
) {
  validate_line_scalar(armijo_constant, "armijo_constant")
  validate_line_scalar(max_fn, "max_fn")
  validate_line_scalar(progress_tolerance, "progress_tolerance")
  if (
    is.na(armijo_constant) ||
      !is.finite(armijo_constant) ||
      armijo_constant < 0 ||
      armijo_constant > 1
  ) {
    stop("armijo_constant must be between zero and one")
  }
  if (is.na(max_fn) || max_fn < 0) {
    stop("max_fn must be nonnegative")
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
  proposed_alpha <- quadratic_interpolate_step(initial_point, trial_point)
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

  cubic_shape <- lower_point$d +
    upper_point$d -
    3 *
      (lower_point$f - upper_point$f) /
      (lower_point$alpha - upper_point$alpha)
  discriminant <- cubic_shape^2 - lower_point$d * upper_point$d
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return((lower_alpha + upper_alpha) / 2)
  }

  cubic_root <- sqrt(discriminant)
  proposed_alpha <- upper_point$alpha -
    (upper_point$alpha - lower_point$alpha) *
      ((upper_point$d + cubic_root - cubic_shape) /
        (upper_point$d - lower_point$d + 2 * cubic_root))

  min(max(proposed_alpha, lower_alpha), upper_alpha)
}
