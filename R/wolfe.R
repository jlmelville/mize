# Functions for line searches

# p62 of Nocedal & Wright defines a "loose" line search as c1 = 1.e-4, c2 = 0.9
# But note that CG and SD methods are not considered suitable for loose line
# search because of the search directions are not well-scaled. c2 = 0.1 is
# suggested for CG on p34. With the Strong Wolfe conditions, reducing c2 makes
# the line search stricter (i.e. forces it closer to a minimum).

# More-Thuente ------------------------------------------------------------
more_thuente_ls <- function(
  c1 = c2 / 2,
  c2 = 0.1,
  max_alpha = Inf,
  max_alpha_mult = Inf,
  initializer = "s",
  initial_step_length = 1,
  try_newton_step = FALSE,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  approx_armijo = FALSE,
  strong_curvature = TRUE,
  safeguard_cubic = FALSE
) {
  max_line_evaluations <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(
    make_wolfe_line_search(
      core = more_thuente_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = max_line_evaluations,
      strong_curvature = strong_curvature,
      approximate_armijo = approx_armijo,
      method_policy = make_more_thuente_policy(
        alpha_max = max_alpha,
        safeguard_cubic = safeguard_cubic
      )
    ),
    name = "more-thuente",
    max_alpha_mult = max_alpha_mult,
    initializer = initializer,
    initial_step_length = initial_step_length,
    try_newton_step = try_newton_step
  )
}


# Rasmussen ---------------------------------------------------------------

rasmussen_ls <- function(
  c1 = c2 / 2,
  c2 = 0.1,
  interior_fraction = 0.1,
  expansion_factor = 3,
  max_alpha_mult = Inf,
  initializer = "s",
  initial_step_length = 1,
  try_newton_step = FALSE,
  eps = .Machine$double.eps,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  strong_curvature = TRUE,
  approx_armijo = FALSE
) {
  max_line_evaluations <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(
    make_rasmussen_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = max_line_evaluations,
      strong_curvature = strong_curvature,
      approximate_armijo = approx_armijo,
      interior_fraction = interior_fraction,
      expansion_factor = expansion_factor
    ),
    name = "rasmussen",
    max_alpha_mult = max_alpha_mult,
    initializer = initializer,
    initial_step_length = initial_step_length,
    try_newton_step = try_newton_step,
    eps = eps
  )
}


# Schmidt (minfunc) -------------------------------------------------------

schmidt_ls <- function(
  c1 = c2 / 2,
  c2 = 0.1,
  max_alpha_mult = Inf,
  initializer = "s",
  initial_step_length = "schmidt",
  try_newton_step = FALSE,
  eps = .Machine$double.eps,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  strong_curvature = TRUE,
  approx_armijo = FALSE
) {
  max_line_evaluations <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(
    make_schmidt_wolfe_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = max_line_evaluations,
      strong_curvature = strong_curvature,
      approximate_armijo = approx_armijo
    ),
    name = "schmidt",
    max_alpha_mult = max_alpha_mult,
    initializer = initializer,
    initial_step_length = initial_step_length,
    try_newton_step = try_newton_step,
    eps = eps
  )
}


schmidt_armijo_ls <- function(
  c1 = 0.005,
  max_alpha_mult = Inf,
  initializer = "s",
  initial_step_length = "schmidt",
  try_newton_step = FALSE,
  step_down = NULL,
  eps = .Machine$double.eps,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf
) {
  max_line_evaluations <- if (is.null(step_down)) {
    min(max_fn, max_gr, floor(max_fg / 2))
  } else {
    min(max_fn, max_fg)
  }

  line_search(
    make_schmidt_armijo_search(
      armijo_constant = c1,
      max_evaluations = max_line_evaluations,
      step_down = step_down
    ),
    name = "schmidt_armijo",
    max_alpha_mult = max_alpha_mult,
    initializer = initializer,
    initial_step_length = initial_step_length,
    try_newton_step = try_newton_step,
    eps = eps
  )
}


# Hager-Zhang -------------------------------------------------------------

hager_zhang_ls <- function(
  c1 = c2 / 2,
  c2 = 0.1,
  max_alpha_mult = Inf,
  initializer = "hz",
  initial_step_length = "hz",
  try_newton_step = FALSE,
  eps = .Machine$double.eps,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf,
  strong_curvature = FALSE,
  approx_armijo = TRUE,
  scale_rescue = FALSE,
  initialization_observer = NULL
) {
  max_line_evaluations <- min(max_fn, max_gr, floor(max_fg / 2))

  line_search(
    make_hager_zhang_search(
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = max_line_evaluations,
      strong_curvature = strong_curvature,
      approximate_armijo = approx_armijo
    ),
    name = "hager-zhang",
    max_alpha_mult = max_alpha_mult,
    initializer = initializer,
    initial_step_length = initial_step_length,
    local_function_evaluation_limit = max_fn,
    local_gradient_evaluation_limit = max_gr,
    local_combined_evaluation_limit = max_fg,
    try_newton_step = try_newton_step,
    eps = eps,
    hager_zhang_scale_rescue = scale_rescue,
    hager_zhang_initialization_observer = initialization_observer
  )
}

# Line Search -------------------------------------------------------------

line_search <- function(
  search_line,
  name,
  initializer = "slope ratio",
  try_newton_step = FALSE,
  initial_step_length = 1,
  local_function_evaluation_limit = Inf,
  local_gradient_evaluation_limit = Inf,
  local_combined_evaluation_limit = Inf,
  max_alpha_mult = Inf,
  eps = .Machine$double.eps,
  hager_zhang_scale_rescue = FALSE,
  hager_zhang_initialization_observer = NULL
) {
  if (!is.numeric(initializer)) {
    initializer <- match.arg(
      tolower(initializer),
      c("slope ratio", "quadratic", "hz", "hager-zhang")
    )
    if (initializer == "hager-zhang") {
      initializer <- "hz"
    }
  }

  if (!is.numeric(initial_step_length)) {
    first_alpha_initializer <- match.arg(
      tolower(initial_step_length),
      c(
        "rasmussen",
        "scipy",
        "schmidt",
        "hz",
        "hager-zhang"
      )
    )
    if (first_alpha_initializer == "hager-zhang") {
      first_alpha_initializer <- "hz"
    }
  } else {
    first_alpha_initializer <- initial_step_length
  }

  make_step_size(list(
    name = name,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      if (!is_first_stage(opt, stage)) {
        # Requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn on eager updating
        # message("Wolfe: setting stage updating to eager")
        opt$eager_update <- TRUE
      }

      sub_stage$value <- NULL
      sub_stage$initial_alpha <- NULL
      sub_stage$previous_slope <- NULL
      sub_stage$previous_value <- NULL
      sub_stage$gradient <- NULL
      sub_stage$alpha_init <- NULL
      sub_stage$slope_init <- NULL
      sub_stage$ls_reason <- NULL
      sub_stage$ls_outcome <- NULL
      sub_stage$ls_nf <- NULL
      sub_stage$ls_ng <- NULL

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$alpha_init <- NULL
      sub_stage$slope_init <- NULL
      sub_stage$ls_reason <- NULL
      sub_stage$ls_outcome <- NULL
      sub_stage$ls_nf <- NULL
      sub_stage$ls_ng <- NULL

      search_direction <- stage$direction$value
      if (norm2(search_direction) < .Machine$double.eps) {
        sub_stage$value <- 0
        if (is_last_stage(opt, stage) && has_fn_curr(opt, iter)) {
          opt <- set_fn_new(opt, opt$cache$fn_curr, iter)
        }
        return(list(opt = opt, sub_stage = sub_stage))
      }

      if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
        initial_value <- opt$cache$fn_curr
      } else {
        opt <- calc_fn(opt, par, fg$fn)
        if (!is.null(opt$terminate)) {
          sub_stage$value <- 0
          return(list(opt = opt, sub_stage = sub_stage))
        }
        initial_value <- opt$fn
      }

      initial_point <- list(
        alpha = 0,
        value = initial_value,
        gradient = get_gr_curr(opt, iter),
        slope = dot(get_gr_curr(opt, iter), search_direction),
        parameters = par
      )
      sub_stage$slope_init <- initial_point$slope

      previous_alpha <- sub_stage$value

      evaluate_line <- make_line_function(
        par,
        fg,
        search_direction,
        calc_gradient_default = TRUE
      )

      remaining_function_evaluations <- local_function_evaluation_limit
      remaining_gradient_evaluations <- local_gradient_evaluation_limit
      remaining_combined_evaluations <- local_combined_evaluation_limit
      if (
        !is.null(opt$convergence$max_fn) && is.finite(opt$convergence$max_fn)
      ) {
        remaining_function_evaluations <- min(
          remaining_function_evaluations,
          opt$convergence$max_fn - opt$counts$fn
        )
      }
      if (
        !is.null(opt$convergence$max_gr) && is.finite(opt$convergence$max_gr)
      ) {
        remaining_gradient_evaluations <- min(
          remaining_gradient_evaluations,
          opt$convergence$max_gr - opt$counts$gr
        )
      }
      if (
        !is.null(opt$convergence$max_fg) && is.finite(opt$convergence$max_fg)
      ) {
        remaining_combined_evaluations <- min(
          remaining_combined_evaluations,
          opt$convergence$max_fg - (opt$counts$fn + opt$counts$gr)
        )
      }

      proposed_initial_alpha <- 0
      probe_function_evaluations <- 0L
      if (is.numeric(initializer)) {
        proposed_initial_alpha <- initializer
      } else if (
        initializer == "slope ratio" &&
          !is.null(sub_stage$previous_slope)
      ) {
        # described on p59 of Nocedal and Wright
        proposed_initial_alpha <- propose_slope_ratio_alpha(
          previous_alpha,
          sub_stage$previous_slope,
          initial_point,
          eps
        )
      } else if (
        initializer == "quadratic" &&
          !is.null(sub_stage$previous_value)
      ) {
        # quadratic interpolation
        proposed_initial_alpha <- propose_quadratic_initial_alpha(
          sub_stage$previous_value,
          initial_point,
          try_newton_step = try_newton_step
        )
      } else if (initializer == "hz" && !is.null(previous_alpha)) {
        probe_is_affordable <-
          remaining_function_evaluations >= 2 &&
          remaining_gradient_evaluations >= 1 &&
          remaining_combined_evaluations >= 3
        scale_estimate <- if (is.null(sub_stage$previous_slope)) {
          NA_real_
        } else {
          propose_slope_ratio_alpha(
            previous_alpha,
            sub_stage$previous_slope,
            initial_point,
            eps
          )
        }
        trial_evaluation_budget <- min(
          remaining_function_evaluations - as.integer(probe_is_affordable),
          remaining_gradient_evaluations,
          floor(
            (remaining_combined_evaluations -
              as.integer(probe_is_affordable)) /
              2
          )
        )
        initial_alpha_result <- propose_next_hager_zhang_alpha(
          evaluate_line,
          previous_alpha,
          initial_point,
          evaluate_quadratic_probe = probe_is_affordable,
          scale_estimate = scale_estimate,
          trial_evaluation_budget = trial_evaluation_budget,
          scale_rescue = hager_zhang_scale_rescue
        )
        proposed_initial_alpha <- initial_alpha_result$alpha
        probe_function_evaluations <-
          initial_alpha_result$function_evaluations
        opt$counts$fn <- opt$counts$fn + probe_function_evaluations
        remaining_function_evaluations <-
          remaining_function_evaluations - probe_function_evaluations
        remaining_combined_evaluations <-
          remaining_combined_evaluations - probe_function_evaluations
        if (is.function(hager_zhang_initialization_observer)) {
          hager_zhang_initialization_observer(list(
            iter = iter,
            previous_alpha = previous_alpha,
            previous_slope = sub_stage$previous_slope,
            current_slope = initial_point$slope,
            unrescued_alpha = initial_alpha_result$unrescued_alpha,
            scale_estimate = initial_alpha_result$scale_estimate,
            trial_evaluation_budget = initial_alpha_result$trial_evaluation_budget,
            required_contractions = initial_alpha_result$required_contractions,
            available_contractions = initial_alpha_result$available_contractions,
            selected_initial_alpha = initial_alpha_result$alpha,
            rescue_applied = initial_alpha_result$rescue_applied,
            probe_function_evaluations = probe_function_evaluations
          ))
        }
      }

      # Prevent the next step initial guess being too large
      if (
        !is.null(previous_alpha) &&
          isTRUE(proposed_initial_alpha / previous_alpha > max_alpha_mult)
      ) {
        proposed_initial_alpha <- previous_alpha * max_alpha_mult
      }
      sub_stage$value <- proposed_initial_alpha

      initializer_is_nonfinite <- !is.null(sub_stage$value) &&
        !isTRUE(is.finite(sub_stage$value))
      if (
        !initializer_is_nonfinite &&
          (is.null(sub_stage$value) || sub_stage$value <= 0)
      ) {
        sub_stage$value <- guess_initial_alpha(
          first_alpha_initializer,
          par,
          initial_point$value,
          initial_point$gradient,
          initial_point$slope,
          try_newton_step
        )
      }
      initializer_is_nonfinite <- !isTRUE(is.finite(sub_stage$value))

      sub_stage$initial_alpha <- sub_stage$value
      sub_stage$previous_slope <- initial_point$slope
      sub_stage$previous_value <- initial_point$value

      if (initializer_is_nonfinite) {
        sub_stage$value <- 0
        sub_stage$initial_alpha <- 0
        sub_stage$ls_reason <- "nonfinite_initial_alpha"
        sub_stage$ls_outcome <- "no_step"
        sub_stage$ls_nf <- probe_function_evaluations
        sub_stage$ls_ng <- 0L
        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, initial_point$value, iter)
          sub_stage$gradient <- initial_point$gradient
          sub_stage$gradient_is_current <- TRUE
        }
      } else {
        sub_stage$alpha_init <- sub_stage$value
        search_result <- search_line(
          evaluate_line = evaluate_line,
          initial_point = initial_point,
          initial_alpha = sub_stage$value,
          remaining_function_evaluations = remaining_function_evaluations,
          remaining_gradient_evaluations = remaining_gradient_evaluations,
          remaining_combined_evaluations = remaining_combined_evaluations,
          search_direction = search_direction
        )
        opt$counts$fn <- opt$counts$fn + search_result$function_evaluations
        opt$counts$gr <- opt$counts$gr + search_result$gradient_evaluations

        sub_stage$ls_reason <- search_result$termination_reason
        sub_stage$ls_outcome <- search_result$outcome
        sub_stage$ls_nf <-
          probe_function_evaluations + search_result$function_evaluations
        sub_stage$ls_ng <- search_result$gradient_evaluations

        sub_stage$gradient_is_current <-
          !is.null(search_result$line_point$gradient)
        sub_stage$value <- search_result$line_point$alpha

        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, search_result$line_point$value, iter)
          if (is.null(search_result$line_point$gradient)) {
            sub_stage$gradient <- rep(sub_stage$eps, length(par))
          } else {
            sub_stage$gradient <- search_result$line_point$gradient
          }
        }
      }

      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0, update) {
      if (opt$ok && is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
        opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
      }

      # Armijo LS does not necessarily calculate gradients
      if (
        opt$ok &&
          is_single_stage(opt) &&
          (is.null(sub_stage$gradient_is_current) ||
            sub_stage$gradient_is_current)
      ) {
        opt <- set_gr_curr(opt, sub_stage$gradient, iter + 1)
      }

      list(opt = opt)
    },
    depends = c("gradient")
  ))
}

make_line_function <- function(
  parameters,
  fg,
  search_direction,
  calc_gradient_default = FALSE
) {
  # LS functions are responsible for updating fn and gr count
  function(alpha, calc_gradient = calc_gradient_default) {
    trial_parameters <- parameters + (alpha * search_direction)

    if (!is.null(fg$fg) && calc_gradient) {
      fg_res <- mize_validate_combined_result(
        fg$fg(trial_parameters),
        length(trial_parameters),
        "fg$fg(par)"
      )
      value <- fg_res[["fn", exact = TRUE]]
      gradient <- fg_res[["gr", exact = TRUE]]

      line_point <- list(
        alpha = alpha,
        value = value,
        gradient = gradient,
        slope = dot(gradient, search_direction),
        parameters = trial_parameters
      )
    } else {
      value <- mize_validate_objective_result(
        fg$fn(trial_parameters),
        "fg$fn(par)"
      )
      if (calc_gradient) {
        gradient <- mize_validate_gradient_result(
          fg$gr(trial_parameters),
          length(trial_parameters),
          "fg$gr(par)"
        )
        line_point <- list(
          alpha = alpha,
          value = value,
          gradient = gradient,
          slope = dot(gradient, search_direction),
          parameters = trial_parameters
        )
      } else {
        line_point <- list(
          alpha = alpha,
          value = value,
          parameters = trial_parameters
        )
      }
    }
    line_point
  }
}

# Ensure Valid Step Size
#
# Given an initial step size, if either the function value or the directional
# derivative is non-finite (NaN or infinite), reduce the step size until
# finite values are found.
#
# @param evaluator Managed line evaluator.
# @param alpha Initial step size.
# @param min_alpha Minimum step size.
# @param max_evaluations Maximum number of function evaluations allowed.
# @return List containing:
#
# * `line_point`: Valid point or the last point evaluated, or NULL if the
#     evaluation allowance is zero.
# * `succeeded`: Whether a usable point was found.
recover_finite_line_point <- function(
  evaluator,
  alpha,
  min_alpha = 0,
  max_evaluations = 20
) {
  evaluator_state <- environment(evaluator)
  final_evaluation_count <- min(
    evaluator_state$max_evaluations,
    evaluator_state$evaluation_count + max_evaluations
  )
  succeeded <- FALSE
  line_point <- NULL
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      !isTRUE(is.finite(alpha)) ||
      !is.numeric(min_alpha) ||
      length(min_alpha) != 1L ||
      !isTRUE(is.finite(min_alpha)) ||
      !isTRUE(alpha >= min_alpha)
  ) {
    return(list(
      line_point = line_point,
      succeeded = succeeded
    ))
  }
  while (
    evaluator_state$evaluation_count < final_evaluation_count &&
      alpha >= min_alpha
  ) {
    line_point <- evaluator(alpha)
    if (wolfe_trial_point_is_usable(line_point)) {
      succeeded <- TRUE
      if (
        !is.null(evaluator_state$initial_point) &&
          line_point$value < evaluator_state$initial_point$value &&
          (is.null(evaluator_state$best_decreasing_point) ||
            line_point$value < evaluator_state$best_decreasing_point$value)
      ) {
        evaluator_state$best_decreasing_point <- line_point
      }
      break
    }
    next_alpha <- min_alpha + (alpha - min_alpha) / 2
    if (
      !isTRUE(is.finite(next_alpha)) ||
        !isTRUE(next_alpha > min_alpha) ||
        !isTRUE(next_alpha < alpha)
    ) {
      break
    }
    alpha <- next_alpha
  }
  list(
    line_point = line_point,
    succeeded = succeeded
  )
}

# Initial Step Length -----------------------------------------------------

# Set the initial step length. If initial_step_length is a numeric scalar,
# then use that as-is. Otherwise, use one of several variations based around
# the only thing we know (the directional derivative)
guess_initial_alpha <- function(
  initializer,
  parameters,
  value,
  gradient,
  slope,
  try_newton_step = FALSE
) {
  if (is.numeric(initializer)) {
    return(initializer)
  }

  proposed_alpha <- switch(
    initializer,
    rasmussen = propose_rasmussen_first_alpha(slope),
    scipy = propose_scipy_first_alpha(gradient, slope),
    schmidt = propose_schmidt_first_alpha(gradient),
    hz = propose_first_hager_zhang_alpha(
      parameters,
      value,
      gradient,
      initial_scale = 0.01
    )
  )

  if (try_newton_step) {
    proposed_alpha <- min(1, 1.01 * proposed_alpha)
  }
  proposed_alpha
}

# Rasmussen scales the first trial by the initial directional derivative.
propose_rasmussen_first_alpha <- function(initial_slope) {
  1 / (1 - initial_slope)
}

# found in _minimize_bfgs in optimize.py with this comment:
# # Sets the initial step guess to dx ~ 1
# actually sets f_old to f0 + 0.5 * ||g||2 then uses f_old in the quadratic
# update formula. Algebraically this is the gradient norm divided by the
# negative directional derivative.
# Assuming steepest descent for initial_point, this can be simplified further to
# the reciprocal square root of the negative slope.
propose_scipy_first_alpha <- function(initial_gradient, initial_slope) {
  -norm2(initial_gradient) / initial_slope
}

# Schmidt scales the first trial by the reciprocal gradient one-norm.
propose_schmidt_first_alpha <- function(initial_gradient) {
  1 / norm1(initial_gradient)
}

# Hager and Zhang's first-search scaling uses the relative size of the
# parameters and gradient, or the objective and squared Euclidean gradient norm
# when the parameter vector is zero.
propose_first_hager_zhang_alpha <- function(
  parameters,
  objective,
  gradient,
  initial_scale = 0.01
) {
  proposed_alpha <- 1
  if (is.null(parameters)) {
    return(proposed_alpha)
  }
  gradient_inf_norm <- norm_inf(gradient)
  if (isTRUE(is.finite(gradient_inf_norm)) && gradient_inf_norm != 0) {
    parameter_inf_norm <- norm_inf(parameters)
    if (isTRUE(is.finite(parameter_inf_norm)) && parameter_inf_norm != 0) {
      proposed_alpha <-
        initial_scale * (parameter_inf_norm / gradient_inf_norm)
    } else if (is_finite_numeric(objective) && objective != 0) {
      gradient_norm_squared <- sqnorm2(gradient)
      if (
        is_finite_numeric(gradient_norm_squared) &&
          gradient_norm_squared != 0
      ) {
        proposed_alpha <-
          initial_scale * (abs(objective) / gradient_norm_squared)
      }
    }
  }
  safeguard_initial_alpha(proposed_alpha)
}

# Next Step Length --------------------------------------------------------

# described on p59 of Nocedal and Wright
# slope ratio method
propose_slope_ratio_alpha <- function(
  previous_alpha,
  previous_slope,
  initial_point,
  eps
) {
  # NB the p vector must be a descent direction or the directional
  # derivative will be positive => a negative initial step size!
  slope_ratio <- previous_slope / (initial_point$slope + eps)
  proposed_alpha <- previous_alpha * slope_ratio
  max(proposed_alpha, eps)
}

# quadratic interpolation
propose_quadratic_initial_alpha <- function(
  previous_value,
  initial_point,
  try_newton_step = FALSE
) {
  proposed_alpha <-
    2 * (initial_point$value - previous_value) / initial_point$slope
  if (try_newton_step) {
    proposed_alpha <- min(1, 1.01 * proposed_alpha)
  }
  max(.Machine$double.eps, proposed_alpha)
}

# On later searches, Hager and Zhang first try a quadratic model from one
# objective-only probe. If that optional model is unavailable or unusable, the
# previous alpha is enlarged by a fixed multiplier.
propose_next_hager_zhang_alpha <- function(
  evaluate_line,
  previous_alpha,
  initial_point,
  quadratic_probe_fraction = 0.1,
  previous_alpha_multiplier = 2,
  evaluate_quadratic_probe = TRUE,
  scale_estimate = NA_real_,
  trial_evaluation_budget = Inf,
  scale_rescue = FALSE
) {
  if (previous_alpha < .Machine$double.eps) {
    return(list(
      alpha = .Machine$double.eps,
      function_evaluations = 0L,
      unrescued_alpha = .Machine$double.eps,
      scale_estimate = scale_estimate,
      trial_evaluation_budget = trial_evaluation_budget,
      required_contractions = NA_real_,
      available_contractions = NA_real_,
      rescue_applied = FALSE
    ))
  }

  proposed_alpha <- safeguard_initial_alpha(
    previous_alpha * previous_alpha_multiplier
  )
  function_evaluations <- 0L
  if (evaluate_quadratic_probe) {
    probe_point <- evaluate_line(
      previous_alpha * quadratic_probe_fraction,
      calc_gradient = FALSE
    )
    function_evaluations <- 1L
    if (
      isTRUE(is.finite(probe_point$value)) &&
        isTRUE(probe_point$value <= initial_point$value)
    ) {
      quadratic_alpha <- propose_quadratic_alpha(initial_point, probe_point)
      quadratic_is_usable <-
        isTRUE(is.finite(initial_point$slope)) &&
        initial_point$slope < 0 &&
        isTRUE(is.finite(quadratic_alpha)) &&
        quadratic_alpha > 0
      if (quadratic_is_usable) {
        proposed_alpha <- safeguard_initial_alpha(quadratic_alpha)
      }
    }
  }

  unrescued_alpha <- proposed_alpha
  required_contractions <- NA_real_
  available_contractions <- NA_real_
  rescue_applied <- FALSE
  scale_is_usable <- is.numeric(scale_estimate) &&
    length(scale_estimate) == 1L &&
    isTRUE(is.finite(scale_estimate)) &&
    scale_estimate > 0
  budget_is_usable <- is.numeric(trial_evaluation_budget) &&
    length(trial_evaluation_budget) == 1L &&
    isTRUE(is.finite(trial_evaluation_budget)) &&
    trial_evaluation_budget >= 1
  if (
    isTRUE(scale_rescue) &&
      scale_is_usable &&
      budget_is_usable &&
      unrescued_alpha > scale_estimate
  ) {
    required_contractions <- max(
      0,
      ceiling(log2(unrescued_alpha) - log2(scale_estimate))
    )
    available_contractions <- max(0, floor(trial_evaluation_budget) - 1)
    if (required_contractions > available_contractions) {
      proposed_alpha <- scale_estimate
      rescue_applied <- TRUE
    }
  }
  list(
    alpha = proposed_alpha,
    function_evaluations = function_evaluations,
    unrescued_alpha = unrescued_alpha,
    scale_estimate = scale_estimate,
    trial_evaluation_budget = trial_evaluation_budget,
    required_contractions = required_contractions,
    available_contractions = available_contractions,
    rescue_applied = rescue_applied
  )
}

safeguard_initial_alpha <- function(alpha) {
  if (isTRUE(alpha == Inf)) {
    return(.Machine$double.xmax)
  }
  if (!isTRUE(is.finite(alpha)) || alpha <= 0) {
    return(.Machine$double.eps)
  }
  max(.Machine$double.eps, alpha)
}

# Line Search Checks -------------------------------------------------------

# Armijo Rule (or Sufficient Decrease Condition)
#
# Line search test.
#
# The sufficient decrease condition is met if the line search step length yields
# a decrease in the function value that is sufficiently large (relative to the
# size of the step).
#
# This test prevents large step sizes that, while representing a function value
# decrease, don't reduce it by very much, which could indicate that the
# function minimum has been stepped over and you're now going back up the slope.
# Also, this condition can always be met by taking a sufficiently small step,
# so line searches involving this condition can always terminate. The downside
# is that you can end up taking very small steps, so it's usual to combine this
# condition with one that encourages larger step sizes.
#
# @param initial_value Objective value at the start of the line search.
# @param initial_slope Directional derivative at the start of the line search.
# @param alpha the step length.
# @param trial_value Objective value at alpha.
# @param c1 the sufficient decrease constant. Should take a value between 0 and
#   1.
# @return `TRUE` if the step `alpha` represents a sufficient decrease.
armijo_condition_is_met <- function(
  initial_value,
  initial_slope,
  alpha,
  trial_value,
  c1
) {
  trial_value <= initial_value + c1 * alpha * initial_slope
}

# Armijo Rule (or Sufficient Decrease Condition)
#
# Line search test.
#
# The sufficient decrease condition is met if the line search step length yields
# a decrease in the function value that is sufficiently large (relative to the
# size of the step).
#
# This test prevents large step sizes that, while representing a function value
# decrease, don't reduce it by very much, which could indicate that the
# function minimum has been stepped over and you're now going back up the slope.
# Also, this condition can always be met by taking a sufficiently small step,
# so line searches involving this condition can always terminate. The downside
# is that you can end up taking very small steps, so it's usual to combine this
# condition with one that encourages larger step sizes.
#
# @param initial_point Line search values at starting point.
# @param trial_point Line search values at a trial point.
# @param c1 the sufficient decrease constant. Should take a value between 0 and
#   1.
# @return `TRUE` if the step represents a sufficient decrease.
line_point_satisfies_armijo <- function(initial_point, trial_point, c1) {
  armijo_condition_is_met(
    initial_point$value,
    initial_point$slope,
    trial_point$alpha,
    trial_point$value,
    c1
  )
}


# Curvature Condition
#
# Line search test.
#
# Ensures that the directional derivative of the line search direction at a
# candidate step size is greater than a specified fraction of the slope of the
# line at the starting point of the search. This condition is used to stop step
# sizes being too small.
#
# In combination with the sufficient decrease condition [armijo_condition_is_met()]
# these conditions make up the Wolfe conditions.
#
# @param initial_slope Directional derivative at the starting point.
# @param trial_slope Directional derivative at the trial point.
# @param c2 Curvature condition constant. Should take a value between `c1`
#  (the constant used in the sufficient decrease condition check) and 1.
# @return `TRUE` if the curvature condition is met.
weak_curvature_condition_is_met <- function(initial_slope, trial_slope, c2) {
  trial_slope >= c2 * initial_slope
}

# Curvature Condition
#
# Line search test.
#
# Ensures that the directional derivative of the line search direction at a
# candidate step size is greater than a specified fraction of the slope of the
# line at the starting point of the search. This condition is used to stop step
# sizes being too small.
#
# In combination with the sufficient decrease condition [armijo_condition_is_met()]
# these conditions make up the Wolfe conditions.
#
# @param initial_point Line search values at starting point.
# @param trial_point Line search values at a trial point.
# @param c2 Curvature condition constant. Should take a value between `c1`
#  (the constant used in the sufficient decrease condition check) and 1.
# @return `TRUE` if the curvature condition is met.
line_point_satisfies_weak_curvature <- function(
  initial_point,
  trial_point,
  c2
) {
  weak_curvature_condition_is_met(initial_point$slope, trial_point$slope, c2)
}
