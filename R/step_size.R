# Step Size ---------------------------------------------------------------

# Constructor -------------------------------------------------------------

make_step_size <- function(sub_stage) {
  make_sub_stage(sub_stage, "step_size")
}

# Constant ----------------------------------------------------------------

# A constant step size
constant_step_size <- function(value = 1) {
  make_step_size(list(
    name = "constant",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      list(sub_stage = sub_stage)
    },
    value = value
  ))
}


# Bold Driver -------------------------------------------------------------
# Performs a back tracking line search, but rather than use the Armijo
# (sufficient decrease) condition, accepts the first step size that provides
# any reduction in the function. On the next iteration, the first candidate
# step size is a multiple of accepted step size at the previous iteration.
# inc_mult - the accepted step size at the previous time step will be multiplied
#   by this amount to generate the first candidate step size at the next
#   time step.
# dec_mult - the candidate step sizes will be multiplied by this value (and
#   hence should be a value between 0 and 1 exclusive) while looking for an
#   an acceptable step size.
# init_step_size - the initial candidate step size for the first line search.
bold_driver <- function(
  inc_mult = 1.1,
  dec_mult = 0.5,
  inc_fn = partial(`*`, inc_mult),
  dec_fn = partial(`*`, dec_mult),
  init_step_size = 1,
  min_step_size = sqrt(.Machine$double.eps),
  max_step_size = NULL,
  max_fn = Inf
) {
  make_step_size(list(
    name = "bold_driver",
    init = function(opt, stage, sub_stage, par, fg, iter) {
      if (!is_first_stage(opt, stage)) {
        # Bold driver requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn eager updating
        opt$eager_update <- TRUE
      }
      sub_stage$value <- sub_stage$init_value
      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      pm <- stage$direction$value
      sub_stage$alpha_init <- sub_stage$value
      sub_stage$ls_nf <- 0L
      sub_stage$ls_ng <- 0L
      sub_stage$ls_reason <- NULL
      sub_stage$ls_outcome <- NULL
      if (has_gr_curr(opt, iter)) {
        sub_stage$slope_init <- dot(opt$cache$gr_curr, pm)
      } else {
        sub_stage$slope_init <- NULL
      }

      # Optionally use the gradient if it's available to give up early
      # if we're not going downhill
      if (isTRUE(all(pm == 0))) {
        sub_stage$value <- 0
        return(list(opt = opt, sub_stage = sub_stage))
      }
      if (
        stage$type == "gradient_descent" &&
          has_gr_curr(opt, iter) &&
          dot(opt$cache$gr_curr, pm) > 0
      ) {
        sub_stage$value <- 0
        sub_stage$ls_reason <- "non_descent_direction"
        sub_stage$ls_outcome <- "no_step"
        return(list(opt = opt, sub_stage = sub_stage))
      }

      remaining_fn <- max_fn_per_ls(opt, max_fn)
      if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
        f0 <- opt$cache$fn_curr
      } else {
        if (remaining_fn <= 0) {
          sub_stage$value <- 0
          return(list(opt = opt, sub_stage = sub_stage))
        }
        opt <- calc_fn(opt, par, fg$fn)
        if (!is.null(opt$terminate)) {
          sub_stage$value <- 0
          return(list(opt = opt, sub_stage = sub_stage))
        }
        sub_stage$ls_nf <- sub_stage$ls_nf + 1L
        f0 <- opt$fn
        remaining_fn <- remaining_fn - 1
      }

      alpha <- sub_stage$value
      accepted <- FALSE
      last_fn <- NULL
      previous_parameters <- NULL
      termination_reason <- "budget_exhausted"
      while (remaining_fn > 0) {
        para <- project_line_parameters(par, alpha, pm)
        if (line_parameters_have_same_values(para, par)) {
          termination_reason <- "rounding_stagnation"
          break
        }
        if (
          !is.null(previous_parameters) &&
            line_parameters_have_same_values(para, previous_parameters)
        ) {
          alpha_new <- sclamp(
            sub_stage$dec_fn(alpha),
            min = sub_stage$min_value,
            max = sub_stage$max_value
          )
          if (identical(alpha_new, alpha)) {
            termination_reason <- "step_size_floor"
            break
          }
          alpha <- alpha_new
          next
        }
        opt <- calc_fn(opt, para, fg$fn)
        if (!is.null(opt$terminate)) {
          sub_stage$value <- 0
          return(list(opt = opt, sub_stage = sub_stage))
        }
        sub_stage$ls_nf <- sub_stage$ls_nf + 1L
        remaining_fn <- remaining_fn - 1
        last_fn <- opt$fn
        previous_parameters <- para

        if (is.finite(last_fn) && isTRUE(last_fn < f0)) {
          accepted <- TRUE
          termination_reason <- "objective_decrease"
          break
        }
        if (remaining_fn <= 0) {
          termination_reason <- "budget_exhausted"
          break
        }
        if (alpha <= sub_stage$min_value) {
          termination_reason <- if (is.finite(last_fn)) {
            "step_size_floor"
          } else {
            "nonfinite_trial"
          }
          break
        }
        alpha_new <- sclamp(
          sub_stage$dec_fn(alpha),
          min = sub_stage$min_value,
          max = sub_stage$max_value
        )
        if (identical(alpha_new, alpha)) {
          termination_reason <- "step_size_floor"
          break
        }
        alpha <- alpha_new
      }

      if (!accepted) {
        sub_stage$value <- 0
        sub_stage$ls_reason <- termination_reason
        sub_stage$ls_outcome <- "no_step"
        return(list(opt = opt, sub_stage = sub_stage))
      }

      sub_stage$value <- alpha
      sub_stage$ls_reason <- termination_reason
      sub_stage$ls_outcome <- "objective_decrease"
      if (is_last_stage(opt, stage)) {
        opt <- set_fn_new(opt, last_fn, iter)
      }
      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0, update) {
      alpha_old <- sub_stage$value
      # increase the step size for the next step
      if (opt$ok) {
        alpha_new <- sub_stage$inc_fn(alpha_old)
      } else {
        alpha_new <- alpha_old
      }

      sub_stage$value <- sclamp(
        alpha_new,
        min = sub_stage$min_value,
        max = sub_stage$max_value
      )

      if (opt$ok && is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
        opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
      }

      list(opt = opt, sub_stage = sub_stage)
    },
    inc_fn = inc_fn,
    dec_fn = dec_fn,
    init_value = init_step_size,
    min_value = min_step_size,
    max_value = max_step_size
  ))
}

max_fn_per_ls <- function(opt, max_ls_fn = Inf) {
  max_fn <- max_ls_fn
  if (!is.null(opt$convergence$max_fn)) {
    max_fn <- min(max_fn, opt$convergence$max_fn - opt$counts$fn)
  }
  if (!is.null(opt$convergence$max_fg)) {
    max_fn <- min(
      max_fn,
      opt$convergence$max_fg - (opt$counts$fn + opt$counts$gr)
    )
  }
  max_fn
}
