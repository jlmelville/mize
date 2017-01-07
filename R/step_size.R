# Step Size ---------------------------------------------------------------

1
# Constructor -------------------------------------------------------------

make_step_size <- function(sub_stage) {
  make_sub_stage(sub_stage, 'step_size')
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
bold_driver <- function(inc_mult = 1.1, dec_mult = 0.5,
                inc_fn = partial(`*`, inc_mult),
                dec_fn = partial(`*`, dec_mult),
                init_step_size = 1,
                min_step_size = sqrt(.Machine$double.eps),
                max_step_size = NULL,
                max_fn = Inf) {
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

      # Optionally use the gradient if it's available to give up early
      # if we're not going downhill
      if (stage == "gradient_descent"
          && has_gr_curr(opt, iter)
          && dot(opt$cache$gr_curr, pm) > 0) {
        sub_stage$value <- sub_stage$min_value
      }
      else {
        if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
          f0 <- opt$cache$fn_curr
        }
        else {
          opt <- calc_fn(opt, par, fg$fn)
          f0 <- opt$fn
        }

        max_fn <- max_fn_per_ls(opt, max_fn)
      alpha <- sub_stage$value
        para <- par + pm * alpha
        opt <- calc_fn(opt, para, fg$fn)
        num_steps <- 0
        while ((!is.finite(opt$fn) || opt$fn > f0)
               && alpha > sub_stage$min_value
               && num_steps < max_fn) {
          alpha <- sclamp(sub_stage$dec_fn(alpha),
                          min = sub_stage$min_value,
                          max = sub_stage$max_value)
          para <- par + pm * alpha
          opt <- calc_fn(opt, para, fg$fn)
          num_steps <- num_steps + 1
        }
        sub_stage$value <- alpha
        if (!is.finite(opt$fn)) {
          message(stage$type, " ", sub_stage$name,
                  " non finite cost found at iter ", iter)
          sub_stage$value <- sub_stage$min_value
          return(list(opt = opt, sub_stage = sub_stage))
        }

        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, opt$fn, iter)
        }
      }
      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {
      alpha_old <- sub_stage$value
      # increase the step size for the next step
      if (opt$ok) {
        alpha_new <- sub_stage$inc_fn(alpha_old)
      }
      else {
        alpha_new <- alpha_old
      }

      sub_stage$value <- sclamp(alpha_new,
                                min = sub_stage$min_value,
                                max = sub_stage$max_value)


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


# Backtracking Line Search ------------------------------------------------

# At each stage, starts the line search at init_step_size, and then back tracks
# reducing, the step size by a factor of rho each time, until the Armijo
# sufficient decrease condition is satisfied.
backtracking <- function(rho = 0.5,
                        init_step_size = 1,
                        min_step_size = sqrt(.Machine$double.eps),
                        max_step_size = NULL,
                        c1 = 1e-4,
                        max_fn = Inf) {
  make_step_size(list(
    name = "backtracking",
    init = function(opt, stage, sub_stage, par, fg, iter) {

      if (!is_first_stage(opt, stage)) {
        # Requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn eager updating
        opt$eager_update <- TRUE
      }

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      pm <- stage$direction$value
      sub_stage$value <- sub_stage$init_value

      # Optionally use the gradient if it's available to give up early
      # if we're not going downhill
      if (stage == "gradient_descent"
          && has_gr_curr(opt, iter)
          && dot(opt$cache$gr_curr, pm) > 0) {
        sub_stage$value <- sub_stage$min_value
      }
      else {

        if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
          f0 <- opt$cache$fn_curr
        }
        else {
          opt <- calc_fn(opt, par, fg$fn)
          f0 <- opt$fn
        }

        d0 = dot(opt$cache$gr_curr, pm)

        alpha <- sub_stage$value
        para <- par + pm * alpha
        opt <- calc_fn(opt, para, fg$fn)

        max_fn <- max_fn_per_ls(opt, max_fn)

        while ((!is.finite(opt$fn) || !armijo_ok(f0, d0, alpha, opt$fn, c1))
               && alpha > sub_stage$min_value
               && opt$counts$fn < max_fn) {
          alpha <- sclamp(alpha * rho,
                          min = sub_stage$min_value,
                          max = sub_stage$max_value)

          para <- par + pm * alpha
          opt <- calc_fn(opt, para, fg$fn)
        }
        sub_stage$value <- alpha
        if (!is.finite(opt$fn)) {
          message(stage$type, " ", sub_stage$name,
                  " non finite cost found at iter ", iter)
          sub_stage$value <- sub_stage$min_value
          return(list(opt = opt, sub_stage = sub_stage))
        }

        if (is_last_stage(opt, stage)) {
          opt <- set_fn_new(opt, opt$fn, iter)
        }
      }
      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {

      if (opt$ok && is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
        opt <- set_fn_curr(opt, opt$cache$fn_new, iter + 1)
      }

      list(opt = opt, sub_stage = sub_stage)
    },
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
    max_fn <- min(max_fn,
                  opt$convergence$max_fg - (opt$counts$fn + opt$counts$gr))
  }
  max_fn
}
