# Step Size ---------------------------------------------------------------


# Constructor -------------------------------------------------------------

make_step_size <- function(sub_stage) {
  sub_stage$type <- "step_size"
  if (!is.null(sub_stage$init)) {
    attr(sub_stage$init, 'event') <- 'init'
    attr(sub_stage$init, 'name') <- paste0(sub_stage$type,' init')
  }
  if (!is.null(sub_stage$calculate)) {
    attr(sub_stage$calculate, 'event') <- sub_stage$type
    attr(sub_stage$calculate, 'name') <- paste0(sub_stage$type, ' calculate')
  }
  if (!is.null(sub_stage$after_step)) {
    attr(sub_stage$after_step, 'event') <- 'after step'
    attr(sub_stage$after_step, 'name') <- paste0(sub_stage$type, ' after step')
  }

  sub_stage
}

# Constant ----------------------------------------------------------------


constant_step_size <- function(value = 1) {
  make_step_size(
    list(
      name = "constant",
      calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
        message("Calculating constant step size")
        list(stage = stage)
      },
      value = value
    )
  )
}



# Bold Driver -------------------------------------------------------------



# Wants: gradient (optional)
bold_driver <- function(inc_mult = 1.1, dec_mult = 0.5,
                inc_fn = partial(`*`, inc_mult),
                dec_fn = partial(`*`, dec_mult),
                init_step_size = 1,
                min_step_size = sqrt(.Machine$double.eps),
                max_step_size = NULL) {
  make_step_size(list(
    name = "bold_driver",
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      message("Initializing bold driver for ", stage$type)

      if (!is_first_stage(opt, stage)) {
        # Bold driver requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn eager updating
        message("bold driver for ", stage$type, ": setting stage updating to eager")
        opt$eager_update <- TRUE
      }
      sub_stage$value <- sub_stage$init_value
      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      pm <- stage$direction$value

      # Optionally use the gradient if it's available to give up early
      # if we're not going downhill
      if (stage == "gradient_descent"
          && has_gr_curr(opt, iter)
          && dot(opt$cache$gr_curr, pm) > 0) {
        message(stage$type, " ", sub_stage$name, ": direction is not descent, setting to min value")
        sub_stage$value <- sub_stage$min_value
      }
      else {
        message(stage$type, " ", sub_stage$name, " calculating")

        if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
          message(stage$type, " ", sub_stage$name, ": fetching fn_curr from cache")
          f0 <- opt$cache$fn_curr
        }
        else {
          message(stage$type, " ", sub_stage$name, ": calculating f0")
          opt <- calc_fn(opt, par, fn)
          f0 <- opt$fn
        }

        alpha <- sub_stage$value
        para <- par + pm * alpha
        opt <- calc_fn(opt, para, fn)
        while ((!is.finite(opt$fn) || opt$fn > f0)
               && alpha > sub_stage$min_value) {
          alpha <- sclamp(sub_stage$dec_fn(alpha),
                          min = sub_stage$min_value,
                          max = sub_stage$max_value)

          para <- par + pm * alpha
          opt <- calc_fn(opt, para, fn)

        }
        sub_stage$value <- alpha
        if (!is.finite(opt$fn)) {
          stop(stage$type, " ", sub_stage$name, " non finite cost found at iter ", iter)
        }
        message(stage$type, " ", sub_stage$name,
                " step size = ", formatC(sub_stage$value),
                " cost = ", formatC(opt$fn))

        if (is_last_stage(opt, stage)) {
          message(stage$type, " ", sub_stage$name, " setting fn_step for iter ", iter)
          opt <- set_fn_new(opt, opt$fn, iter)
        }
      }
      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
      message(stage$type, " ", sub_stage$name, " after step")

      alpha_old <- sub_stage$value
      # increase the step size for the next step
      alpha_new <- sub_stage$inc_fn(alpha_old)
      sub_stage$value <- sclamp(alpha_new,
                                min = sub_stage$min_value,
                                max = sub_stage$max_value)
      message(stage$type, " ", sub_stage$name, ": step size is now = ", formatC(sub_stage$value))

      if (is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
        message(stage$type, " ", sub_stage$name,  ": setting next fn_curr from fn_new for ", iter)
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


# Jacobs ------------------------------------------------------------------



# Wants: gradient, update_old
jacobs <- function(inc_mult = 1.1, dec_mult = 0.5,
                   inc_fn = partial(`*`, inc_mult),
                   dec_fn = partial(`*`, dec_mult),
                   init_gain = 1,
                   min_gain = 0.01,
                   epsilon = 1) {


  make_step_size(list(
    name = "jacob",
    inc_fn = inc_fn,
    dec_fn = dec_fn,
    init_gain = init_gain,
    min_gain = min_gain,
    epsilon = epsilon,
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      message("Initializing Jacobs")
      sub_stage$old_gain <- rep(1, length(par))
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      message("Jacobs calculate")
      gm <- opt$cache$gr_curr
      old_gain <- sub_stage$old_gain
      inc_fn <- sub_stage$inc_fn
      dec_fn <- sub_stage$dec_fn
      update <- opt$cache$update_old
      min_gain <- sub_stage$min_gain

      new_gain <-
        inc_fn(old_gain) * abs(sign(gm) != sign(update)) +
        dec_fn(old_gain) * abs(sign(gm) == sign(update))

      new_gain <- clamp(new_gain, min_gain)

      sub_stage$gain <- new_gain
      sub_stage$value <- new_gain * sub_stage$epsilon
      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
      message("After step Jacobs")
      sub_stage$old_gain <- sub_stage$gain
      list(opt = opt, sub_stage = sub_stage)
    },
    depends = c("gradient", "update_old")
  ))
}



