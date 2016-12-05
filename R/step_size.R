# Step Size ---------------------------------------------------------------


# Constructor -------------------------------------------------------------

make_step_size <- function(sub_stage) {
  make_sub_stage(sub_stage, 'step_size')
}

# Constant ----------------------------------------------------------------


constant_step_size <- function(value = 1) {
  make_step_size(list(
      name = "constant",
      calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
        #message("Constant step size: ", formatC(sub_stage$value))
        list()
      },
      value = value
    ))
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
      #message("Initializing bold driver for ", stage$type)

      if (!is_first_stage(opt, stage)) {
        # Bold driver requires knowing f at the current location
        # If this step size is part of any stage other than the first
        # we have to turn eager updating
        #message("bold driver for ", stage$type, ": setting stage updating to eager")
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
        #message(stage$type, " ", sub_stage$name, ": direction is not descent, setting to min value")
        sub_stage$value <- sub_stage$min_value
      }
      else {
        #message(stage$type, " ", sub_stage$name, " calculating")

        if (is_first_stage(opt, stage) && has_fn_curr(opt, iter)) {
          #message(stage$type, " ", sub_stage$name, ": fetching fn_curr from cache")
          f0 <- opt$cache$fn_curr
        }
        else {
          #message(stage$type, " ", sub_stage$name, ": calculating f0")
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
          #message(stage$type, " ", sub_stage$name, " calculating cost for candidate step size")
          opt <- calc_fn(opt, para, fn)
        }
        sub_stage$value <- alpha
        if (!is.finite(opt$fn)) {
          stop(stage$type, " ", sub_stage$name, " non finite cost found at iter ", iter)
        }
        #message(stage$type, " ", sub_stage$name,
         #       " step size = ", formatC(sub_stage$value),
          #      " cost = ", formatC(opt$fn))

        if (is_last_stage(opt, stage)) {
          #message(stage$type, " ", sub_stage$name, " setting fn_step for iter ", iter)
          opt <- set_fn_new(opt, opt$fn, iter)
        }
      }
      list(opt = opt, sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
      #message(stage$type, " ", sub_stage$name, " after step")

      alpha_old <- sub_stage$value
      # increase the step size for the next step
      alpha_new <- sub_stage$inc_fn(alpha_old)
      sub_stage$value <- sclamp(alpha_new,
                                min = sub_stage$min_value,
                                max = sub_stage$max_value)
      #message(stage$type, " ", sub_stage$name, ": step size is now = ", formatC(sub_stage$value))

      if (opt$ok && is_last_stage(opt, stage) && has_fn_new(opt, iter)) {
        # message(stage$type, " ", sub_stage$name,  ": setting next fn_curr from fn_new for ", iter)
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

# Delta-Bar-Delta ---------------------------------------------------------


delta_bar_delta <- function(kappa = 1.1, kappa_fun = `*`,
                                      phi = 0.5, epsilon = 1,
                                min_eps = 0,
                                theta = 0.1,
                                use_momentum = FALSE) {
  make_step_size(list(
    name = "delta_bar_delta",
    kappa = kappa,
    kappa_fun = kappa_fun,
    phi = phi,
    min_eps = min_eps,
    theta = theta,
    epsilon = epsilon,
    use_momentum = use_momentum,
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      sub_stage$delta_bar_old <- rep(0, length(par))
      sub_stage$gamma_old <- rep(1, length(par))
      sub_stage$value <- rep(sub_stage$init_eps, length(par))
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      delta <- opt$cache$gr_curr

      if (use_momentum && !is.null(opt$cache$update_old)) {
        # previous update includes -eps*grad_old, so reverse sign
        delta_bar_old <- -opt$cache$update_old
      }
      else {
        delta_bar_old <- sub_stage$delta_bar_old
      }
      # technically delta_bar_delta = delta_bar * delta
      # but only its sign matters, so just compare signs of delta_bar and delta
      # Force step size increase on first stage to be like the t-SNE
      # implementation
      if (all(delta_bar_old == 0)) {
        delta_bar_delta <- TRUE
      }
      else {
        delta_bar_delta <- sign(delta_bar_old) == sign(delta)
      }
      kappa <- sub_stage$kappa
      phi <- sub_stage$phi
      gamma_old <- sub_stage$gamma_old

      # signs of delta_bar and delta are the same, increase step size
      # if they're not, decrease.
      gamma <-
        (kappa_fun(gamma_old,kappa)) * abs(delta_bar_delta) +
        (gamma_old * phi) * abs(!delta_bar_delta)

      sub_stage$value <- clamp(epsilon * gamma, min_val = sub_stage$min_eps)
      if (!use_momentum || is.null(opt$cache$update_old)) {
        theta <- sub_stage$theta
        sub_stage$delta_bar_old <- ((1 - theta) * delta) + (theta * delta_bar_old)
      }
      sub_stage$gamma_old <- gamma
      list(opt = opt, sub_stage = sub_stage)
    },
    depends = c("gradient")
  ))
}
