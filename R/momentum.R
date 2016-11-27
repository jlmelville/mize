# Momentum ----------------------------------------------------------------

opt_with_momentum <- function(opt, momentum_stage, nesterov = FALSE,
                              adaptive = FALSE, dec_mult = 0,
                              linear_weight = FALSE) {
  if (adaptive) {
    #    momentum_stage <- adaptive_restart2(momentum_stage, dec_mult = dec_mult)
  }
  # linear weighting and nesterov both need the momentum stage to be calculated
  # first (but for different reasons) so we'll deal with them in the same block
  if (nesterov || linear_weight) {
    opt <- prepend_stage(opt, momentum_stage)

    if (linear_weight) {
      # because momentum stage happens first, we can just look at the step
      # size value for that stage when calculating the gradient descent
      # step size
      if (is.null(opt$stages[["gradient_descent"]])) {
        stop("Can't linear weight a momentum optimizer without a gradient ",
             "descent stage")
      }

      gd_step_calc <- opt$stages[["gradient_descent"]]$step_size$calculate
      new_calc <- function(opt, stage, par, fn, gr, iter) {
        res <- gd_step_calc(opt, stage, par, fn, gr, iter)
        stage <- res$stage
        if (!is.null(res$opt)) {
          opt <- res$opt
        }

        mu <- opt$stages[["momentum"]]$step_size$value
        wt <- 1 - mu
        stage$step_size$value <- wt * stage$step_size$value

        list(opt = opt, stage = stage)
      }
      opt$stages[["gradient_descent"]]$step_size$calculate <- new_calc
    }

    if (nesterov) {
      opt$eager_update <- TRUE
    }
  }
  else {
    opt <- append_stage(opt, momentum_stage)
  }
  opt
}

# Wants: update_old
momentum_direction <- function(normalize = FALSE) {
  make_direction(list(
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      message("Calculating momentum direction")

      sub_stage$value <- opt$cache$update_old

      if (sub_stage$normalize) {
        sub_stage$value <- normalize(sub_stage$value)
      }
      message("momentum pm = ", vec_formatC(sub_stage$value))
      list(sub_stage = sub_stage)

    },
    depends = c("update_old"),
    normalize = normalize
  ))
}

make_momentum_step <- function(mu_fn,
                               init_momentum = 0,
                               min_momentum = 0,
                               max_momentum = 1,
                               verbose = FALSE) {
  list(
    init = function(opt, stage, par, fn, gr, iter) {
      stage$step_size$t <- 0
      stage$step_size$value <- stage$step_size$init_value
      list(opt = opt, stage = stage)
    },
    calculate = function(opt, stage, par, fn, gr, iter) {
      stage$step_size$value <-
        sclamp(stage$step_size$mu_fn(stage$step_size$t),
               min = stage$step_size$min_value,
               max = stage$step_size$max_value)
      list(stage = stage)
    },
    after_step = function(opt, stage, par, fn, gr, iter, par0, update) {
      stage$step_size$t <- stage$step_size$t + 1
      message("stage$step_size$t = ", formatC(stage$step_size$t))

      list(stage = stage)
    },
    mu_fn = mu_fn,
    init_value = init_momentum,
    min_value = min_momentum,
    max_value = max_momentum,
    t = 0
  )
}

make_switch <- function(init_value = 0.5, final_value = 0.8,
                        switch_iter = 250) {
  function(iter) {
    if (iter >= switch_iter) {
      return(final_value)
    }
    else {
      return(init_value)
    }
  }
}

make_ramp <- function(max_iter,
                      init_value = 0,
                      final_value = 0.9, ...) {
  function(iter) {
    ds <- (final_value - init_value) / max_iter
    (ds * iter) + init_value
  }
}

make_nesterov_nsc <- function(burn_in = 0) {
  function(iter) {
    if (iter < burn_in) {
      return(0)
    }
    1 - (3 / ((iter - burn_in) + 5))
  }
}

make_constant <- function(value) {
  function(iter) {
    value
  }
}


# Adaptive Restart --------------------------------------------------------

append_depends <- function(opt, stage = NULL, sub_stage = NULL, new_depends) {
  if (!is.null(sub_stage)) {
    if (is.null(stage)) {
      stop("Must provide stage for sub_stage '", sub_stage, "'")
    }
    depends <- opt$stages[[stage]][[sub_stage]]$depends
  }
  else if (!is.null(stage)) {
    depends <- opt$stages[[stage]]$depends
  }
  else {
    depends <- opt$depends
  }
  if (is.null(depends)) {
    depends <- c()
  }

  depends <- c(depends, new_depends)

  if (!is.null(sub_stage)) {
    opt$stages[[stage]][[sub_stage]]$depends <- depends
  }
  else if (!is.null(stage)) {
    opt$stages[[stage]]$depends <- depends
  }
  else {
    opt$depends <- depends
  }

  opt
}

adaptive_restart <- function(opt, validation_type) {
  append_depends(opt, "momentum", "direction",
                  c('adaptive_restart', paste0('validate_', validation_type)))
}

adaptive_restart_fn <- function(opt) {
  adaptive_restart(opt, "fn")
}

adaptive_restart_gr <- function(opt) {
  adaptive_restart(opt, "gr")
}



