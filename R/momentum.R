# Momentum ----------------------------------------------------------------

momentum_direction <- function(normalize = FALSE) {
  make_direction(list(
    name = "classical_momentum",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Calculating momentum direction")

      sub_stage$value <- opt$cache$update_old

      if (sub_stage$normalize) {
        sub_stage$value <- normalize(sub_stage$value)
      }
      #message("momentum pm = ", vec_formatC(sub_stage$value))
      list(sub_stage = sub_stage)

    },
    depends = c("update_old"),
    normalize = normalize
  ))
}

make_momentum_step <- function(mu_fn,
                               min_momentum = 0,
                               max_momentum = 1,
                               verbose = FALSE) {
  make_step_size(list(
    name = "momentum_step",
    init = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$t <- 1
      sub_stage$value <- sub_stage$mu_fn(0)
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      #message("calc mu step")
      sub_stage$value <-
        sclamp(sub_stage$mu_fn(sub_stage$t),
               min = sub_stage$min_value,
               max = sub_stage$max_value)
      #message("mu step " = formatC(sub_stage$value))
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {
      sub_stage$t <- sub_stage$t + 1
      #message("sub_stage$t = ", formatC(sub_stage$t))

      list(sub_stage = sub_stage)
    },
    mu_fn = mu_fn,
    min_value = min_momentum,
    max_value = max_momentum,
    t = 0
  ))
}


# Function Factories ------------------------------------------------------


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

make_constant <- function(value) {
  function(iter) {
    value
  }
}

# Momentum Correction -----------------------------------------------------

# If you want linear weighting on momentum, add an extra stage to subtract
# a fraction (mu worth) of the gradient descent

momentum_correction_direction <- function() {
  make_direction(list(
    name = "momentum_correction_direction",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Calculating momentum correction direction")

      grad_stage <- opt$stages[["gradient_descent"]]
      sub_stage$value <- -grad_stage$direction$value

      list(sub_stage = sub_stage)
    }
  ))
}

momentum_correction_step <- function() {
  make_step_size(list(
    name = "momentum_correction_step",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      #message("correcting momentum step")

      grad_stage <- opt$stages[["gradient_descent"]]
      grad_step <- grad_stage$step_size$value

      mom_stage <- opt$stages[["momentum"]]
      mom_step <- mom_stage$step_size$value

      #message("grad_step = ", formatC(grad_step), " mom_step = ", formatC(mom_step))
      sub_stage$value <- grad_step * mom_step
      list(sub_stage = sub_stage)
    }
  ))
}

# Momentum Dependencies ------------------------------------------------------------

# After step
require_update_old <- function(opt, par, fg, iter, par0, update) {
  #message("update_old: saving old update")
  opt$cache$update_old <- update
  opt
}
attr(require_update_old, 'event') <- 'after step'
attr(require_update_old, 'name') <- 'update_old'
attr(require_update_old, 'depends') <- 'update_old_init'

require_update_old_init <- function(opt, stage, sub_stage, par, fg, iter) {
  # message("Initializing update_old")
  opt$cache$update_old <- rep(0, length(par))
  list(opt = opt)
}
attr(require_update_old_init, 'event') <- 'init momentum direction'
attr(require_update_old_init, 'name') <- 'update_old_init'


