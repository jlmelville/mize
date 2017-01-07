# Momentum ----------------------------------------------------------------

# Create a direction sub stage for momentum
momentum_direction <- function(normalize = FALSE) {
  make_direction(list(
    name = "classical_momentum",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {

      sub_stage$value <- opt$cache$update_old

      if (sub_stage$normalize) {
        sub_stage$value <- normalize(sub_stage$value)
      }
      list(sub_stage = sub_stage)
    },
    depends = c("update_old"),
    normalize = normalize
  ))
}

# Creates a step size sub stage for momentum
# mu_fn a function that takes an iteration number and returns the momentum.
#  Adaptive restart can restart the momentum in which case the function will
#  be passed an "effective" iteration number which may be smaller than the
#  actual iteration value.
# use_init_mom If TRUE, then always use the momentum coefficient specified by
#  mu_fn even when the effective iteration is 1 (first iteration or restart).
#  In some cases using non-standard momentum (e.g. Nesterov or linear-weighted),
#  this could result in the resulting step length being shorter or longer
#  than would be otherwise expected. If FALSE, then the momentum coefficient
#  is always zero.
make_momentum_step <- function(mu_fn,
                               min_momentum = 0,
                               max_momentum = 1,
                               use_init_mom = FALSE,
                               verbose = FALSE) {
  make_step_size(list(
    name = "momentum_step",
    init = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$value <- 0
      sub_stage$t <- 1
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      if (!use_init_mom && sub_stage$t <= 1) {
        sub_stage$value <- 0
      }
      else {
        sub_stage$value <-
          sclamp(sub_stage$mu_fn(sub_stage$t, opt$convergence$max_iter),
                 min = sub_stage$min_value,
                 max = sub_stage$max_value)
      }
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {
      sub_stage$t <- sub_stage$t + 1

      list(sub_stage = sub_stage)
    },
    mu_fn = mu_fn,
    min_value = min_momentum,
    max_value = max_momentum,
    t = 0
  ))
}


# Function Factories ------------------------------------------------------

# A function that switches from one momentum value to another at the
# specified iteration.
make_switch <- function(init_value = 0.5, final_value = 0.8,
                        switch_iter = 250) {
  function(iter, max_iter) {
    if (iter >= switch_iter) {
      return(final_value)
    }
    else {
      return(init_value)
    }
  }
}

# A function that increases from init_value to final_value over
# max_iter iterations. Iter 0 will always return a value of zero, iter 1
# begins with init_value.
#
# wait - if set to a non-zero value, recalculates the values so that
# the init_value is used for 'wait' extra iterations, but with final_value
# still reached after max_iter iterations. Set to 1 for momentum calculations
# where in most cases the momentum on the first iteration would be either
# ignored or the value overridden and set to zero anyway. Stops a larger than
# expected jump on iteration 2.
make_ramp <- function(init_value = 0,
                      final_value = 0.9,
                      wait = 0) {

  function(iter, max_iter) {
    # actual number of iterations
    iters <- max_iter - 1 - wait
    # denominator of linear scaling
    n <- max(iters, 1)
    m <- (final_value - init_value) / n
    t <- iter - 1 - wait
    if (t < 0) {
      return(init_value)
    }

    (m * t) + init_value
  }
}

# A function that returns a constant momentum value
make_constant <- function(value) {
  function(iter, max_iter) {
    value
  }
}

# Momentum Correction -----------------------------------------------------

# Normally, momentum schemes are given as eps*grad + mu*old_update, but
# some momentum schemes define the update as: (1-mu)*eps*grad + mu*old_update
# which can easily be expanded as: eps*grad + mu*old_update - mu*eps*grad
# i.e. add an extra stage to substract a fraction (mu worth) of the gradient
# descent

# The momentum correction direction: the opposite direction the gradient
# descent.
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

# The momentum correction step size: mu times the gradient descent step size.
momentum_correction_step <- function() {
  make_step_size(list(
    name = "momentum_correction_step",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {

      grad_stage <- opt$stages[["gradient_descent"]]
      grad_step <- grad_stage$step_size$value

      mom_stage <- opt$stages[["momentum"]]
      mom_step <- mom_stage$step_size$value

      sub_stage$value <- grad_step * mom_step
      list(sub_stage = sub_stage)
    }
  ))
}

# Momentum Dependencies ------------------------------------------------------------

# Save this update for use in the next step
require_update_old <- function(opt, par, fg, iter, par0, update) {
  opt$cache$update_old <- update
  opt
}
attr(require_update_old, 'event') <- 'after step'
attr(require_update_old, 'name') <- 'update_old'
attr(require_update_old, 'depends') <- 'update_old_init'

# Initialize the old update vector
require_update_old_init <- function(opt, stage, sub_stage, par, fg, iter) {
  opt$cache$update_old <- rep(0, length(par))
  list(opt = opt)
}
attr(require_update_old_init, 'event') <- 'init momentum direction'
attr(require_update_old_init, 'name') <- 'update_old_init'


