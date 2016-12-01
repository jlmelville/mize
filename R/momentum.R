# Momentum ----------------------------------------------------------------

opt_with_momentum <- function(opt, momentum_stage, nesterov = FALSE,
                              adaptive = FALSE, dec_mult = 0,
                              linear_weight = FALSE) {

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

momentum_direction <- function(normalize = FALSE) {
  make_direction(list(
    name = "classical_momentum",
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
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
                               init_momentum = 0,
                               min_momentum = 0,
                               max_momentum = 1,
                               verbose = FALSE) {
  make_step_size(list(
    name = "momentum_step",
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      sub_stage$t <- 0
      sub_stage$value <- sub_stage$init_value
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      #message("calc mu step")
      sub_stage$value <-
        sclamp(sub_stage$mu_fn(sub_stage$t),
               min = sub_stage$min_value,
               max = sub_stage$max_value)
      #message("mu step " = formatC(sub_stage$value))
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
      sub_stage$t <- sub_stage$t + 1
      #message("sub_stage$t = ", formatC(sub_stage$t))

      list(sub_stage = sub_stage)
    },
    mu_fn = mu_fn,
    init_value = init_momentum,
    min_value = min_momentum,
    max_value = max_momentum,
    t = 0
  ))
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

make_nesterov_convex_approx <- function(burn_in = 0) {
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


# Nesterov Momentum -------------------------------------------------------

# Can be considered a momentum scheme where:
# mu * [v + (v_grad - v_grad_old)]
# where v is the update vector, and v_grad and v_grad_old are the
# gradient components of the current and previous update, respectively
# i.e. replaces the gradient component of the previous velocity with the
# gradiemt velocity of the current iteration
nesterov_momentum_direction <- function() {
  make_direction(list(
    name = "nesterov",
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      sub_stage$value <- rep(0, length(par))
      #sub_stage$grad_update_old <- rep(0, length(par))
      sub_stage$update <- rep(0, length(par))
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      #message("Calculating nesterov momentum direction")

      # update_old <- opt$cache$update_old
      # grad_update_old <- sub_stage$grad_update_old
      # grad_update <- opt$stages[["gradient_descent"]]$result
      # sub_stage$value <- update_old + (grad_update - grad_update_old)



      grad_update <- opt$stages[["gradient_descent"]]$result
      sub_stage$value <- grad_update + sub_stage$update


      list(sub_stage = sub_stage)

    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
 #     sub_stage$update <- stage$result
 #     sub_stage$grad_update_old <- opt$stages[["gradient_descent"]]$result
      sub_stage$update <- update - opt$stages[["gradient_descent"]]$result
      list(sub_stage = sub_stage)
    }
 # ,
#    depends = c("update_old")
  ))
}

# q is inversely proportional to how strongly convex the function is
# 0 gives the highest momentum, 1 gives zero momentum
nesterov_convex_step <- function(start_at = 0, q = 0) {
  if (q == 0) {
    # Use the expression for momentum from the Sutskever paper appendix
    nesterov_strong_convex_step(start_at = start_at)
  }
  else {
    # Use the expression for momentum from Candeis paper which includes q term
    nesterov_convex_step_q(q = q, start_at = start_at)
  }
}


nesterov_strong_convex_step <- function(start_at) {
  make_step_size(list(
    start_at = start_at,
    name = "nesterov_convex",
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      #message("Nesterov convex init")
      sub_stage$a_old <- 1
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      if (iter < start_at) {
        sub_stage$value <- 0
        sub_stage$a <- 1
      }
      else {
        a_old <- sub_stage$a_old
        a <- (1 + sqrt(4 * a_old * a_old + 1)) / 2
        sub_stage$value <- (a_old - 1) / a
        sub_stage$a <- a
        #message("Nesterov momentum = ", formatC(sub_stage$value))
      }
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
      #message("nesterov_convex: after step")
      if (!opt$ok) {
        sub_stage$a_old <- 1
      }
      else {
        sub_stage$a_old <- sub_stage$a
      }
      list(sub_stage = sub_stage)
    }
  ))
}

nesterov_convex_step_q <- function(q, start_at = 0) {
  make_step_size(list(
    start_at = start_at,
    name = "nesterov_convex",
    init = function(opt, stage, sub_stage, par, fn, gr, iter) {
      #message("Nesterov convex init")
      sub_stage$theta_old <- 1
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
      if (iter < start_at) {
        sub_stage$value <- 0
        sub_stage$theta_old <- 1
      }
      else {
        theta_old <- sub_stage$theta_old
        thetas <- solve_theta(theta_old, q)
        theta <- max(thetas)
        # Step 4 of algorithm 1 in https://arxiv.org/abs/1204.3982
        # Calculates beta, effectively the momentum
        # q (taking a value between 0 and 1) is related to the strong convexity
        # parameter (which has the symbol mu in the paper, it's not the momentum!).
        # A q of 1 causes momentum to be zero. A q of 0 gives the same results as
        # the Sutskever momentum (not the approximation he gives, but the actual
        # expression given in the appendix/thesis).
        sub_stage$value <- theta_old * (1 - theta_old) / (theta_old * theta_old + theta)
        sub_stage$theta_old <- theta
        #message("Nesterov momentum = ", formatC(sub_stage$value))
      }
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fn, gr, iter, par0,
                          update) {
      #message("nesterov_convex: after step")
      if (!opt$ok) {
        sub_stage$theta_old <- 1
      }
      list(sub_stage = sub_stage)
    }
  ))
}


# solves quadratic equation. Returns either two roots (even if coincident)
# or NULL if there's no solution
solve_quad <- function(a, b, c) {
  disc <- b * b - 4 * a * c
  res <- c()
  if (disc > 0) {
    root_pos = (-b + sqrt(disc)) / (2 * a)
    root_neg = (-b - sqrt(disc)) / (2 * a)
    res <- c(root_pos, root_neg)
  }
}

# Step 3 of algorithm 1 in https://arxiv.org/abs/1204.3982
solve_theta <- function(theta_old, q = 0) {
  theta2 <- theta_old * theta_old
  solve_quad(1, theta2 - q, -theta2)
}

# Momentum Correction -----------------------------------------------------


momentum_correction_direction <- function() {
  make_direction(list(
    name = "momentum_correction_direction",
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
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
    calculate = function(opt, stage, sub_stage, par, fn, gr, iter) {
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



