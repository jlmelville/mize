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
    init = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$value <- rep(0, length(par))
      #sub_stage$grad_update_old <- rep(0, length(par))
      sub_stage$update <- rep(0, length(par))
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Calculating nesterov momentum direction")

      # update_old <- opt$cache$update_old
      # grad_update_old <- sub_stage$grad_update_old
      # grad_update <- opt$stages[["gradient_descent"]]$result
      # sub_stage$value <- update_old + (grad_update - grad_update_old)

      grad_update <- opt$stages[["gradient_descent"]]$result
      sub_stage$value <- grad_update + sub_stage$update


      list(sub_stage = sub_stage)

    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
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

# Approximate Nesterov Convex Momentum Function Factory
# Can be applied to make_momentum_step
make_nesterov_convex_approx <- function(burn_in = 0) {
  function(iter) {
    if (iter < burn_in) {
      return(0)
    }
    1 - (3 / ((iter - burn_in) + 5))
  }
}

# Sutskever's approximation to Nesterov Momentum Scheme
nesterov_convex_approx_step <- function(burn_in = 0) {
  make_momentum_step(mu_fn =
                       make_nesterov_convex_approx(burn_in = burn_in),
                     min_momentum = 0,
                     max_momentum = 1)
}

# q is inversely proportional to how strongly convex the function is
# 0 gives the highest momentum, 1 gives zero momentum
nesterov_convex_step <- function(burn_in = 0, q = 0) {
  if (q == 0) {
    # Use the expression for momentum from the Sutskever paper appendix
    nesterov_strong_convex_step(burn_in = burn_in)
  }
  else {
    # Use the expression for momentum from Candeis paper which includes q term
    nesterov_convex_step_q(q = q, burn_in = burn_in)
  }
}


nesterov_strong_convex_step <- function(burn_in) {
  make_step_size(list(
    burn_in = burn_in,
    name = "nesterov_convex",
    init = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Nesterov convex init")
      sub_stage$a_old <- 1
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      if (iter < burn_in) {
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
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
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

nesterov_convex_step_q <- function(q, burn_in = 0) {
  make_step_size(list(
    burn_in = burn_in,
    name = "nesterov_convex",
    init = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Nesterov convex init")
      sub_stage$theta_old <- 1
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      if (iter < burn_in) {
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
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
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
