# Nesterov Accelerated Gradient ------------------------------------------------

# This is the actual Nesterov Accelerated Gradient scheme, rather than
# the version discussed by Sutskever and popular in the deep learning community,
# although that is also available in Mize.
#
# NAG can be considered to be an optimization consisting of:
# 1. A steepest descent step;
# 2. A pseudo-momentum step, using a specific schedule, depending on how
#    convex the function being optimized is.
# The pseudo-momentum step is:
# mu * [v + (v_grad - v_grad_old)]
# where v is the update vector, and v_grad and v_grad_old are the
# gradient components of the current and previous update, respectively.
# Overall, it replaces the gradient component of the previous velocity with the
# gradient velocity of the current iteration.
nesterov_momentum_direction <- function() {
  make_direction(list(
    name = "nesterov",
    init = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$value <- rep(0, length(par))
      sub_stage$update <- rep(0, length(par))
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      grad_update <- opt$stages[["gradient_descent"]]$result
      sub_stage$value <- grad_update + sub_stage$update
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {
      sub_stage$update <- update - opt$stages[["gradient_descent"]]$result
      list(sub_stage = sub_stage)
    }
  ))
}

# Encapsulates all ways to create a Nesterov schedule
# burn_in Lags the calculation by this number of iterations. By setting this
#   to 2, you get the same "pattern" of results as if you were using the
#   Sutskever Nesterov Momentum approach (i.e. applying a classical momentum
#   step before a steepest descent step)
# q is inversely proportional to how strongly convex the function is
#   0 gives the highest momentum, 1 gives zero momentum. Often, q is assumed to
#   be zero. Ignored if use_approx is TRUE.
# use_approx Use the approximation to the momentum schedule given by
#   Sutskever and co-workers.
# use_init_mu If TRUE, then the momentum calculated on the first iteration uses
#   the calculated non-zero value, otherwise use zero. Because velocity is
#   normally zero initially, this rarely has an effect, unless linear weighting
#   of the momentum is being used. Ignored if use_approx is FALSE.
nesterov_step <- function(burn_in = 0, q = 0, use_approx = FALSE,
                          use_init_mu = FALSE) {
  if (!is_in_range(q, 0, 1)) {
    stop("q must be between 0 and 1")
  }
  if (burn_in < 0) {
    stop("burn_in must be non-negative")
  }

  if (use_approx) {
    nesterov_convex_approx_step(burn_in = burn_in,
                                use_init_mu = use_init_mu)
  }
  else {
    nesterov_convex_step(q = q, burn_in = burn_in)
  }
}

# Approximate Nesterov Convex Momentum Function Factory
#
# Instead of using the exact momentum schedule specified for NAG, use the
# approximation given by Sutskever. NB: this produces much larger momentums
# than the exact result for the first few iterations.
#
# burn_in Lags the calculation by this number of iterations. By setting this
#  to 2, you get the same "pattern" of results as if you were using the
#  Sutskever Nesterov Momentum approach (i.e. applying a classical momentum
#  step before a steepest descent step).
# use_init_mu If TRUE, then return a non-zero momentum on the first iteration.
#  Otherwise use zero. Although a velocity of zero normally enforces
#  steepest descent on the first iteration, for some methods (e.g.
#  NAG or linearly weighted classical momentum), this can have an effect.
#  Set this to TRUE to always get steepest decent.
make_nesterov_convex_approx <- function(burn_in = 0, use_init_mu = FALSE) {
  function(iter, max_iter) {
    # if we haven't waited long enough or we always use zero on the first
    # iteration, return 0
    if (iter < burn_in || (iter == burn_in && !use_init_mu)) {
      return(0)
    }

    1 - (3 / ((iter - burn_in) + 5))
  }
}

# Create a momentum step size sub stage using
# Sutskever's approximation to the NAG pseudo-momentum schedule.
# burn_in Lags the calculation by this number of iterations. By setting this
#  to 2, you get the same "pattern" of results as if you were using the
#  Sutskever Nesterov Momentum approach (i.e. applying a classical momentum
#  step before a steepest descent step).
# use_init_mu if TRUE, then on the first iteration, the momentum
#  uses the equation, which produces a momentum of 0.4. Otherwise, use a
#  momentum coefficient of zero. From reading various papers, a coefficient
#  of zero is probably the intended behavior.
#  Normally, the velocity vector is also zero on the first iteration, so this
#  makes no difference, but if you have linearly weighted the momentum,
#  you will get only 60% of the gradient step you might have been expecting
#  on the first step, and you will get a 60% longer step size if using NAG.
nesterov_convex_approx_step <- function(burn_in = 0, use_init_mu = FALSE) {
  make_momentum_step(mu_fn =
                       make_nesterov_convex_approx(burn_in = burn_in,
                                                   use_init_mu = use_init_mu),
                     min_momentum = 0,
                     max_momentum = 1,
                     use_init_mom = use_init_mu)
}

# The NAG pseudo-momentum schedule.
# q is inversely proportional to how strongly convex the function is
# 0 gives the highest momentum, 1 gives zero momentum. Often, q is assumed to be
# zero.
# burn_in Lags the calculation by this number of iterations. By setting this
#  to 2, you get the same "pattern" of results as if you were using the
#  Sutskever Nesterov Momentum approach (i.e. applying a classical momentum
#  step before a steepest descent step)
nesterov_convex_step <- function(burn_in = 0, q = 0) {
  if (q == 0) {
    # Use the expression for momentum from the Sutskever paper appendix
    nesterov_strong_convex_step(burn_in = burn_in)
  }
  else {
    # Use the expression for momentum from Candes paper which includes q term
    nesterov_convex_step_q(q = q, burn_in = burn_in)
  }
}

# The NAG pseudo-momentum schedule for strongly convex functions.
# This expression is missing the parameter "q" that measures how strongly convex
# the function is. It is implicitly zero, which gives the largest momentum
# values. This uses an expression in the appendix of the Sutskever paper, which
# is a bit simpler to calculate than the version given by O'Donoghue and Candes.
nesterov_strong_convex_step <- function(burn_in) {
  make_step_size(list(
    burn_in = burn_in,
    name = "nesterov_convex",
    init = function(opt, stage, sub_stage, par, fg, iter) {
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
      }
      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                          update) {
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

# The NAG pseudo-momentum schedule. This expression includes the parameter "q"
# that measures how strongly convex This uses the algorithm given by
# O'Donoghue and Candes, which is a bit more complex than the one in the
# appendix of the Sutskever paper (but that one assumes q = 0).
# See https://arxiv.org/abs/1204.3982 for more.
nesterov_convex_step_q <- function(q, burn_in = 0) {
  make_step_size(list(
    burn_in = burn_in,
    name = "nesterov_convex",
    init = function(opt, stage, sub_stage, par, fg, iter) {
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
  res
}

# Step 3 of algorithm 1 in https://arxiv.org/abs/1204.3982
solve_theta <- function(theta_old, q = 0) {
  theta2 <- theta_old * theta_old
  solve_quad(1, theta2 - q, -theta2)
}
