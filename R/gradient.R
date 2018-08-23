# Various Gradient-based optimization routines: steepest descent, conjugate
# gradient, BFGS etc.

# Gradient Direction -----------------------------------------------------------

# Creates a direction sub stage
make_direction <- function(sub_stage) {
  make_sub_stage(sub_stage, 'direction')
}

# Steepest Descent
#
# normalize - If TRUE, then the returned direction vector is normalized to unit
# length. This can be useful for some adaptive line search methods, so that the
# total step length is purely determined by the line search value, rather than
# the product of the line search value and the magnitude of the direction.
sd_direction <- function(normalize = FALSE) {

  make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      opt$cache$gr_curr <- NULL
      list(opt = opt)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$value <- -opt$cache$gr_curr

      if (sub_stage$normalize) {
        sub_stage$value <- normalize(sub_stage$value)
      }

      list(sub_stage = sub_stage)
    },
    normalize = normalize
  ))
}


# Conjugate Gradient ------------------------------------------------------

# Conjugate gradient
#
# ortho_check - If TRUE, check successive direction are sufficiently orthogonal.
#   If the orthogonality check is failed, then the next step is steepest descent.
# nu - the orthogonality threshold. Used only if ortho_check is TRUE. Compared
#   with g_old . g_new / g_new . g_new
# cg_update - Function to generate the next direction using a method of e.g.
#   Fletcher-Reeves or Polak-Ribiere. Pass one of the cg_update functions
#   below, e.g. pr_plus_update
cg_direction <- function(ortho_check = FALSE, nu = 0.1,
                         cg_update = pr_plus_update,
                         eps = .Machine$double.eps) {
  make_direction(list(
    ortho_check = ortho_check,
    nu = nu,
    cg_update = cg_update,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      sub_stage$value <- rep(0, length(par))
      sub_stage$pm_old <- rep(0, length(par))

      opt$cache$gr_curr <- NULL
      opt$cache$gr_old <- NULL

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      gm_old <- opt$cache$gr_old
      pm_old <- sub_stage$pm_old

      # direction is initially steepest descent
      pm <- -gm

      if (!is.null(gm_old)
          && (!sub_stage$ortho_check
              || !sub_stage$cg_restart(gm, gm_old, sub_stage$nu))) {
        beta <- sub_stage$cg_update(gm, gm_old, pm_old, sub_stage$eps)
        pm <- pm + (beta * pm_old)
        descent <- dot(gm, pm)
        if (descent >= 0) {
          #message("Next CG direction is not a descent direction, resetting to SD")
          pm <- -gm
        }
      }

      sub_stage$value <- pm

      list(sub_stage = sub_stage)
    }
   , after_step = function(opt, stage, sub_stage, par, fg, iter, par0,
                         update) {
     sub_stage$pm_old <- sub_stage$value
     list(sub_stage = sub_stage)
   }
    , depends = c("gradient_old")
  ))
}


# CG update formulae, grouped according to their numerators similar to the
# discussion in Hager and Zhang's survey paper
# The FR, CD and DY updates are all susceptible to "jamming": they can end up
# with very small step sizes and make little progress.
# The Fletcher-Reeves update.
fr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  dot(gm) / (dot(gm_old) + eps)
}

# Conjugate Descent update due to Fletcher
cd_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  dot(gm) / (dot(pm_old, (gm - gm_old)) + eps)
}

# The Dai-Yuan update.
dy_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  -cd_update(gm, gm_old, pm_old, eps)
}

# HS, PR and LS share a numerator. According to Hager and Zhang, they
# perform better in practice than the FR, CD and DY updates, despite less
# being known about their provable global convergence properties.

# The Hestenes-Stiefel update.
hs_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  -(dot(gm, gm_old) - dot(gm, gm_old)) / (dot(pm_old, (gm - gm_old)) + eps)
}

# An "HS+" modification of Hestenes-Stiefel, in analogy to the "PR+" variant of
# Polak-Ribiere suggested by Powell. As far as I can tell, Hager and Zhang
# suggested this modification.
# Hager, W. W., & Zhang, H. (2006).
# A survey of nonlinear conjugate gradient methods.
# \emph{Pacific journal of Optimization}, \emph{2}(1), 35-58.
hs_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  beta <- hs_update(gm, gm_old, pm_old, eps)
  max(0, beta)
}

# The Polak-Ribiere method for updating the CG direction. Also known as
# Polak-Ribiere-Polyak (PRP)
pr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  dot(gm, gm - gm_old) / (dot(gm_old) + eps)
}

# The "PR+" update due to Powell. Polak-Ribiere update, but if negative,
# restarts the CG from steepest descent. Prevents a possible lack of
# convergence when using a Wolfe line search.
pr_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  beta <- pr_update(gm, gm_old, pm_old, eps)
  max(0, beta)
}

# Liu-Storey update
ls_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  -hs_update(gm, gm_old, pm_old, eps)
}

# Hager-Zhang update as used in CG_DESCENT
hz_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  ym <- gm - gm_old
  py <- dot(pm_old, ym)
  dot(ym - 2 * pm_old * (dot(ym) / (py + eps)), (gm / (py + eps)))
}

# "Restricted" Hager-Zhang update as used in CG_DESCENT to ensure
# convergence. Analogous to the PR+ and HS+ updates, but dynamically adjusts
# the lower bound as convergence occurs. Choice of eta is from the CG_DESCENT
# paper
hz_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  beta <- hz_update(gm, gm_old, pm_old, eps)
  eta <- 0.01
  eta_k <- -1 / (dot(pm_old) * min(eta, dot(gm_old)))
  max(eta_k, beta)
}


# Restart criteria due to Powell
# Checks that successive gradient vectors are sufficiently orthogonal
# g_new . g_old / g_new . g_new  must be greater than or equal to nu.
cg_restart <- function(g_new, g_old, nu = 0.1) {
  # could only happen on first iteration
  if (is.null(g_old)) {
    return(TRUE)
  }
  ortho_test <- abs(dot(g_new, g_old)) / dot(g_new)
  should_restart <- ortho_test >= nu
  should_restart
}


# BFGS --------------------------------------------------------------------

# The Broyden Fletcher Goldfarb Shanno method
# scale_inverse - if TRUE, scale the inverse Hessian approximation on the first
#   step.
bfgs_direction <- function(eps = .Machine$double.eps,
                           scale_inverse = FALSE) {
  make_direction(list(
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      n <- length(par)
      sub_stage$value <- rep(0, n)

      if (!is.null(fg$hi)) {
        hm <- fg$hi(par)
        if (methods::is(hm, "numeric")) {
          # Allow for just a vector to be passed, representing a diagonal
          # approximation to the Hessian inverse. But we'll store it as the
          # full matrix.
          hm <- diag(hm)
        }
        sub_stage$hm <- hm
      }
      else {
        sub_stage$hm <- diag(1, n)
      }

      opt$cache$gr_curr <- NULL
      opt$cache$gr_old <- NULL

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      gm_old <- opt$cache$gr_old
      if (!is.null(gm_old)) {
        sm <- opt$cache$update_old
        hm <- sub_stage$hm

        ym <- gm - gm_old

        # Nocedal suggests this heuristic for scaling the first
        # approximation in Chapter 6 Section "Implementation"
        # Also used in the definition of L-BFGS
        if (iter == 2 && scale_inverse) {
          gamma <- dot(sm, ym) / dot(ym)
          hm <- gamma * hm
        }

        sub_stage$hm <- bfgs_update(hm, sm, ym, eps = sub_stage$eps)
      }
      pm <- as.vector(-sub_stage$hm %*% gm)

      descent <- dot(gm, pm)
      if (descent >= 0) {
        pm <- -gm
      }

      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
    , depends = c("gradient_old", "update_old")
  ))
}

# H - inverse Hessian
# s - difference in iterates x
# y - difference in gradient of x
bfgs_update <- function(hm, sm, ym, eps) {
  rho <- 1 / (dot(ym, sm) + eps)
  im <- diag(1, nrow(hm))

  rss <- rho * outer(sm, sm)
  irsy <- im - rho * outer(sm, ym)
  irys <- im - rho * outer(ym, sm)

  (irsy %*% (hm %*% irys)) + rss
}

# SR1 ---------------------------------------------------------------------

# The Symmetric-Rank-1 Update. See section 6.2 of Nocedal and Wright.
# A Broyden-class Quasi Newton method, like BFGS, but using a rank-1 matrix
# to update the inverse Hessian approximation (i.e. the outer product of a
# vector).
# The good news: SR1 updates are sometimes better approximations to the Hessian
#  than BFGS (according to Nocedal and Wright).
# The bad news: SR1 updates can contain a division by zero.
# The less-bad news: we can detect this and use the last Hessian approximation
#  in these cases.
# The bad news 2: SR1 updates are not guaranteed to result in a positive
#  definite Hessian, i.e. you won't always get a descent direction. As a result,
#  Nocedal & Wright suggest using it with a trust-region approach rather than
#  line search. Here we simply replace the SR1 update with the BFGS version to
#  ensure a descent direction.
sr1_direction <- function(eps = .Machine$double.eps,
                          scale_inverse = FALSE, skip_bad_update = TRUE,
                          tol = 1e-8, try_bfgs = TRUE) {
  make_direction(list(
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      n <- length(par)
      sub_stage$value <- rep(0, n)

      if (!is.null(fg$hi)) {
        hm <- fg$hi(par)
        if (methods::is(hm, "numeric")) {
          # Allow for just a vector to be passed, representing a diagonal
          # approximation to the Hessian inverse. But we'll store it as the
          # full matrix.
          hm <- diag(hm)
        }
        sub_stage$hm <- hm
      }
      else {
        sub_stage$hm <- diag(1, n)
      }

      opt$cache$gr_curr <- NULL
      opt$cache$gr_old <- NULL

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      gm_old <- opt$cache$gr_old
      if (!is.null(gm_old)) {
        sm <- opt$cache$update_old
        hm <- sub_stage$hm

        ym <- gm - gm_old

        if (iter == 2 && scale_inverse) {
          # Nocedal suggests this heuristic for scaling the first
          # approximation in Chapter 6 Section "Implementation" for BFGS
          # Also used in the definition of L-BFGS
          gamma <- dot(sm, ym) / dot(ym)
          hm <- gamma * hm
        }

        shy <- as.vector(sm - hm %*% ym)
        up_den <- dot(shy, ym)

        # Skip rule based on 6.26 in Nocedal & Wright
        # Their version uses strict < rather than <=, but we want to catch
        # the case where one of shy or ym is exactly zero
        if (skip_bad_update && abs(up_den) <= tol * norm2(shy) * norm2(ym)) {
          # message("SR1: skipping bad update")
          sub_stage$hm <- hm
        }
        else {
          sub_stage$hm <- hm + outer(shy, shy) / up_den
        }
      }
      pm <- as.vector(-sub_stage$hm %*% gm)

      descent <- dot(gm, pm)
      if (descent >= 0) {
        if (try_bfgs) {
          # message("SR1 iter ", iter, " not a descent direction, trying BFGS")
          hm_bfgs <- bfgs_update(hm, sm, ym, sub_stage$eps)
          pm <- as.vector(-hm_bfgs %*% gm)
          descent <- dot(gm, pm)
          if (descent >= 0) {
            # This should never happen: BFGS H should always be pd.
            pm <- -gm
          }
          else {
            # Only store the BFGS H for next iteration if it's a descent
            sub_stage$hm <- hm_bfgs
          }
        }
        else {
          pm <- -gm
        }
      }

      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
    , depends = c("gradient_old", "update_old")
  ))
}

# L-BFGS ------------------------------------------------------------------

# Initial guess for solving Bp = q for p by using an approximation
# for H, the inverse of B, and calculating Hq directly.
# You can provide the inverse hessian function fg$hi and par, which can return
# either a matrix (although if it's dense that's not a great idea), or a vector
# in which case it's assumed to represent a diagonal matrix.
# Or you can set scale_inverse to TRUE and the usual L-BFGS guess for H as
# given by Nocedal and Wright is used. rho and eps may not be NULL in that case.
# Or you can set scale_inverse = FALSE and you will get q back directly, i.e.
# Assume use H = I, which is also used when scale_inverse = TRUE on the first
# iteration anyway
lbfgs_guess <- function(qm, scale_inverse = TRUE,
                        rho = NULL, ym = NULL, eps = NULL,
                        fg = NULL, par = NULL) {
  # User-defined hessian inverse approximations
  if (!is.null(fg) && !is.null(fg$hi)) {
    hm <- fg$hi(par)
    if (methods::is(hm, "matrix")) {
      # Full matrix: not necessarily a great idea for memory usage
      pm <- hm %*% qm
    }
    else {
      # It's a vector representing a diagonal matrix
      pm <- hm * qm
    }
  }
  else {
    # The usual L-BFGS guess
    if (!is.null(rho) && !is.null(ym) && scale_inverse) {
      # Eqn 7.20 in Nocedal & Wright
      gamma <- 1 / (rho * (dot(ym) + eps))
      hm <- rep(gamma, length(par))
      pm <- hm * qm
    }
    else {
      # Effectively H = I
      # This will happen on the first iteration
      pm <- qm
    }
  }
  pm
}

# Solve Bp = q for p by the two-loop recursion. If fn and par are non NULL
# and fg has the Hessian inverse function hi defined, it will be used to
# initialize the inital guess for p. Otherwise, if scale_inverse = TRUE, the
# usual L-BFGS guess is used. Otherwise, q is used.
lbfgs_solve <- function(qm, lbfgs_state, scale_inverse, eps,
                        fg = NULL, par = NULL) {
  sms <- lbfgs_state$sms
  yms <- lbfgs_state$yms
  rhos <- lbfgs_state$rhos

  alphas <- rep(0, length(rhos))
  # loop backwards latest values first
  for (i in length(rhos):1) {
    alphas[i] <- rhos[[i]] * dot(sms[[i]], qm)
    qm <- qm - alphas[[i]] * yms[[i]]
  }

  # Choose estimate of inverse Hessian approximation
  pm <- lbfgs_guess(qm, scale_inverse, rhos[[length(rhos)]], yms[[length(yms)]],
                    eps, fg, par)

  # loop forwards
  for (i in 1:length(rhos)) {
    beta <- rhos[[i]] * dot(yms[[i]], pm)
    pm <- pm + sms[[i]] * (alphas[[i]] - beta)
  }

  pm
}

# Update lbfgs_state with new ym, sm, rho, removing old values if we've
# reached the memory limit
# One real calculation (rho), otherwise just lots of boring housekeeping
lbfgs_memory_update <- function(lbfgs_state, ym, sm, eps) {
  rho <- 1 / (dot(ym, sm) + eps)

  rhos <- lbfgs_state$rhos
  sms <- lbfgs_state$sms
  yms <- lbfgs_state$yms

  # discard oldest values if we've reached memory limit
  if (length(sms) == lbfgs_state$memory) {
    sms <- sms[2:length(sms)]
    yms <- yms[2:length(yms)]
    rhos <- rhos[2:length(rhos)]
  }

  # append latest values to memory
  sms <- c(sms, list(sm))
  yms <- c(yms, list(ym))
  rhos <- c(rhos, list(rho))

  lbfgs_state$sms <- sms
  lbfgs_state$yms <- yms
  lbfgs_state$rhos <- rhos

  lbfgs_state
}

# The Limited Memory BFGS method
#
# memory - The number of previous updates to store.
# scale_inverse - if TRUE, scale the inverse Hessian approximation at each step.
lbfgs_direction <- function(memory = 5, scale_inverse = FALSE,
                            eps = .Machine$double.eps) {
  if (memory < 1) {
    stop("memory must be > 0")
  }
  make_direction(list(
    memory = memory,
    k = 0,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      opt$cache$gr_curr <- NULL
      opt$cache$gr_old <- NULL

      n <- length(par)
      sub_stage$value <- rep(0, n)

      sub_stage$rhos <- c()
      sub_stage$sms <- c()
      sub_stage$yms <- c()

      list(opt = opt, sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      gm_old <- opt$cache$gr_old

      if (is.null(gm_old)) {
        # First iteration do steepest descent unless we have home-brew
        # H approximation in fg
        pm <- lbfgs_guess(-gm, fg = fg, par = par)
      }
      else {
        # Update the memory
        # y_{k-1}, s_{k-1} using notation in Nocedal and Wright
        ym <- gm - gm_old
        sm <- opt$cache$update_old
        sub_stage <- lbfgs_memory_update(sub_stage, ym, sm, sub_stage$eps)

        # Solve Bp = -g with updated memory
        pm <- lbfgs_solve(-gm, sub_stage, scale_inverse, sub_stage$eps, fg, par)
      }

      descent <- dot(gm, pm)
      if (descent >= 0) {
        pm <- -gm
      }

      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
    , depends = c("gradient_old", "update_old")
  ))
}

# Newton Method -----------------------------------------------------------

# Newton method. Requires the Hessian to be calculated, via a function hs in fg.
newton_direction <- function(try_safe_chol = FALSE) {
  make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      opt$cache$gr_curr <- NULL
      list(opt = opt)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr

      if (!is.null(fg$hs)) {
        # B, an approximation to the (or the exact) Hessian
        # We now need to solve Bp = -g for p
        bm <- fg$hs(par)

        if (methods::is(bm, "matrix")) {
          if (try_safe_chol) {
            rm <- safe_chol(bm)
          }
          else {
            chol_result <- try({
              # O(N^3)
              rm <- chol(bm)
            },
            silent = TRUE)
            if (class(chol_result) == "try-error") {
              rm <- NULL
            }
          }

          if (is.null(rm)) {
            # message("Hessian is not positive-definite, resetting to SD")
            pm <- -gm
          }
          else {
            # Forward and back solving is "only" O(N^2)
            pm <- hessian_solve(rm, gm)
          }
        }
        else {
          # vector: assume it's a diagonal approximation of B
          pm <- -(1 / bm) * gm
        }
      }
      else if (!is.null(fg$hi)) {
        # H, an approximation to (or exact) inverse of the Hessian
        hm <- fg$hs(par)
        if (methods::is(hm, "matrix")) {
          pm <- -hm %*% gm
        }
        else {
          # vector: assume it's a diagonal approximation of H
          pm <- -hm * gm
        }
      }
      else {
        stop("No hi or hs function available for Hessian information")
      }

      descent <- dot(gm, pm)
      if (descent >= 0) {
        pm <- -gm
      }
      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
  ))
}

# A Partial Hessian approach: calculates the Cholesky decomposition of the
# Hessian (or some approximation) on the first iteration only. Future steps
# solve using this Hessian and the current gradient.
partial_hessian_direction <- function(hessian_every = 0) {
  make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      if (hessian_every == 0) {
        hm <- fg$hs(par)
        sub_stage$rm <- chol(hm)
      }
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      if (hessian_every > 0 && iter %% hessian_every == 0) {
        hm <- fg$hs(par)
        sub_stage$rm <- chol(hm)
      }

      gm <- opt$cache$gr_curr
      rm <- sub_stage$rm
      pm <- hessian_solve(rm, gm)

      descent <- dot(gm, pm)
      if (descent >= 0) {
        pm <- -gm
      }
      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
  ))
}

# Solves U'Ux = b
# U is upper triangular; need to solve U'(Ux) = b
# Step 1: Solve U'y = b which is a forwardsolve
# Step 2: Solve Ub = y which is a backsolve
# Can avoid explicit transpose in Step 1
# by passing the transpose argument to backsolve
upper_solve <- function(um, bm) {
  backsolve(um, backsolve(um, bm, transpose = TRUE))
}

# Given the upper triangular cholesky decomposition of the hessian (or an
# approximation), U, and the gradient vector, g, solve Up = -g
hessian_solve <- function(um, gm) {
  # This monkeying with dimensions is a tiny hack for the case of a repeated
  # block diagonal Hessian, so that you only need to provide the N x N block.
  nucol <- ncol(um)
  ngcol <- length(gm) / nucol
  dim(gm) <- c(nucol, ngcol)
  pm <- upper_solve(um, -gm)
  dim(pm) <- NULL
  pm
}


# Attempts to ensure a safe Cholesky decomposition of a Hessian by detecting
# a failure, rebuilding the Hessian by setting negative eigenvalues to a small
# positive value and then trying again. Not a fast procedure!
#
# Suggested by
# https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/
# Refs:
# Brissette, F. P., Khalili, M., & Leconte, R. (2007).
# Efficient stochastic generation of multi-site synthetic precipitation data.
# Journal of Hydrology, 345(3), 121-133.
# https://www.etsmtl.ca/getattachment/Unites-de-recherche/Drame/Publications/Brissette_al07---JH.pdf
#
# Rebonato, R., & Jaeckel, P. (2011).
# The most general methodology to create a valid correlation matrix for risk
# management and option pricing purposes.
# doi 10.21314/JOR.2000.023
safe_chol <- function(hm, eps = 1e-10) {
  rm <- NULL
  chol_result <- try({
    # O(N^3)
    rm <- chol(hm)
  },
  silent = TRUE)
  if (class(chol_result) == "try-error") {

    # Also O(N^3)
    eig <- eigen(hm)
    eig$values[eig$values < 0] <- 1e-10
    hm <- eig$vectors %*% (eig$values * diag(nrow(hm))) %*% t(eig$vectors)
    chol_result <- try({
      rm <- chol(hm)
    }, silent = TRUE)
  }
  rm
}


# Truncated Newton --------------------------------------------------------

# Truncated Newton, aka "Line Search Newton-CG"
# See Nocedal & Wright Chapter 7, algorithm 7.1
tn_direction <- function() {
  make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      inner_res <- tn_inner_cg(opt, fg, par, gm)
      opt <- inner_res$opt
      pm <- inner_res$zm

      descent <- dot(gm, pm)
      if (descent >= 0) {
        # message("TN: Not a descent direction")
        pm <- -gm
      }
      sub_stage$value <- pm
      list(opt = opt, sub_stage = sub_stage)
    }
  ))
}

# Linear CG inner loop of truncated Newton
# Inner loop from Nocedal & Wright Chapter 7, Algo 7.1
# Gives up when negative curvature or convergence occurs
# returns a list with opt (containing updated gradient count)
# and zm, the final estimate of pm solving Bp = -g (or -g)
tn_inner_cg <- function(opt, fg, par, gm) {
  gn <- norm2(gm)
  eps <- min(0.5, sqrt(gn)) * gn

  j <- 0
  zm <- 0
  rm <- gm
  dm <- -rm
  dot_rm <- dot(rm)

  # Safeguard for pathological situations.
  # In exact arithmetic CG will converge in N iterations.
  # The point of TN is to stop way earlier, but I've seen it get stuck,
  # presumably due to numerical issues.
  max_j <- length(par)
  while (j < max_j) {
    if (opt$counts$gr >= opt$convergence$max_gr) {
      if (j == 0) {
        zm <- -gm
      }
      break
    }

    Bd <- bd_approx(fg, par, dm, gm)
    opt$counts$gr <- opt$counts$gr + 1

    dBd <- dot(dm, Bd)
    if (dBd <= 0) {
      # -ve curvature, bail out
      if (j == 0) {
        # Use steepest descent if the initial Bd estimate is unusable
        zm <- -gm
      }
      break
    }

    alpha <- dot_rm / dBd
    zm <- zm + alpha * dm
    rm <- rm + alpha * Bd
    if (norm2(rm) < eps) {
      break
    }

    dot_rm_new <- dot(rm)
    beta <- dot_rm_new / dot_rm

    dm <- beta * dm - rm
    dot_rm <- dot_rm_new
    j <- j + 1
  }

  list(
    opt = opt,
    zm = zm
  )
}

# Finite difference approximation of the Hessian-vector product, Bd
# The default step size, h, is suggested by
# Andrei, N. (2009). Accelerated conjugate gradient algorithm with finite
# difference Hessian/vector product approximation for unconstrained
# optimization. Journal of Computational and Applied Mathematics, 230(2),
# 570-582.
# Found in:
# http://timvieira.github.io/blog/post/2014/02/10/gradient-vector-product/
# Something similar is used in minfunc
bd_approx <- function(fg, par, dm, gm,
                      h = 2 * sqrt(.Machine$double.eps) *
                        (1 + norm2(par)) / norm2(dm)) {
  g_fwd <- fg$gr(par + h * dm)
  (g_fwd - gm) / h
}

# Gradient Dependencies ------------------------------------------------------------

# Calculate the gradient at par.
require_gradient <- function(opt, stage, par, fg, iter) {
  if (!has_gr_curr(opt, iter)) {
    opt <- calc_gr_curr(opt, par, fg$gr, iter)

    if (any(!is.finite(opt$cache$gr_curr))) {
      opt$terminate <- list(
        what = "gr_inf",
        val = Inf
      )
      opt$is_terminated <- TRUE
    }
  }

  list(opt = opt)
}
attr(require_gradient, 'event') <- 'before gradient_descent'
attr(require_gradient, 'name') <- 'gradient'

# Caches the gradient at the current step.
require_gradient_old <- function(opt, par, fg, iter, par0, update) {
  opt$cache$gr_old <- opt$cache$gr_curr
  opt
}
attr(require_gradient_old, 'event') <- 'after step'
attr(require_gradient_old, 'name') <- 'gradient_old'
