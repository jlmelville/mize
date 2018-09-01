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
                         preconditioner = "", memory = 5,
                         eps = .Machine$double.eps) {
  cg <- make_direction(list(
    ortho_check = ortho_check,
    nu = nu,
    cg_update = cg_update,
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      if (preconditioner == "l-bfgs") {
        res <- lbfgs_init(opt, stage, sub_stage, par, fg, iter)
        res$sub_stage$memory <- memory
        res$sub_stage$eps <- eps
        sub_stage$preconditioner <- res$sub_stage
      }
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
      # wm = Pg or just g if we're not preconditioning
      wm <- gm
      if (!is.null(gm_old)
          && (!sub_stage$ortho_check
              || !sub_stage$cg_restart(gm, gm_old, sub_stage$nu))) {

        precondition_fn <- NULL
        if (preconditioner == "l-bfgs" && !is.null(opt$cache$gr_old)) {
          lbfgs <- sub_stage$preconditioner
          ym <- gm - opt$cache$gr_old
          sm <- opt$cache$update_old
          lbfgs <- lbfgs_memory_update(lbfgs, ym, sm, lbfgs$eps)
          precondition_fn <- function(rm) {
            lbfgs_solve(rm, lbfgs, scale_inverse = TRUE, eps = lbfgs$eps)
          }
          wm <- precondition_fn(gm)
          # p <- -Pg
          pm <- -wm

          sub_stage$preconditioner <- lbfgs
        }
        beta <- sub_stage$cg_update(gm, gm_old, pm_old, sub_stage$eps,
                                    wm = wm,
                                    preconditioner = precondition_fn)
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
  )

  )
  if (preconditioner == "l-bfgs") {
    cg$depends <- append(cg$depends, c("update_old"))
  }
  cg
}


# CG update formulae, grouped according to their numerators similar to the
# discussion in Hager and Zhang's survey paper
#
# The FR, CD and DY updates are all susceptible to "jamming": they can end up
# with very small step sizes and make little progress.
#
# Preconditioned expressions are based off those given by the Hager-Zhang
# CG survey (2006), L-CG_DESCENT paper (2013) and damped CG paper by
# Al-Baali, Caliciotti, Fasano and Roma (2017), except for preconditioned
# LS, DY and PR-FR which I haven't seen explicitly written down, but which
# it is easy to derive from the other examples.

# Fletcher-Reeves update.
fr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL,
                      preconditioner = NULL) {
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
    wm_old <- preconditioner(gm_old)
  }
  else {
    if (is.null(wm)) {
      wm <- gm
    }
    wm_old <- gm_old
  }
  dot(gm, wm) / (dot(gm_old, wm_old) + eps)
}

# Conjugate Descent update due to Fletcher.
cd_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL,
                      preconditioner = NULL) {
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  }
  else {
    if (is.null(wm)) {
      wm <- gm
    }
  }
  dot(-gm, wm) / (dot(pm_old, gm_old) + eps)
}

# The Dai-Yuan update.
dy_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL,
                      preconditioner = NULL) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  }
  else {
    if (is.null(wm)) {
      wm <- gm
    }
  }
  dot(gm, wm) / (dot(pm_old, ym) + eps)
}

# HS, PR and LS share a numerator. According to Hager and Zhang, they
# perform better in practice than the FR, CD and DY updates, despite less
# being known about their provable global convergence properties.

# The Hestenes-Stiefel update.
hs_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL,
                      preconditioner = NULL) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  }
  else {
    wm <- gm
  }
  dot(ym, wm) / (dot(pm_old, ym) + eps)
}

# An "HS+" modification of Hestenes-Stiefel, in analogy to the "PR+" variant of
# Polak-Ribiere suggested by Powell. As far as I can tell, Hager and Zhang
# suggested this modification.
# Hager, W. W., & Zhang, H. (2006).
# A survey of nonlinear conjugate gradient methods.
# \emph{Pacific journal of Optimization}, \emph{2}(1), 35-58.
hs_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                           wm = NULL,
                           preconditioner = NULL) {
  beta <- hs_update(gm, gm_old, pm_old, eps,
                    wm = wm, preconditioner = preconditioner)
  max(0, beta)
}

# The Polak-Ribiere method for updating the CG direction. Also known as
# Polak-Ribiere-Polyak (PRP)
pr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL,
                      preconditioner = NULL) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
    wm_old <- preconditioner(gm_old)
  }
  else {
    if (is.null(wm)) {
      wm <- gm
    }
    wm_old <- gm_old
  }
  dot(wm, ym) / (dot(gm_old, wm_old) + eps)
}

# The "PR+" update due to Powell. Polak-Ribiere update, but if negative,
# restarts the CG from steepest descent. Prevents a possible lack of
# convergence when using a Wolfe line search.
pr_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                           wm = NULL,
                           preconditioner = NULL) {
  beta <- pr_update(gm, gm_old, pm_old, eps,
                    wm = wm, preconditioner = preconditioner)
  max(0, beta)
}

# Liu-Storey update
ls_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL, preconditioner = NULL) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  }
  else {
    if (is.null(wm)) {
      wm <- gm
    }
  }
  dot(-ym, wm) / (dot(pm_old, gm_old) + eps)
}

# Hager-Zhang update as used in CG_DESCENT, theta = 2
hz_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                      wm = NULL,
                      preconditioner = NULL) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    vm <- preconditioner(ym)
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  }
  else {
    vm <- ym
    if (is.null(wm)) {
      wm <- gm
    }
  }
  ipy <- 1 / (dot(pm_old, ym) + eps)
  (dot(ym, wm) * ipy) - 2 * dot(ym, vm) * ipy * dot(pm_old, gm) * ipy
}

# "Restricted" Hager-Zhang update as used in CG_DESCENT to ensure
# convergence. Analogous to the PR+ and HS+ updates, but dynamically adjusts
# the lower bound as convergence occurs. Choice of eta is from the CG_DESCENT
# paper
hz_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                           wm = wm, preconditioner = NULL) {
  beta <- hz_update(gm, gm_old, pm_old, eps,
                    wm = wm, preconditioner = preconditioner)
  eta <- 0.01
  eta_k <- -1 / (dot(pm_old) * min(eta, dot(gm_old)))
  max(eta_k, beta)
}

# The PR-FR update suggested by Gilbert and Nocedal (1992)
prfr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps,
                        wm = wm, preconditioner = NULL) {
  bpr <- pr_update(gm, gm_old, pm_old, eps = eps,
                   wm = wm, preconditioner = preconditioner)
  bfr <- fr_update(gm, gm_old, pm_old, eps = eps,
                   wm = wm, preconditioner = preconditioner)
  if (bpr < -bfr) {
    beta <- -bfr
  }
  else if (abs(bpr) <= bfr) {
    beta <- bpr
  }
  else if (bpr > bfr) {
    beta <- bfr
  }
  else {
    stop("Problem in PR-FR Update")
  }
  beta
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

lbfgs_init <- function(opt, stage, sub_stage, par, fg, iter) {
  opt$cache$gr_curr <- NULL
  opt$cache$gr_old <- NULL

  n <- length(par)
  sub_stage$value <- rep(0, n)

  sub_stage$rhos <- c()
  sub_stage$sms <- c()
  sub_stage$yms <- c()

  list(opt = opt, sub_stage = sub_stage)
}
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
      hm <- rep(gamma, length(ym))
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
# Note that we usually want to find p = -Hg, so usually take the negative of
# the return value
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
    init = lbfgs_init,
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
        # lbfgs_solve returns Hg, so take negative of result
        pm <- -lbfgs_solve(gm, sub_stage, scale_inverse, sub_stage$eps, fg, par)
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
# L-BFGS preconditioner suggested e.g. by Hager & Zhang (2006)
# also implemented in minfunc
tn_direction <- function(init = 0, exit_criterion = "curvature",
                         preconditioner = "", memory = 5,
                         eps = .Machine$double.eps) {
  tn <- make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      if (preconditioner == "l-bfgs") {
        res <- lbfgs_init(opt, stage, sub_stage, par, fg, iter)
        res$sub_stage$memory <- memory
        res$sub_stage$eps <- eps
        sub_stage$preconditioner <- res$sub_stage
      }

      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr

      precondition_fn <- NULL
      # In the context of non-linear CG, Hager and Zhang (2006) mention
      # using the approximation of the inverse Hessian as a preconditioner P
      # replacing g with Pg and using L-BFGS specifically.
      # In linear CG, use r instead of g
      if (preconditioner == "l-bfgs" && !is.null(opt$cache$gr_old)) {
        lbfgs <- sub_stage$preconditioner
        ym <- gm - opt$cache$gr_old
        sm <- opt$cache$update_old
        lbfgs <- lbfgs_memory_update(lbfgs, ym, sm, lbfgs$eps)
        precondition_fn <- function(rm) {
          # This solves Bp = r for p => p = Hr which is what we want,
          # so don't take negative (as needed for L-BFGS solution of pm)
          lbfgs_solve(rm, lbfgs, scale_inverse = TRUE, eps = lbfgs$eps)
        }

        sub_stage$preconditioner <- lbfgs
      }

      if (init == 0) {
        # Standard initialization. Saves on a finite difference gradient
        # calculation because we know the residual is g
        zm <- 0
      }
      else if (init == "l-bfgs" && preconditioner == "l-bfgs") {
        # Use the L-BFGS guess as we've done a lot of the work already
        # A potentially bad idea of my own invention
        if (!is.null(opt$cache$gr_old)) {
          lbfgs <- sub_stage$preconditioner
          zm <- -lbfgs_solve(gm, lbfgs, scale_inverse = TRUE, eps = lbfgs$eps)
        }
        else {
          zm <- 0
        }
      }
      else {
        # "previous", i.e. use the result from the last iteration
        # Martens (2010)
        if (!is.null(sub_stage$value)) {
          zm <- sub_stage$value
        }
        else {
          zm <- 0
        }
      }
      inner_res <- tn_inner_cg(opt, fg, par, gm, precondition_fn, zm = zm,
                               exit_criterion = exit_criterion)
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

  if (preconditioner == "l-bfgs") {
    tn$depends <- append(tn$depends, c("gradient_old", "update_old"))
  }

  tn
}

# Linear CG inner loop of truncated Newton
# Inner loop from Nocedal & Wright Chapter 7, Algo 7.1
# with modification for preconditioner as in Chapter 5, Algo 5.3
# Gives up when negative curvature or convergence occurs
# returns a list with opt (containing updated gradient count)
# and zm, the final estimate of pm solving Bp = -g
# max_iter Maximum number of iterations for the inner CG loop. With exact
#  arithmetic, we would expect the loop to terminate in N steps at most,
#  but rounding errors can prevent that and we want to stop way before that
#  anyway. Default is from TNPACK (Xie and Schlick 1999).
# exit_criterion: either "curvature" (check negative curvature) or "descent",
#  a comparison of descent direction used in TNPACK (Xie and Schlick 1999)
tn_inner_cg <- function(opt, fg, par, gm, preconditioner = NULL, zm = 0,
                        exit_criterion = "curvature", max_iter = 40) {
  gn <- norm2(gm)
  eps <- min(0.5, sqrt(gn)) * gn

  if (length(zm) == 1 && zm == 0) {
    rm <- gm
  }
  else {
    if (opt$counts$gr >= opt$convergence$max_gr) {
      zm <- -gm
      return(list(
        opt = opt,
        zm = zm)
      )
    }
    Bd <- bd_approx(fg, par, zm, gm)
    opt$counts$gr <- opt$counts$gr + 1

    rm <- Bd + gm
  }
  if (!is.null(preconditioner)) {
    ym <- preconditioner(rm)
  }
  else {
    ym <- rm
  }

  dm <- -ym
  dot_rm_ym <- dot(rm, ym)

  j <- 0
  while (j < max_iter) {
    if (opt$counts$gr >= opt$convergence$max_gr) {
      if (j == 0) {
        zm <- -gm
      }
      return(list(
        opt = opt,
        zm = zm)
      )
    }

    Bd <- bd_approx(fg, par, dm, gm)
    opt$counts$gr <- opt$counts$gr + 1

    dBd <- dot(dm, Bd)
    if (exit_criterion == "curvature") {
      if (dBd <= 0) {
        # -ve curvature, bail out
        if (j == 0) {
          # Use steepest descent if the initial Bd estimate is unusable
          zm <- -gm
        }
        break
      }
    }

    alpha <- dot_rm_ym / dBd
    # only used in the strong curvature criterion
    zm_old <- zm
    zm <- zm + alpha * dm
    if (exit_criterion == "strong") {
      # alternative exit criterion used in TNPACK
      if (dot(gm, zm) >= dot(gm, zm_old) + 1e-15) {
        if (j == 0) {
          zm <- -gm
        }
        else {
          zm <- zm_old
        }

        break
      }
    }

    rm <- rm + alpha * Bd
    if (norm2(rm) < eps) {
      break
    }

    if (!is.null(preconditioner)) {
      ym <- preconditioner(rm)
    }
    else {
      ym <- rm
    }

    dot_rm_ym_new <- dot(rm, ym)
    beta <- dot_rm_ym_new / dot_rm_ym

    dm <- beta * dm - ym
    dot_rm_ym <- dot_rm_ym_new
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
