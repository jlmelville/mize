# BFGS --------------------------------------------------------------------

# The Broyden Fletcher Goldfarb Shanno method
# scale_inverse - if TRUE, scale the inverse Hessian approximation on the first
#   step.
bfgs_direction <- function(eps = .Machine$double.eps, scale_inverse = FALSE) {
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
      } else {
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
        if (iter == 2 && scale_inverse && has_bfgs_curvature(ym, sm)) {
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
    },
    depends = c("gradient_old", "update_old")
  ))
}

# H - inverse Hessian
# s - difference in iterates x
# y - difference in gradient of x
bfgs_update <- function(hm, sm, ym, eps, eps_curv = 1e-8) {
  # Nocedal & Wright caution that skipping can discard useful curvature
  # information and prefer damped BFGS for this situation. This guard is a
  # narrow fallback for steps that do not satisfy Wolfe curvature.
  if (!has_bfgs_curvature(ym, sm, eps_curv)) {
    return(hm)
  }

  rho <- 1 / (dot(ym, sm) + eps)
  im <- diag(1, nrow(hm))

  rss <- rho * outer(sm, sm)
  irsy <- im - rho * outer(sm, ym)
  irys <- im - rho * outer(ym, sm)

  (irsy %*% (hm %*% irys)) + rss
}

# Check if the curvature condition is satisfied for BFGS updates.
has_bfgs_curvature <- function(ym, sm, eps_curv = 1e-8) {
  ys <- dot(ym, sm)
  curvature_scale <- eps_curv * norm2(ym) * norm2(sm)

  is.finite(ys) && is.finite(curvature_scale) && ys > curvature_scale
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
sr1_direction <- function(
  eps = .Machine$double.eps,
  scale_inverse = FALSE,
  skip_bad_update = TRUE,
  tol = 1e-8,
  try_bfgs = TRUE
) {
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
      } else {
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
        } else {
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
          } else {
            # Only store the BFGS H for next iteration if it's a descent
            sub_stage$hm <- hm_bfgs
          }
        } else {
          pm <- -gm
        }
      }

      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    },
    depends = c("gradient_old", "update_old")
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
lbfgs_guess <- function(
  qm,
  scale_inverse = TRUE,
  rho = NULL,
  ym = NULL,
  eps = NULL,
  fg = NULL,
  par = NULL
) {
  # User-defined hessian inverse approximations
  if (!is.null(fg) && !is.null(fg$hi)) {
    hm <- fg$hi(par)
    if (methods::is(hm, "matrix")) {
      # Full matrix: not necessarily a great idea for memory usage
      pm <- hm %*% qm
    } else {
      # It's a vector representing a diagonal matrix
      pm <- hm * qm
    }
  } else {
    # The usual L-BFGS guess
    if (!is.null(rho) && !is.null(ym) && scale_inverse) {
      # Eqn 7.20 in Nocedal & Wright
      gamma <- 1 / (rho * (dot(ym) + eps))
      hm <- rep(gamma, length(ym))
      pm <- hm * qm
    } else {
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
lbfgs_solve <- function(
  qm,
  lbfgs_state,
  scale_inverse,
  eps,
  fg = NULL,
  par = NULL
) {
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
  pm <- lbfgs_guess(
    qm,
    scale_inverse,
    rhos[[length(rhos)]],
    yms[[length(yms)]],
    eps,
    fg,
    par
  )

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
lbfgs_memory_update <- function(lbfgs_state, ym, sm, eps, eps_curv = 1e-8) {
  if (!has_bfgs_curvature(ym, sm, eps_curv)) {
    return(lbfgs_state)
  }

  rho <- 1 / (dot(ym, sm) + eps)

  rhos <- lbfgs_state$rhos
  sms <- lbfgs_state$sms
  yms <- lbfgs_state$yms

  # append latest values to memory
  sms <- c(sms, list(sm))
  yms <- c(yms, list(ym))
  rhos <- c(rhos, list(rho))

  # discard oldest values if we've exceeded the memory limit
  if (length(sms) > lbfgs_state$memory) {
    keep <- seq.int(length(sms) - lbfgs_state$memory + 1L, length(sms))
    sms <- sms[keep]
    yms <- yms[keep]
    rhos <- rhos[keep]
  }

  lbfgs_state$sms <- sms
  lbfgs_state$yms <- yms
  lbfgs_state$rhos <- rhos

  lbfgs_state
}

# The Limited Memory BFGS method
#
# memory - The number of previous updates to store.
# scale_inverse - if TRUE, scale the inverse Hessian approximation at each step.
lbfgs_direction <- function(
  memory = 5,
  scale_inverse = FALSE,
  eps = .Machine$double.eps
) {
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
      } else {
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
    },
    depends = c("gradient_old", "update_old")
  ))
}
