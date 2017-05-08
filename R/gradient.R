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
      list(sub_stage = sub_stage)
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
  dot(gm, gm) / (dot(gm_old, gm_old) + eps)
}

# Conjugate Descent update due to Fletcher
cd_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  dot(gm, gm) / (dot(pm_old, (gm - gm_old)) + eps)
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
  dot(gm, gm - gm_old) / (dot(gm_old, gm_old) + eps)
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
  dot(ym - 2 * pm_old * (dot(ym, ym) / (py + eps)), (gm / (py + eps)))
}

# "Restricted" Hager-Zhang update as used in CG_DESCENT to ensure
# convergence. Analogous to the PR+ and HS+ updates, but dynamically adjusts
# the lower bound as convergence occurs. Choice of eta is from the CG_DESCENT
# paper
hz_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  beta <- hz_update(gm, gm_old, pm_old, eps)
  eta <- 0.01
  eta_k <- -1 / (dot(pm_old, pm_old) * min(eta, dot(gm_old, gm_old)))
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
  ortho_test <- abs(dot(g_new, g_old)) / dot(g_new, g_new)
  should_restart <- ortho_test >= nu
  should_restart
}


# BFGS --------------------------------------------------------------------

# The Broyden Fletcher Goldfarb Shanno method
# scale_inverse - if TRUE, scale the inverse Hessian approximation on the first
#   step.
bfgs_direction <- function(eps =  .Machine$double.eps,
                           scale_inverse = FALSE) {
  make_direction(list(
    eps = eps,
    init = function(opt, stage, sub_stage, par, fg, iter) {
      n <- length(par)
      sub_stage$value <- rep(0, n)
      sub_stage$hm <- diag(1, n)
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      gm_old <- opt$cache$gr_old
      if (is.null(gm_old)) {
        pm <- -gm
      }
      else {
        sm <- opt$cache$update_old
        hm <- sub_stage$hm

        ym <- gm - gm_old

        if (iter == 2 && scale_inverse) {
          # Nocedal suggests this heuristic for scaling the first
          # approximation in Chapter 6 Section "Implementation"
          # Also used in the definition of L-BFGS
          gamma <- dot(sm, ym) / dot(ym, ym)
          hm <- gamma * hm
        }

        rho <- 1 / (dot(ym, sm) + sub_stage$eps)
        im <- diag(1, nrow(hm))

        rss <- rho * outer(sm, sm)
        irsy <- im - rho * outer(sm, ym)
        irys <- im - rho * outer(ym, sm)

        sub_stage$hm <- (irsy %*% (hm %*% irys)) + rss

        pm <- as.vector(-sub_stage$hm %*% gm)

        descent <- dot(gm, pm)
        if (descent >= 0) {
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
      n <- length(par)
      sub_stage$value <- rep(0, n)

      sub_stage$rhos <- c()
      sub_stage$sms <- c()
      sub_stage$yms <- c()

      list(sub_stage = sub_stage)

    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      gm_old <- opt$cache$gr_old

      if (is.null(gm_old)) {
        pm <- -gm
      }
      else {
        rhos <- sub_stage$rhos
        sms <- sub_stage$sms
        yms <- sub_stage$yms

        # discard oldest values if we've reached memory limit
        if (length(sms) == sub_stage$memory) {
          sms <- sms[2:length(sms)]
          yms <- yms[2:length(yms)]
          rhos <- rhos[2:length(rhos)]
        }

        # y_{k-1}, s_{k-1}, rho_{k-1} using notation in Nocedal
        ym <- gm - gm_old
        sm <- opt$cache$update_old
        rho <- 1 / (dot(ym, sm) + sub_stage$eps)

        # append latest values to memory
        sms <- c(sms, list(sm))
        yms <- c(yms, list(ym))
        rhos <- c(rhos, list(rho))

        qm <- gm
        alphas <- rep(0, length(rhos))
        # loop backwards latest values first
        for (i in length(rhos):1) {
          alphas[i] <- rhos[[i]] * dot(sms[[i]], qm)
          qm <- qm - alphas[[i]] * yms[[i]]
        }

        if (scale_inverse) {
          gamma <- dot(sm, ym) / (dot(ym, ym) + sub_stage$eps)
        }
        else {
          gamma <- 1
        }

        hm <- rep(gamma, length(par))
        pm <- hm * qm
        # loop forwards
        for (i in 1:length(rhos)) {
          beta <- rhos[[i]] * dot(yms[[i]], pm)
          pm <- pm + sms[[i]] * (alphas[[i]] - beta)
        }
        pm <- -pm

        descent <- dot(gm, pm)
        if (descent >= 0) {
          pm <- -gm
        }

        sub_stage$value <- pm
        sub_stage$sms <- sms
        sub_stage$yms <- yms
        sub_stage$rhos <- rhos

      }
      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
    , depends = c("gradient_old", "update_old")
  ))
}


# Newton Method -----------------------------------------------------------

# Newton method. Requires the Hessian to be calculated, via a function hs in fg.
newton_direction <- function() {
  make_direction(list(
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      if (is.null(fg$hs)) {
        stop("No Hessian function 'hs', defined for fg")
      }
      hm <- fg$hs(par)

      chol_result <- try({
        # O(N^3)
          rm <- chol(hm)
        },
        silent = TRUE)
      if (class(chol_result) == "try-error") {
        # Suggested by https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/
        # Refs:
        # FP Brissette, M Khalili, R Leconte, Journal of Hydrology, 2007,
        # Efficient stochastic generation of multi-site synthetic precipitation data
        # https://www.etsmtl.ca/getattachment/Unites-de-recherche/Drame/Publications/Brissette_al07---JH.pdf
        # Rebonato, R., & Jäckel, P. (2011).
        # The most general methodology to create a valid correlation matrix for risk management and option pricing purposes.
        # doi 10.21314/JOR.2000.023
        # Also O(N^3)
        eig <- eigen(hm)
        eig$values[eig$values < 0] <- 1e-10
        hm <- eig$vectors %*% (eig$values * diag(nrow(hm))) %*% t(eig$vectors)
        chol_result <- try({
          rm <- chol(hm)
        }, silent = TRUE)
      }
      if (class(chol_result) == "try-error") {
        # we gave it a good go, but let's just do steepest descent this time
        #message("Hessian is not positive-definite, resetting to SD")
        pm <- -gm
      }
      else {
        # Forward and back solving is "only" O(N^2)
        pm <- hessian_solve(rm, gm)

        descent <- dot(gm, pm)
        if (descent >= 0) {
          pm <- -gm
        }
      }
      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
  ))
}

# A Partial Hessian approach: calculates the Cholesky decomposition of the
# Hessian (or some approximation) on the first iteration only. Future steps
# solve using this Hessian and the current gradient.
partial_hessian_direction <- function() {
  make_direction(list(
    init = function(opt, stage, sub_stage, par, fg, iter) {
      hm <- fg$hs(par)
      sub_stage$rm <- chol(hm)
      list(sub_stage = sub_stage)
    },
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
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
  nucol <- ncol(um)
  ngcol <- length(gm) / nucol
  dim(gm) <- c(nucol, ngcol)
  pm <- upper_solve(um, -gm)
  dim(pm) <- NULL
  pm
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
