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
cg_direction <- function(
  ortho_check = FALSE,
  nu = 0.1,
  cg_update = pr_plus_update,
  preconditioner = "",
  memory = 5,
  eps = .Machine$double.eps
) {
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
      if (
        !is.null(gm_old) &&
          (!sub_stage$ortho_check ||
            !sub_stage$cg_restart(gm, gm_old, sub_stage$nu))
      ) {
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
        beta <- sub_stage$cg_update(
          gm,
          gm_old,
          pm_old,
          sub_stage$eps,
          wm = wm,
          preconditioner = precondition_fn
        )
        pm <- pm + (beta * pm_old)
        descent <- dot(gm, pm)
        if (descent >= 0) {
          # message("Next CG direction is not a descent direction, resetting to SD")
          pm <- -gm
        }
      }

      sub_stage$value <- pm

      list(sub_stage = sub_stage)
    },
    after_step = function(opt, stage, sub_stage, par, fg, iter, par0, update) {
      sub_stage$pm_old <- sub_stage$value
      list(sub_stage = sub_stage)
    },
    depends = c("gradient_old")
  ))
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
fr_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
    wm_old <- preconditioner(gm_old)
  } else {
    if (is.null(wm)) {
      wm <- gm
    }
    wm_old <- gm_old
  }
  dot(gm, wm) / (dot(gm_old, wm_old) + eps)
}

# Conjugate Descent update due to Fletcher.
cd_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  } else {
    if (is.null(wm)) {
      wm <- gm
    }
  }
  dot(-gm, wm) / (dot(pm_old, gm_old) + eps)
}

# The Dai-Yuan update.
dy_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  } else {
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
hs_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  } else {
    wm <- gm
  }
  dot(ym, wm) / (dot(pm_old, ym) + eps)
}

# An "HS+" modification of Hestenes-Stiefel, in analogy to the "PR+" variant of
# Polak-Ribiere suggested by Powell. As far as I can tell, Hager and Zhang
# suggested this modification.
# Hager, W. W., & Zhang, H. (2006).
# A survey of nonlinear conjugate gradient methods.
# *Pacific journal of Optimization*, *2*(1), 35-58.
hs_plus_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  beta <- hs_update(
    gm,
    gm_old,
    pm_old,
    eps,
    wm = wm,
    preconditioner = preconditioner
  )
  max(0, beta)
}

# The Polak-Ribiere method for updating the CG direction. Also known as
# Polak-Ribiere-Polyak (PRP)
pr_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
    wm_old <- preconditioner(gm_old)
  } else {
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
pr_plus_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  beta <- pr_update(
    gm,
    gm_old,
    pm_old,
    eps,
    wm = wm,
    preconditioner = preconditioner
  )
  max(0, beta)
}

# Liu-Storey update
ls_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  } else {
    if (is.null(wm)) {
      wm <- gm
    }
  }
  dot(-ym, wm) / (dot(pm_old, gm_old) + eps)
}

# Hager-Zhang update as used in CG_DESCENT, theta = 2
hz_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  ym <- gm - gm_old
  if (!is.null(preconditioner)) {
    vm <- preconditioner(ym)
    if (is.null(wm)) {
      wm <- preconditioner(gm)
    }
  } else {
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
hz_plus_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  beta <- hz_update(
    gm,
    gm_old,
    pm_old,
    eps,
    wm = wm,
    preconditioner = preconditioner
  )
  eta <- 0.01
  eta_k <- -1 / (dot(pm_old) * min(eta, dot(gm_old)))
  max(eta_k, beta)
}

# The PR-FR update suggested by Gilbert and Nocedal (1992)
prfr_update <- function(
  gm,
  gm_old,
  pm_old,
  eps = .Machine$double.eps,
  wm = NULL,
  preconditioner = NULL
) {
  bpr <- pr_update(
    gm,
    gm_old,
    pm_old,
    eps = eps,
    wm = wm,
    preconditioner = preconditioner
  )
  bfr <- fr_update(
    gm,
    gm_old,
    pm_old,
    eps = eps,
    wm = wm,
    preconditioner = preconditioner
  )
  if (bpr < -bfr) {
    beta <- -bfr
  } else if (abs(bpr) <= bfr) {
    beta <- bpr
  } else if (bpr > bfr) {
    beta <- bfr
  } else {
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
