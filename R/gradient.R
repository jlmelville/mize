# Gradient Direction -----------------------------------------------------------
make_direction <- function(sub_stage) {
  make_sub_stage(sub_stage, 'direction')
}

sd_direction <- function(normalize = FALSE) {

  make_direction(list(
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      #message("Calculating steepest descent direction")

      sub_stage$value <- -opt$cache$gr_curr

      if (sub_stage$normalize) {
        sub_stage$value <- normalize(sub_stage$value)
      }

      #message("sd pm = ", vec_formatC(sub_stage$value))

      list(sub_stage = sub_stage)
    },
    normalize = normalize
  ))
}


# Conjugate Gradient ------------------------------------------------------

# Wants: gradient, gradient_old, direction_old
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
        #message("beta = ", formatC(beta), " pm = ", vec_formatC(pm))
        descent <- dot(gm, pm)
        if (descent > 0) {
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

pr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  (dot(gm, gm) - dot(gm, gm_old)) / (dot(gm_old, gm_old) + eps)
}

pr_plus_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  beta <- pr_update(gm, gm_old, pm_old, eps)
#  if (beta < 0) {
#      message("PR+: restarting")
#  }
  max(0, beta)
}

fr_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  dot(gm, gm) / (dot(gm_old, gm_old) + eps)
}

hs_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  -(dot(gm, gm_old) - dot(gm, gm_old)) / (dot(pm_old, (gm - gm_old)) + eps)
}

dy_update <- function(gm, gm_old, pm_old, eps = .Machine$double.eps) {
  -dot(gm, gm) / (dot(pm_old, (gm - gm_old)) + eps)
}

cg_restart <- function(g_new, g_old, nu = 0.1) {
  # could only happen on first iteration
  if (is.null(g_old)) {
    return(TRUE)
  }
  ortho_test <- abs(dot(g_new, g_old)) / dot(g_new, g_new)
  should_restart <- ortho_test >= nu
  #  if (should_restart) {
  #    message("New CG direction not sufficiently orthogonal: ",
  #            formatC(ortho_test), " >= ", formatC(nu), " restarting")
  #  }
  should_restart
}


# BFGS --------------------------------------------------------------------

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
        if (descent > 0) {
          #message("BFGS direction is not a descent direction, resetting to SD")
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

lbfgs_direction <- function(memory = 100, scale_inverse = FALSE,
                            eps = .Machine$double.eps) {
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
          gamma <- dot(sm, ym) / dot(ym, ym)
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
        if (descent > 0) {
          #message("L-BFGS direction is not a descent direction, resetting to SD")
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

newton_direction <- function(scale_hess = FALSE) {
  make_direction(list(
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      gm <- opt$cache$gr_curr
      hm <- fg$hs(par)

      chol_result <- try({
        # O(N^3)
          rm <- chol(hm)
        },
        silent = TRUE)
      if (class(chol_result) == "try-error") {
        # Suggested by https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/
        # Refs:
        # FP Brissette, M Khalili, R Leconte, Journal of Hydrology, 2007, “Efficient stochastic generation of multi-site synthetic precipitation data”
        # https://www.etsmtl.ca/getattachment/Unites-de-recherche/Drame/Publications/Brissette_al07---JH.pdf
        # Rebonato, R., & Jäckel, P. (2011). The most general methodology to create a valid correlation matrix for risk management and option pricing purposes.
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

        #hm <- hm + (min(hm[hm > 0]) * 0.1 * diag(ncol(hm)))
        #rm <- chol(hm)
        # Forward and back solving is "only" O(N^2)
        ym <- forwardsolve(t(rm), -gm)
        pm <- backsolve(rm, ym)
        descent <- dot(gm, pm)
        if (descent > 0) {
          message("Newton direction is not a descent direction, resetting to SD")
          pm <- -gm
        }
      }
      sub_stage$value <- pm
      list(sub_stage = sub_stage)
    }
  ))
}

# Gradient Dependencies ------------------------------------------------------------

require_gradient <- function(opt, stage, par, fg, iter) {
  if (!has_gr_curr(opt, iter)) {
    #message("require gradient: calculating gr_curr ", iter)
    opt <- calc_gr_curr(opt, par, fg$gr, iter)

    if (any(is.nan(opt$cache$gr_curr))) {
      stop("NaN in grad. descent at iter ", iter)
    }
  }
  else {
    #message("require gradient: already have gr_curr")
  }
  list(opt = opt)
}
attr(require_gradient, 'event') <- 'before gradient_descent'
attr(require_gradient, 'name') <- 'gradient'

require_gradient_old <- function(opt, par, fg, iter, par0, update) {
  #  message("saving old gradient")
  opt$cache$gr_old <- opt$cache$gr_curr
  opt
}
attr(require_gradient_old, 'event') <- 'after step'
attr(require_gradient_old, 'name') <- 'gradient_old'
