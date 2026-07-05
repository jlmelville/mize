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
          } else {
            chol_result <- try(
              {
                # O(N^3)
                rm <- chol(bm)
              },
              silent = TRUE
            )
            if (methods::is(chol_result, "try-error")) {
              rm <- NULL
            }
          }

          if (is.null(rm)) {
            # message("Hessian is not positive-definite, resetting to SD")
            pm <- -gm
          } else {
            # Forward and back solving is "only" O(N^2)
            pm <- hessian_solve(rm, gm)
          }
        } else {
          # vector: assume it's a diagonal approximation of B
          pm <- -(1 / bm) * gm
        }
      } else if (!is.null(fg$hi)) {
        # H, an approximation to (or exact) inverse of the Hessian
        hm <- fg$hi(par)
        if (methods::is(hm, "matrix")) {
          pm <- -hm %*% gm
        } else {
          # vector: assume it's a diagonal approximation of H
          pm <- -hm * gm
        }
      } else {
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
  chol_result <- try(
    {
      # O(N^3)
      rm <- chol(hm)
    },
    silent = TRUE
  )
  if (methods::is(chol_result, "try-error")) {
    # Also O(N^3)
    eig <- eigen(hm)
    eig$values[eig$values < 0] <- 1e-10
    hm <- eig$vectors %*% (eig$values * diag(nrow(hm))) %*% t(eig$vectors)
    chol_result <- try(
      {
        rm <- chol(hm)
      },
      silent = TRUE
    )
  }
  rm
}


# Truncated Newton --------------------------------------------------------

# Truncated Newton, aka "Line Search Newton-CG"
# See Nocedal & Wright Chapter 7, algorithm 7.1
# L-BFGS preconditioner suggested e.g. by Hager & Zhang (2006)
# also implemented in minfunc
tn_direction <- function(
  init = 0,
  exit_criterion = "curvature",
  preconditioner = "",
  memory = 5,
  eps = .Machine$double.eps
) {
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
      } else if (init == "l-bfgs" && preconditioner == "l-bfgs") {
        # Use the L-BFGS guess as we've done a lot of the work already
        # A potentially bad idea of my own invention
        if (!is.null(opt$cache$gr_old)) {
          lbfgs <- sub_stage$preconditioner
          zm <- -lbfgs_solve(gm, lbfgs, scale_inverse = TRUE, eps = lbfgs$eps)
        } else {
          zm <- 0
        }
      } else {
        # "previous", i.e. use the result from the last iteration
        # Martens (2010)
        if (!is.null(sub_stage$value)) {
          zm <- sub_stage$value
        } else {
          zm <- 0
        }
      }
      inner_res <- tn_inner_cg(
        opt,
        fg,
        par,
        gm,
        precondition_fn,
        zm = zm,
        exit_criterion = exit_criterion
      )
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
tn_inner_cg <- function(
  opt,
  fg,
  par,
  gm,
  preconditioner = NULL,
  zm = 0,
  exit_criterion = "curvature",
  max_iter = 40
) {
  gn <- norm2(gm)
  eps <- min(0.5, sqrt(gn)) * gn

  if (length(zm) == 1 && zm == 0) {
    rm <- gm
  } else {
    if (opt$counts$gr >= opt$convergence$max_gr) {
      zm <- -gm
      return(list(
        opt = opt,
        zm = zm
      ))
    }
    Bd <- bd_approx(fg, par, zm, gm)
    opt$counts$gr <- opt$counts$gr + 1

    rm <- Bd + gm
  }
  if (!is.null(preconditioner)) {
    ym <- preconditioner(rm)
  } else {
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
        zm = zm
      ))
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
        } else {
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
    } else {
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
bd_approx <- function(
  fg,
  par,
  dm,
  gm,
  h = 2 * sqrt(.Machine$double.eps) * (1 + norm2(par)) / norm2(dm)
) {
  g_fwd <- fg$gr(par + h * dm)
  (g_fwd - gm) / h
}
