# Line Search as described by Hager and Zhang in:
#
# Hager, W. W., & Zhang, H. (2005).
# A new conjugate gradient method with guaranteed descent and an efficient line
# search.
# SIAM Journal on Optimization, 16(1), 170-192.
#
# Hager, W. W., & Zhang, H. (2006).
# Algorithm 851: CG_DESCENT, a conjugate gradient method with guaranteed
# descent.
# ACM Transactions on Mathematical Software (TOMS), 32(1), 113-137.
#
# I have tried to indicate in the comments which parts of which routines
# match the notation in the above papers.


# Adapter -----------------------------------------------------------------

hager_zhang <- function(c1 = c2 / 2, c2 = 0.1, max_fn = Inf,
                        strong_curvature = TRUE,
                        approx_armijo = TRUE) {
  if (c2 < c1) {
    stop("hager-zhang line search: c2 < c1")
  }
  function(phi, step0, alpha,
           total_max_fn = Inf, total_max_gr = Inf, total_max_fg = Inf,
           pm) {
    maxfev <- min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    if (maxfev <= 0) {
      return(list(step = step0, nfn = 0, ngr = 0))
    }
    res <- line_search_hz(alpha, step0, phi, c1 = c1, c2 = c2,
                               eps = 1e-6, theta = 0.5, rho = 5,
                               gamma = 0.66,
                               max_fn = maxfev,
                               xtol = 1e-6,
                               strong_curvature = strong_curvature,
                               always_check_convergence = TRUE,
                               approx_armijo = approx_armijo,
                               verbose = FALSE)

    res$ngr = res$nfn
    res
  }
}


# Line Search -------------------------------------------------------------

# Routine 'Line Search Algorithm' L0-3

# alpha initial guess
# step0
# phi
# c1
# c2
# eps - determines when to apply approximate Wolfe condition. Ignored if
# approx_armijo is FALSE.
# theta - bisection weight during bracket update
# rho - factor to increase step size by during bracket phase
# gamma - bisection weight if secant2 step fails
# max_fn - maximum number of function evaluations allowed
# xtol - stop if bracket size ever falls below alpha * xtol
# strong_curvature - use Strong curvature condition
# always_check_convergence - if FALSE, only check for satisfaction of the
#   Wolfe conditions when a step size is produced by interpolation, not
#   bisection or the bracketing phase. May produce a result closer to a
#   minimizer.
# approx_armijo - if TRUE, use the approximate version of the armijo condition
# when testing the Wolfe conditions.
line_search_hz <- function(alpha, step0, phi, c1 = 0.1, c2 = 0.9,
                           eps = 1e-6, theta = 0.5, rho = 5,
                           gamma = 0.66,
                           max_fn = Inf,
                           xtol = 1e-6,
                           strong_curvature = FALSE,
                           always_check_convergence = TRUE,
                           approx_armijo = TRUE,
                           verbose = FALSE) {
  ls_max_fn <- max_fn
  if (max_fn == 0) {
    return(list(step = step0, nfn = 0))
  }
  # For approximate Wolfe error parameters implicitly set delta = 0
  # no C or Q parameters to update
  eps_k <- eps * abs(step0$f)

  # Bisect alpha if needed in case initial alpha guess gives a mad or bad value
  nfn <- 0
  result <- find_finite(phi, alpha, max_fn, min_alpha = 0)
  nfn <- nfn + result$nfn
  if (!result$ok) {
    if (verbose) {
      message("Unable to create finite initial guess")
    }
    return(list(step = step0, nfn = nfn))
  }
  step_c <- result$step

  if (always_check_convergence && hz_ok_step(step_c, step0, c1, c2, eps_k,
                                             strong_curvature = strong_curvature,
                                             approx_armijo = approx_armijo)) {
    if (verbose) {
      message("initial step OK = ", formatC(step_c$alpha))
    }
    return(list(step = step_c, nfn = nfn))
  }
  ls_max_fn <- max_fn - nfn
  if (ls_max_fn <= 0) {
    if (verbose) {
      message("max fn reached after initial step")
    }
    return(list(step = step_c, nfn = nfn))
  }

  # L0
  br_res <- bracket_hz(step_c, step0, phi, eps_k, ls_max_fn, theta, rho,
                       xtol = xtol, verbose = verbose)
  bracket <- br_res$bracket
  nfn <- nfn + br_res$nfn
  if (!br_res$ok) {
    LOpos <- which.min(bracket_props(bracket, 'f'))
    if (verbose) {
      message("Failed to create bracket, aborting, alpha = ",
              formatC(bracket[[LOpos]]$alpha))
    }
    return(list(step = bracket[[LOpos]], nfn = nfn))
  }

  # Check for T1/T2
  if (always_check_convergence) {
    LOpos <- hz_ok_bracket_pos(bracket, step0, c1, c2, eps_k,
                               strong_curvature = strong_curvature,
                               approx_armijo = approx_armijo)
    if (LOpos > 0) {
      if (verbose) {
        message("Bracket OK after bracket phase alpha = ",
                formatC(bracket[[LOpos]]$alpha))
      }
      return(list(step = bracket[[LOpos]], nfn = nfn))
    }
  }
  ls_max_fn <- max_fn - nfn
  if (ls_max_fn <= 0) {
    if (verbose) {
      message("max fn reached after bracket phase")
    }
    LOpos <- which.min(bracket_props(bracket, 'f'))
    return(list(step = bracket[[LOpos]], nfn = nfn))
  }

  old_bracket <- bracket
  # Zoom step
  while (TRUE) {
    if (verbose) {
      message("BRACKET: ", format_bracket(old_bracket))
    }
    # L1

    # do some of S1 here in case it's already an ok step
    alpha_c <- secant_hz(old_bracket[[1]], old_bracket[[2]])

    if (!is_finite_numeric(alpha_c)) {
      # probably only an issue when tolerances are very low and approx armijo
      # is off: can get NaN when bracket size approaches zero
      LOpos <- which.min(bracket_props(bracket, 'f'))
      if (verbose) {
        message("bad secant alpha, aborting line search")
      }
      return(list(step = bracket[[LOpos]], nfn = nfn))
    }
    fsec_res <- find_finite(phi, alpha_c, ls_max_fn,
                            min_alpha = bracket_min_alpha(old_bracket))
    nfn <- nfn + fsec_res$nfn
    ls_max_fn <- max_fn - nfn
    if (!fsec_res$ok) {
      if (verbose) {
        message("No finite alpha during secant bisection, aborting line search")
      }
      break
    }
    step_c <- fsec_res$step

    if (verbose) {
      message("S1: secant step_c alpha = ", formatC(alpha_c))
    }
    if (hz_ok_step(step_c, step0, c1, c2, eps_k,
                   strong_curvature = strong_curvature,
                   approx_armijo = approx_armijo)) {
      if (verbose) {
        message("step OK after secant alpha = ", formatC(step_c$alpha))
      }
      return(list(step = step_c, nfn = nfn))
    }
    ls_max_fn <- max_fn - nfn
    if (ls_max_fn <= 0) {
      if (verbose) {
        message("max fn reached after secant")
      }
      bracket[[3]] <- step_c
      LOpos <- which.min(bracket_props(bracket, 'f'))
      return(list(step = bracket[[LOpos]], nfn = nfn))
    }

    sec2_res <- secant2_hz(old_bracket, step_c, step0, phi, eps_k, ls_max_fn,
                           theta, verbose = verbose)
    bracket <- sec2_res$bracket
    nfn <- nfn + sec2_res$nfn
    if (!sec2_res$ok) {
      break
    }

    if (verbose) {
      message("new bracket: ", format_bracket(bracket))
    }
    # Check for T1/T2
    LOpos <- hz_ok_bracket_pos(bracket, step0, c1, c2, eps_k,
                               strong_curvature = strong_curvature,
                               approx_armijo = approx_armijo)
    if (LOpos > 0) {
      # if we only allow interpolated steps to count as converged,
      # check that any satisfactory result of secant2 came from a secant step
      # (i.e. wasn't already part of the bracket)
      if (always_check_convergence ||
          (bracket[[LOpos]]$alpha != old_bracket[[1]]$alpha &&
           bracket[[LOpos]]$alpha != old_bracket[[2]]$alpha)) {
        if (verbose) {
          message("Bracket OK after secant2 alpha = ", formatC(bracket[[LOpos]]$alpha))
        }
        return(list(step = bracket[[LOpos]], nfn = nfn))
      }
    }
    ls_max_fn <- max_fn - nfn
    if (ls_max_fn <= 0) {
      if (verbose) {
        message("max fn reached after secant2")
      }
      break
    }

    # L2
    # Ensure the bracket size decreased by a factor of gamma
    old_range <- bracket_props(old_bracket, 'alpha')
    old_diff <- abs(old_range[2] - old_range[1])
    new_range <- bracket_props(bracket, 'alpha')
    new_diff <- abs(new_range[2] - new_range[1])

    if (verbose) {
      message("Bracket size = ", formatC(new_diff), " old bracket size = ",
              formatC(old_diff))
    }

    if (new_diff < xtol * max(new_range)) {
      if (verbose) {
        message("Bracket size reduced below tolerance")
      }

      break
    }

    if (new_diff > old_diff * gamma) {
      if (verbose) {
        message("Bracket size did not decrease sufficiently: bisecting")
      }
      # bracket size wasn't decreased sufficiently: bisection step
      alpha_c <- mean(new_range)
      bisec_res <- find_finite(phi, alpha_c, ls_max_fn,
                               min_alpha = bracket_min_alpha(bracket))
      nfn <- nfn + bisec_res$nfn
      ls_max_fn <- max_fn - nfn
      if (!bisec_res$ok) {
        if (verbose) {
          message("No finite alpha during bisection, aborting line search")
        }
        break
      }
      step_c <- bisec_res$step

      if (ls_max_fn <= 0) {
        if (verbose) {
          message("Reached max_fn, returning")
        }

        bracket[[3]] <- step_c
        LOpos <- which.min(bracket_props(bracket, 'f'))
        return(list(step = bracket[[LOpos]], nfn = nfn))
      }

      up_res <- update_bracket_hz(bracket, step_c, step0, phi, eps_k, ls_max_fn,
                                  xtol = xtol, theta)
      bracket <- up_res$bracket
      nfn <- nfn + up_res$nfn
      if (!up_res$ok) {
        break
      }
      # Check for T1/T2
      if (always_check_convergence) {
        LOpos <- hz_ok_bracket_pos(bracket, step0, c1, c2, eps_k,
                                   strong_curvature = strong_curvature,
                                   approx_armijo = approx_armijo)
        if (LOpos > 0) {
          if (verbose) {
            message("Bracket OK after bisection alpha = ", formatC(bracket[[LOpos]]$alpha))
          }
          return(list(step = bracket[[LOpos]], nfn = nfn))
        }
      }
      ls_max_fn <- max_fn - nfn
      if (ls_max_fn <= 0) {
        if (verbose) {
          message("max fn reached after bisection")
        }
        break
      }
    }

    # L3
    old_bracket <- bracket
  }

  LOpos <- which.min(bracket_props(bracket, 'f'))
  list(step = bracket[[LOpos]], nfn = nfn)
}

# Bracket -----------------------------------------------------------------

# Routine 'bracket' B1-3
# Generates an initial bracket satisfying the opposite slope condition
# or if max_fn is reached or a non-finite f/g value generated, returns the best
# two values it can find.
bracket_hz <- function(step_c, step0, phi, eps, max_fn, theta = 0.5, rho = 5,
                       xtol = .Machine$double.eps, verbose = FALSE) {
  step_c_old <- step0
  # used only if bracket step fails (hit max_fn or non-finite f/g)
  step_c_old_old <- step0

  ls_max_fn <- max_fn
  nfn <- 0

  # step_c is the latest attempt at a bracketed point
  # step_c_old is the previous step_c
  while (TRUE) {
    if (verbose) {
      message("Bracketing: step = ", formatC(step_c$alpha))
    }
    if (step_c$d > 0) {
      # B1 slope is +ve: bracketing successful: return [step_c_old, step_c]
      if (verbose) {
        message("B1: slope +ve")
      }
      return(list(bracket = list(step_c_old, step_c), nfn = nfn, ok = TRUE))
    }
    if (step_c$f > step0$f + eps) {
      # B2 slope is -ve but f is higher than starting point
      # we must have stepped too far beyond the minimum and the +ve slope
      # and its maximum
      # find the minimum by weighted bisection
      if (verbose) {
        message("B2: f > phi0 + eps")
      }

      if (ls_max_fn <= 0) {
        # avoid returning step_c in this case, which might be outside the
        # current minimizer "basin"
        return(list(bracket = list(step_c_old_old, step_c_old), nfn = nfn,
                    ok = FALSE))
      }

      # Probably could use step_c_old as LHS of bracket
      # but HZ paper specifies step0
      bracket_sub = list(step0, step_c)
      bisect_res <- update_bracket_bisect_hz(bracket_sub, step0, phi, eps,
                                             ls_max_fn, theta,
                                             xtol = xtol,
                                             verbose = verbose)
      bisect_res$nfn <- bisect_res$nfn + nfn
      # return bisection result: may have failed
      return(bisect_res)
    }

    # B3: slope is -ve and f < f0, so we haven't passed the minimum yet
    if (ls_max_fn <= 0) {
      return(list(bracket = list(step_c_old, step_c), nfn = nfn, ok = FALSE))
    }

    # extrapolate: increase the step size
    step_c_old_old <- step_c_old
    step_c_old <- step_c
    alpha <- step_c$alpha * rho

    fext_res <- find_finite(phi, alpha, ls_max_fn, min_alpha = step_c$alpha)
    nfn <- nfn + fext_res$nfn
    ls_max_fn <- max_fn - nfn
    if (!fext_res$ok) {
      if (verbose) {
        message("No finite alpha during extrapolation bisection, aborting line search")
      }
      return(list(bracket = list(step_c_old_old, step_c_old), nfn = nfn,
                  ok = FALSE))
    }
    step_c <- fext_res$step
  }
}

# Update ------------------------------------------------------------------

# routine 'update' U0-U3
# Given a bracket, create a new bracket with end points which are inside
# the original bracket.
# Returns with ok = TRUE if any of U0-U3 succeed (i.e. U0 is not a failure).
# Returns with ok = FALSE if U3 fails (i.e. exceed max_fn or non-finite
# function/gradient is calculated before bisection succeeds)
update_bracket_hz <- function(bracket, step_c, step0, phi, eps, max_fn,
                              theta = 0.5, xtol = .Machine$double.eps,
                              verbose = FALSE) {
  if (verbose) {
    message("U: alpha = ", formatC(step_c$alpha),
            " bracket = ", format_bracket(bracket))
  }
  nfn <- 0
  ok <- TRUE

  if (!is_in_bracket(bracket, step_c$alpha)) {
    # U0: c is not inside, reject it
    new_bracket <- bracket
    if (verbose) {
      message("U0: step not in bracket, reject")
    }
  }
  else if (step_c$d >= 0) {
    # U1: c is on the +ve slope, make it the new hi
    new_bracket <- list(bracket[[1]], step_c)
    if (verbose) {
      message("U1: step has +ve slope, new hi")
    }
  }
  else if (step_c$f <= step0$f + eps) {
    # U2: c is on the -ve slope and closer to minimum than a
    # make it the new lo
    new_bracket <- list(step_c, bracket[[2]])
    if (verbose) {
      message("U2: step has -ve slope and closer to minimum, new lo")
    }
  }
  else {
    # U3
    # c is on the -ve slope but larger than f0: must have missed the minimum
    # and the +ve slope and the maximum
    # find new hi by weighted bisection
    if (verbose) {
      message("U3: step has -ve slope but not closer to minimum, bisect")
    }
    sub_bracket <- list(bracket[[1]], step_c)
    sub_res <- update_bracket_bisect_hz(sub_bracket, step0, phi, eps, max_fn,
                                        theta, xtol = xtol,
                                        verbose = verbose)
    new_bracket <- sub_res$bracket
    nfn <- sub_res$nfn
    ok <- sub_res$ok
  }

  list(bracket = new_bracket, nfn = nfn, ok = ok)
}

# U3a-c from routine 'update'
#
# Use weighted bisection of the current bracket so that bracket[2] contains a
# step size with a +ve slope. bracket[1] will also be updated if a point
# with -ve slope closer to the minimizer is found.
#
# Also used during the bracket step if the step size gets too large.
# Called when step size leads to a -ve slope but f is > f0, implying that
# step size was so large it missed the local minimum, the +ve slope and the
# local maximum and we are now going downhill to some other minimum.
# Use weighted bisection until the hi of the bracket has a +ve slope
# lo of bracket will also be updated if we find a suitable point during
# bisection.
#
# If bisection succeeds, then this function returns with ok = TRUE.
# If the number of bisections exceeds max_fn, or if any step size contains a
# non-finite slope or function value, the most recent finite-valued bracket is
# returned with ok = FALSE
update_bracket_bisect_hz <- function(bracket, step0, phi, eps, max_fn,
                                     theta = 0.5,
                                     xtol = .Machine$double.eps,
                                     verbose = FALSE) {
  res <- bracket
  nfn <- 0
  ls_max_fn <- max_fn
  ok <- FALSE
  while (TRUE) {
    if (verbose) {
      message("U3: Bracket: ", format_bracket(res), " width = ",
              bracket_width(res))
    }
    if (bracket_width(res) <= xtol * res[[2]]$alpha) {
      if (verbose) {
        message("Relative bracket width reduced below tolerance, aborting")
      }
      break
    }

    ls_max_fn <- max_fn - nfn
    if (ls_max_fn <= 0) {
      if (verbose) {
        message("max_fn reached, aborting bisection bracket update")
      }
      break
    }

    # U3a new point is (weighted) bisection of current bracket
    alpha <- (1 - theta) * res[[1]]$alpha + theta * res[[2]]$alpha
    fwbi_res <- find_finite(phi, alpha, ls_max_fn,
                            min_alpha = bracket_min_alpha(res))
    nfn <- nfn + fwbi_res$nfn
    ls_max_fn <- max_fn - nfn
    if (!fwbi_res$ok) {
      if (verbose) {
        message("No finite alpha during weighted bisection, aborting line search")
      }
      break
    }
    step_d <- fwbi_res$step

    if (step_d$d >= 0) {
      # d is on +ve slope, make it the new hi and return
      res[[2]] <- step_d
      ok <- TRUE
      break
    }
    if (step_d$f <= step0$f + eps) {
      if (verbose) {
        message("U3b: alpha ", formatC(step_d$alpha),
                " f = ", formatC(step_d$f),
                " d  = ", formatC(step_d$d),
                " closer to minimizer: new lo")
      }
      # U3b: d is on -ve slope but closer to minimizer, make it new lo and loop
      res[[1]] <- step_d
    } else {
      # U3c: d has -ve slope but still > f0 so still too large a step,
      # make it the new hi and loop
      if (verbose) {
        message("U3b: alpha ", formatC(step_d$alpha),
                " f = ", formatC(step_d$f),
                " d  = ", formatC(step_d$d),
                " -ve slope but > f0: new hi")
      }
      res[[2]] <- step_d
    }
  }

  list(bracket = res, nfn = nfn, ok = ok)
}

# Secant ------------------------------------------------------------------

# Routine 'secant2'
# Do the secant step to generate c for step S1 outside of this routine because
# it may be an acceptable step without having to update any brackets
secant2_hz <- function(bracket, step_c, step0, phi, eps, max_fn,
                       theta = 0.5, xtol = .Machine$double.eps,
                       verbose = FALSE) {
  nfn <- 0
  ls_max_fn <- max_fn

  if (ls_max_fn <= 0) {
    return(list(bracket = bracket, nfn = nfn, ok = FALSE))
  }

  if (verbose) {
    message("S1: Creating AB")
  }
  bracket_AB_res <- update_bracket_hz(bracket, step_c, step0, phi, eps, theta,
                                      xtol = xtol, verbose = verbose)
  if (!bracket_AB_res$ok) {
    return(list(bracket = bracket, nfn = nfn, ok = FALSE))
  }
  bracket_AB <- bracket_AB_res$bracket

  if (verbose) {
    message("S1: secant alpha = ", formatC(step_c$alpha),
            " ab = ", format_bracket(bracket),
            " AB = ", format_bracket(bracket_AB))
  }

  ok <- TRUE
  # following two if blocks rely on exact floating point comparison
  if (step_c$alpha == bracket_AB[[2]]$alpha) {
    # S2 c == B
    alpha_cbar <- secant_hz(bracket[[2]], bracket_AB[[2]])

    if (verbose) {
      message("S2 c = B: cbar = secant(b, B) = (",
              formatC(bracket[[2]]$alpha), ", ",
              formatC(bracket_AB[[2]]$alpha), ") = ",
              formatC(alpha_cbar))
    }
    # update routine would also check that c_bar is in [A,B] but do it manually
    # here to avoid calculating phi(c_bar) if we don't need to
    if (is.finite(alpha_cbar) && is_in_bracket(bracket_AB, alpha_cbar)
        && max_fn > 0) {
      step_cbar <- phi(alpha_cbar)
      nfn <- nfn + 1
      max_fn <- ls_max_fn - nfn

      res <- update_bracket_hz(bracket_AB, step_cbar, step0, phi, eps, max_fn,
                               theta, xtol = xtol, verbose = verbose)
      new_bracket <- res$bracket
      nfn <- nfn + res$nfn
      max_fn <- ls_max_fn - nfn
      ok <- res$ok
    }
    else {
      new_bracket <- bracket_AB
    }
    if (verbose) {
      message("S2 bracket: ", format_bracket(new_bracket))
    }
  }
  else if (step_c$alpha == bracket_AB[[1]]$alpha) {
    # S3 c == A
    alpha_cbar <- secant_hz(bracket[[1]], bracket_AB[[1]])
    if (verbose) {
      message("S3 c = A: cbar = secant(a, A) = (",
              formatC(bracket[[1]]$alpha), ", ",
              formatC(bracket_AB[[1]]$alpha), ") = ",
              formatC(alpha_cbar))
    }
    if (is.finite(alpha_cbar) && is_in_bracket(bracket_AB, alpha_cbar)
        && max_fn > 0) {

      step_cbar <- phi(alpha_cbar)
      nfn <- nfn + 1
      max_fn <- ls_max_fn - nfn

      res <- update_bracket_hz(bracket_AB, step_cbar, step0, phi, eps, max_fn,
                               theta, xtol = xtol, verbose = verbose)
      new_bracket <- res$bracket
      nfn <- nfn + res$nfn
      max_fn <- ls_max_fn - nfn
      ok <- res$ok
    }
    else {
      new_bracket <- bracket_AB
    }
    if (verbose) {
      message("S3 bracket: ", format_bracket(new_bracket))
    }
  }
  else {
    # S4
    new_bracket <- bracket_AB
    if (verbose) {
      message("S4 bracket: ", format_bracket(new_bracket))
    }
  }

  list(bracket = new_bracket, nfn = nfn, ok = ok)
}

secant_hz <- function(step_a, step_b) {
  (step_a$alpha * step_b$d - step_b$alpha * step_a$d) / (step_b$d - step_a$d)
}


# Termination Conditions --------------------------------------------------

hz_ok_step <- function(step, step0, c1, c2, eps, strong_curvature = FALSE,
                       approx_armijo = TRUE) {
  if (strong_curvature) {
    ok <- strong_curvature_ok_step(step0, step, c2)
  }
  else {
    ok <- curvature_ok_step(step0, step, c2)
  }
  if (!ok) {
    return(ok)
  }

  if (armijo_ok_step(step0, step, c1)) {
    return(ok)
  }

  approx_armijo && (step$f <= step0$f + eps) &&
    approx_armijo_ok_step(step0, step, c1)
}

hz_ok_bracket_pos <- function(bracket, step0, c1, c2, eps,
                          strong_curvature = FALSE,
                          approx_armijo = TRUE) {
  ok_pos <- 0
  if (hz_ok_step(bracket[[1]], step0, c1, c2, eps,
                 strong_curvature = strong_curvature,
                 approx_armijo = approx_armijo)) {
    ok_pos <- 1
  }
  if (hz_ok_step(bracket[[2]], step0, c1, c2, eps,
                 strong_curvature = strong_curvature,
                 approx_armijo = approx_armijo)) {
    # if somehow we've reached a situation where both sides of the bracket
    # meet the conditions, choose the one with the lower function value
    if (ok_pos == 0 || bracket[[2]]$f < bracket[[1]]$f) {
      ok_pos <- 2
    }
  }
  ok_pos
}

