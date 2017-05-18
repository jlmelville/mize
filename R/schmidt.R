# Translation of Mark Schmidt's minFunc line search code for satisfying the
# Strong Wolfe conditions (and also the Armijo conditions)
# http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.

# Adapters ----------------------------------

# Uses the default line search settings: cubic interpolation/extrapolation
# Falling back to Armijo backtracking (also using cubic interpolation) if
# a non-legal value is found
schmidt <- function(c1 = c2 / 2, c2 = 0.1, max_fn = Inf, eps = 1e-6,
                    strong_curvature = TRUE, approx_armijo = FALSE) {
  if (c2 < c1) {
    stop("schmidt line search: c2 < c1")
  }
  function(phi, step0, alpha,
           total_max_fn = Inf, total_max_gr = Inf, total_max_fg = Inf,
           pm) {
    maxfev <- min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    if (maxfev <= 0) {
      return(list(step = step0, nfn = 0, ngr = 0))
    }

    if (approx_armijo) {
      armijo_check_fn <- make_approx_armijo_ok_step(eps)
    }
    else {
      armijo_check_fn <- armijo_ok_step
    }

    if (strong_curvature) {
      curvature_check_fn <- strong_curvature_ok_step
    }
    else {
      curvature_check_fn <- curvature_ok_step
    }

    res <- WolfeLineSearch(alpha = alpha, f = step0$f, g = step0$df,
                           gtd = step0$d,
                           c1 = c1, c2 = c2, LS_interp = 2, LS_multi = 0,
                           maxLS = maxfev,
                           funObj = phi, varargin = NULL,
                           pnorm_inf = max(abs(pm)),
                           progTol = 1e-9,
                           armijo_check_fn = armijo_check_fn,
                           curvature_check_fn = curvature_check_fn,
                           debug = FALSE)
    res$ngr = res$nfn
    res
  }
}

# step_down if non-NULL, multiply the step size by this value when backtracking
# Otherwise, use a cubic interpolation based on previous function and derivative
# values
schmidt_armijo_backtrack <- function(c1 = 0.05, step_down = NULL, max_fn = Inf) {
  function(phi, step0, alpha,
           total_max_fn = Inf, total_max_gr = Inf, total_max_fg = Inf,
           pm) {
    maxfev <- min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    if (maxfev <= 0) {
      return(list(step = step0, nfn = 0, ngr = 0))
    }
    if (!is.null(step_down)) {
      # fixed-size step reduction by a factor of step_down
      LS_interp <- 0
    }
    else {
      # cubic interpolation
      LS_interp <- 2
    }

    res <- ArmijoBacktrack(step = alpha, f = step0$f, g = step0$df,
                           gtd = step0$d,
                           c1 = c1, LS_interp = LS_interp, LS_multi = 0,
                           maxLS = maxfev, step_down = step_down,
                           funObj = phi, varargin = NULL,
                           pnorm_inf = max(abs(pm)),
                           progTol = 1e-9, debug = FALSE)

    res$ngr <- res$nfn
    res
  }
}

# Translated minFunc routines ---------------------------------------------

# Bracketing Line Search to Satisfy Wolfe Conditions
#
# Inputs:
#   x: starting location
#   step: initial step size
#   d: descent direction
#   f: function value at starting location
#   g: gradient at starting location
#   gtd: directional derivative at starting location
#   c1: sufficient decrease parameter
#   c2: curvature parameter
#   debug: display debugging information
#   LS_interp: type of interpolation
#   maxLS: maximum number of funEvals (changed from matlab original)
#   progTol: minimum allowable step length
#   funObj: objective function
#   varargin: parameters of objective function
#
# For the Wolfe line-search, these interpolation strategies are available ('LS_interp'):
#   - 0 : Step Size Doubling and Bisection
#   - 1 : Cubic interpolation/extrapolation using new function and gradient values (default)
#   - 2 : Mixed quadratic/cubic interpolation/extrapolation
# Outputs:
#   step: step length
#   f_new: function value at x+step*d
#   g_new: gradient value at x+step*d
#   funEvals: number function evaluations performed by line search
#   H: Hessian at initial guess (only computed if requested
# @returns [step,f_new,g_new,funEvals,H]
WolfeLineSearch <-
  function(alpha, f, g, gtd,
           c1 = 1e-4, c2 = 0.1, LS_interp = 1, LS_multi = 0, maxLS = 25,
           funObj, varargin = NULL,
           pnorm_inf, progTol = 1e-9, armijo_check_fn = armijo_ok_step,
           curvature_check_fn = strong_curvature_ok_step,
           debug = FALSE) {
    # Bracket an Interval containing a point satisfying the
    # Wolfe criteria
    step0 <- list(alpha = 0, f = f, df = g, d = gtd)

    bracket_res <- schmidt_bracket(alpha, LS_interp, maxLS, funObj, step0,
                                   c1, c2, armijo_check_fn, curvature_check_fn,
                                   debug)
    bracket_step <- bracket_res$bracket
    funEvals <- bracket_res$funEvals
    done <- bracket_res$done

    if (!bracket_is_finite(bracket_step)) {
      if (debug) {
        message('Switching to Armijo line-search')
      }
      alpha <- mean(bracket_props(bracket_step, 'alpha'))

      # Do Armijo
      armijo_res <- ArmijoBacktrack(alpha, step0$f, step0$df, step0$d,
                                    c1 = c1, LS_interp = LS_interp,
                                    LS_multi = LS_multi,
                                    maxLS = maxLS - funEvals,
                                    funObj = funObj, varargin = NULL,
                                    pnorm_inf = pnorm_inf,
                                    progTol = progTol, debug = debug)

      armijo_res$nfn <- armijo_res$nfn + funEvals
      return(armijo_res)
    }

    ## Zoom Phase
    # We now either have a point satisfying the criteria, or a bracket
    # surrounding a point satisfying the criteria
    # Refine the bracket until we find a point satisfying the criteria
    if (!done) {
      maxLS <- maxLS - funEvals
      zoom_res <- schmidt_zoom(bracket_step, LS_interp, maxLS, funObj,
                               step0, c1, c2, pnorm_inf, progTol,
                               armijo_check_fn, curvature_check_fn,
                               debug)

      funEvals <- funEvals + zoom_res$funEvals
      bracket_step <- zoom_res$bracket
    }

    list(step = best_bracket_step(bracket_step), nfn = funEvals)
  }

# Change from original: maxLS refers to maximum allowed funEvals, not LS iters
schmidt_bracket <- function(alpha, LS_interp, maxLS, funObj, step0, c1, c2,
                            armijo_check_fn, curvature_check_fn, debug) {
  # did we find a bracket
  ok <- FALSE
  # did we find a step that already fulfils the line search
  done <- FALSE

  step_prev <- step0

  # Evaluate the Objective and Gradient at the Initial Step
  ff_res <- find_finite(funObj, alpha, maxLS, min_alpha = 0)
  funEvals <- ff_res$nfn
  step_new <- ff_res$step

  LSiter <- 0
  while (funEvals < maxLS) {

    if (!ff_res$ok) {
      if (debug) {
        message('Extrapolated into illegal region, returning')
      }
      bracket_step <- list(step_prev, step_new)

      return(list(bracket = bracket_step, done = done,
                  funEvals = funEvals, ok = FALSE))
    }

    # See if we have found the other side of the bracket
    if (!armijo_check_fn(step0, step_new, c1) ||
        (LSiter > 1 && step_new$f >= step_prev$f)) {
      bracket_step <- list(step_prev, step_new)
      ok <- TRUE
      if (debug) {
        message('Armijo failed or step_new$f >= step_prev$f: bracket is [prev new]')
      }
      break
    }
    else if (curvature_check_fn(step0, step_new, c2)) {
      bracket_step <- list(step_new)
      ok <- TRUE
      done <- TRUE

      if (debug) {
        message('Sufficient curvature found: bracket is [new]')
      }
      break
    }
    else if (step_new$d >= 0) {
      bracket_step <- list(step_prev, step_new)

      if (debug) {
        message('step_new$d >= 0: bracket is [prev, new]')
      }
      break
    }

    minStep <- step_new$alpha + 0.01 * (step_new$alpha - step_prev$alpha)
    maxStep <- step_new$alpha * 10

    if (LS_interp <= 1) {
      if (debug) {
        message('Extending Bracket')
      }
      alpha_new <- maxStep
    }
    else if (LS_interp == 2) {
      if (debug) {
        message('Cubic Extrapolation')
      }
      alpha_new <- polyinterp(point_matrix_step(step_prev, step_new),
                              minStep, maxStep)
    }
    else {
      # LS_interp == 3
      alpha_new <- mixedExtrap_step(step_prev, step_new, minStep, maxStep,
                                    debug)
    }

    step_prev <- step_new

    ff_res <- find_finite(funObj, alpha_new, maxLS - funEvals,
                          min_alpha = step_prev$alpha)
    funEvals <- funEvals + ff_res$nfn
    step_new <- ff_res$step

    LSiter <- LSiter + 1
  }

  # If we ran out of fun_evals, need to repeat finite check for last iteration
  if (!ok && !ff_res$ok) {
    if (debug) {
      message('Extrapolated into illegal region, returning')
    }
  }

  if (funEvals >= maxLS && !ok) {
    if (debug) {
      message("max_fn reached in bracket step")
    }
  }

  list(bracket = bracket_step, done = done, funEvals = funEvals, ok = ok)
}

# Change from original: maxLS refers to max allowed funEvals not LSiters
schmidt_zoom <- function(bracket_step, LS_interp, maxLS, funObj,
                         step0, c1, c2, pnorm_inf, progTol, armijo_check_fn,
                         curvature_check_fn,
                         debug) {
  insufProgress <- FALSE
  Tpos <- 2 # position in the bracket of the current best step
  # mixed interp only: if true, save point from previous bracket
  LOposRemoved <- FALSE

  funEvals <- 0

  done <- FALSE
  while (!done && funEvals < maxLS) {
    # Find High and Low Points in bracket
    LOpos <- which.min(bracket_props(bracket_step, 'f'))
    HIpos <- -LOpos + 3 # 1 or 2, whichever wasn't the LOpos

    # Compute new trial value
    if (LS_interp <= 1 || !bracket_is_finite(bracket_step)) {
      if (!bracket_is_finite(bracket_step)) {
        message("Bad f/g in bracket - bisecting")
      }
      alpha <- mean(bracket_props(bracket_step, 'alpha'))
      if (debug) {
        message('Bisecting: trial step = ', formatC(alpha))
      }
    }
    else if (LS_interp == 2) {
      alpha <- polyinterp(
        point_matrix_step(bracket_step[[1]], bracket_step[[2]]),
        debug = debug)
      if (debug) {
        message('Grad-Cubic Interpolation: trial step = ', formatC(alpha))
      }
    }
    else {
      # Mixed Case #
      nonTpos <- -Tpos + 3
      if (!LOposRemoved) {
        oldLO <- bracket_step[[nonTpos]]
      }
      alpha <- mixedInterp_step(bracket_step, Tpos, oldLO, debug)
      if (debug) {
        message('Mixed Interpolation: trial step = ', formatC(alpha))
      }
    }

    # Ensure that alpha is finite
    if (!is.finite(alpha)) {
      alpha <- mean(bracket_props(bracket_step, 'alpha'))
      if (debug) {
        message("Non-finite trial alpha, bisecting: alpha = ", formatC(alpha))
      }
    }

    # Test that we are making sufficient progress
    bracket_alphas <- bracket_props(bracket_step, 'alpha')
    alpha_max <- max(bracket_alphas)
    alpha_min <- min(bracket_alphas)
    alpha_range <- alpha_max - alpha_min
    if (alpha_range > 0) {
      if (min(alpha_max - alpha, alpha - alpha_min) / alpha_range < 0.1) {
        if (debug) {
          message('Interpolation close to boundary')
        }

        if (insufProgress || alpha >= alpha_max || alpha <= alpha_min) {
          if (debug) {
            message('Evaluating at 0.1 away from boundary')
          }
          if (abs(alpha - alpha_max) < abs(alpha - alpha_min)) {
            alpha <- alpha_max - 0.1 * alpha_range
          }
          else {
            alpha <- alpha_min + 0.1 * alpha_range
          }
          insufProgress <- FALSE
        }
        else {
          insufProgress <- TRUE
        }
      }
      else {
        insufProgress <- FALSE
      }
    }

    # Evaluate new point

    # code attempts to handle non-finite values but this is easier in Matlab
    # where NaN can safely be compared with finite values (returning 0 in all
    # comparisons), whereas R returns NA. Instead, let's attempt to find a
    # finite value by bisecting repeatedly. If we run out of evaluations or
    # hit the bracket, we give up.
    ff_res <- find_finite(funObj, alpha, maxLS - funEvals,
                          min_alpha = bracket_min_alpha(bracket_step))
    funEvals <- funEvals + ff_res$nfn
    if (!ff_res$ok) {
      if (debug) {
        message("Failed to find finite legal step size in zoom phase, aborting")
      }
      break
    }
    step_new <- ff_res$step

    # Update bracket
    if (!armijo_check_fn(step0, step_new, c1) ||
        step_new$f >= bracket_step[[LOpos]]$f) {
      if (debug) {
        message("New point becomes new HI")
      }
      # Armijo condition not satisfied or not lower than lowest point
      bracket_step[[HIpos]] <- step_new
      Tpos <- HIpos
      # [LO, new]
    }
    else {
      if (curvature_check_fn(step0, step_new, c2)) {
        # Wolfe conditions satisfied
        done <- TRUE
        # [new, HI]
      }
      else if (step_new$d * (bracket_step[[HIpos]]$alpha - bracket_step[[LOpos]]$alpha) >= 0) {
        if (debug) {
          message("Old LO becomes new HI")
        }
        # Old HI becomes new LO
        bracket_step[[HIpos]] <- bracket_step[[LOpos]]

        if (LS_interp == 3) {
          if (debug) {
            message('LO Pos is being removed!')
          }
          LOposRemoved <- TRUE
          oldLO <- bracket_step[[LOpos]]
        }
        # [new, LO]
      }
      # else [new, HI]

      if (debug) {
        message("New point becomes new LO")
      }
      # New point becomes new LO
      bracket_step[[LOpos]] <- step_new
      Tpos <- LOpos
    }

    if (!done && bracket_width(bracket_step) * pnorm_inf < progTol) {
      if (debug) {
        message('Line-search bracket has been reduced below progTol')
      }
      break
    }
  } # end of while loop
  if (funEvals >= maxLS) {
    if (debug) {
      message('Line Search Exceeded Maximum Function Evaluations')
    }
  }
  list(bracket = bracket_step, funEvals = funEvals)
}

# Backtracking linesearch to satisfy Armijo condition
#
# Inputs:
#   x: starting location
#   t: initial step size
#   d: descent direction
#   f: function value at starting location
#   gtd: directional derivative at starting location
#   c1: sufficient decrease parameter
#   debug: display debugging information
#   LS_interp: type of interpolation
#   progTol: minimum allowable step length
#   doPlot: do a graphical display of interpolation
#   funObj: objective function
#   varargin: parameters of objective function
#
# For the Armijo line-search, several interpolation strategies are available
# ('LS_interp'):
#   - 0 : Step size halving
#   - 1 : Polynomial interpolation using new function values
#   - 2 : Polynomial interpolation using new function and gradient values (default)
#
# When (LS_interp = 1), the default setting of (LS_multi = 0) uses quadratic
# interpolation, while if (LS_multi = 1) it uses cubic interpolation if more
# than one point are available.
#
# When (LS_interp = 2), the default setting of (LS_multi = 0) uses cubic interpolation,
# while if (LS_multi = 1) it uses quartic or quintic interpolation if more than
# one point are available
#
# Outputs:
#   t: step length
#   f_new: function value at x+t*d
#   g_new: gradient value at x+t*d
#   funEvals: number function evaluations performed by line search
#
# recet change: LS changed to LS_interp and LS_multi

ArmijoBacktrack <-
  function(step, f, g, gtd,
           c1 = 1e-4,
           LS_interp = 2, LS_multi = 0, maxLS = Inf,
           step_down = 0.5,
           funObj,
           varargin = NULL,
           pnorm_inf,
           progTol = 1e-9, debug = FALSE)
  {
    # Evaluate the Objective and Gradient at the Initial Step

    f_prev <- NA
    t_prev <- NA
    g_prev <- NA
    gtd_prev <- NA

    fun_obj_res <- funObj(step)
    f_new <- fun_obj_res$f
    g_new <- fun_obj_res$df
    gtd_new <- fun_obj_res$d

    funEvals <- 1

    while (funEvals < maxLS && (f_new > f + c1 * step * gtd || !is.finite(f_new))) {
      temp <- step
      if (LS_interp == 0 || !is.finite(f_new)) {
        # Ignore value of new point
        if (debug) {
          message('Fixed BT')
        }
        step <- step_down * step
      }
      else if (LS_interp == 1 || !is.finite(g_new)) {
        # Use function value at new point, but not its derivative
        if (funEvals < 2 || LS_multi == 0 || !is.finite(f_prev)) {
          # Backtracking w/ quadratic interpolation based on two points
          if (debug) {
            message('Quad BT')
          }
          step <- polyinterp(point_matrix(c(0, step), c(f, f_new), c(gtd, NA)),
                             0, step)
        }
        else {
          # Backtracking w/ cubic interpolation based on three points
          if (debug) {
            message('Cubic BT')
          }
          step <-
            polyinterp(point_matrix(
              c(0, step, t_prev), c(f, f_new, f_prev), c(gtd, NA, NA)),
              0, step)
        }
      }
      else {
        # Use function value and derivative at new point
        if (funEvals < 2 || LS_multi == 0 || !is.finite(f_prev)) {
          # Backtracking w/ cubic interpolation w/ derivative
          if (debug) {
            message('Grad-Cubic BT')
          }
          step <- polyinterp(
            point_matrix(c(0, step), c(f, f_new), c(gtd, gtd_new)),
            0, step)
        }
        else if (!is.finite(g_prev)) {
          # Backtracking w/ quartic interpolation 3 points and derivative
          # of two
          if (debug) {
            message('Grad-Quartic BT')
          }

          step <- polyinterp(point_matrix(
            c(0, step, t_prev), c(f, f_new, f_prev), c(gtd, gtd_new, NA)),
            0, step)
        }
        else {
          # Backtracking w/ quintic interpolation of 3 points and derivative
          # of three
          if (debug) {
            message('Grad-Quintic BT')
          }

          step <- polyinterp(point_matrix(
            c(0, step, t_prev),
            c(f, f_new, f_prev),
            c(gtd, gtd_new, gtd_prev)),
            0, step)
        }
      }


      if (!is_finite_numeric(step)) {
        step <- temp * 0.6
      }
      # Adjust if change in step is too small/large
      if (step < temp * 1e-3) {
        if (debug) {
          message('Interpolated Value Too Small, Adjusting')
        }
        step <- temp * 1e-3

      } else if (step > temp * 0.6) {
        if (debug) {
          message('Interpolated Value Too Large, Adjusting')
        }
        step <- temp * 0.6
      }

      # Store old point if doing three-point interpolation
      if (LS_multi) {
        f_prev <- f_new
        t_prev <- temp

        if (LS_interp == 2) {
          g_prev <- g_new
          gtd_prev <- gtd_new
        }
      }

      fun_obj_res <- funObj(step)
      f_new <- fun_obj_res$f
      g_new <- fun_obj_res$df
      gtd_new <- fun_obj_res$d

      funEvals <- funEvals + 1

      # Check whether step size has become too small
      if (pnorm_inf * step <= progTol) {
        if (debug) {
          message('Backtracking Line Search Failed')
        }
        step <- 0
        f_new <- f
        g_new <- g
        gtd_new <- gtd
        break
      }
    }

    list(
      step = list(alpha = step, f = f_new, df = g_new, d = gtd_new),
      nfn = funEvals
    )
  }

mixedExtrap_step <- function(step0, step1, minStep, maxStep, debug) {
  mixedExtrap(step0$alpha, step0$f, step0$d, step1$alpha, step1$f, step1$d,
              minStep, maxStep, debug)
}

mixedExtrap <- function(x0, f0, g0, x1, f1, g1, minStep, maxStep, debug) {
  alpha_c <- polyinterp(point_matrix(c(x0, x1), c(f0, f1), c(g0, g1)),
                        minStep, maxStep, debug = debug)
  alpha_s <- polyinterp(point_matrix(c(x0, x1), c(f0, NA), c(g0, g1)),
                        minStep, maxStep, debug = debug)
  if (debug) {
    message("cubic ext = ", formatC(alpha_c), " secant ext = ", formatC(alpha_s),
            " minStep = ", formatC(minStep),
            " alpha_c > minStep ? ", alpha_c > minStep,
            " |ac - x1| = ", formatC(abs(alpha_c - x1)),
            " |as - x1| = ", formatC( abs(alpha_s - x1))
    )
  }
  if (alpha_c > minStep && abs(alpha_c - x1) < abs(alpha_s - x1)) {
    if (debug) {
      message('Cubic Extrapolation ', formatC(alpha_c))
    }
    res <- alpha_c
  }
  else {
    if (debug) {
      message('Secant Extrapolation ', formatC(alpha_s))
    }
    res <- alpha_s
  }
  res
}

mixedInterp_step <- function(bracket_step,
                        Tpos,
                        oldLO,
                        debug) {

  bracket <- c(bracket_step[[1]]$alpha, bracket_step[[2]]$alpha)
  bracketFval <- c(bracket_step[[1]]$f, bracket_step[[2]]$f)
  bracketDval <- c(bracket_step[[1]]$d, bracket_step[[2]]$d)

  mixedInterp(bracket, bracketFval, bracketDval, Tpos,
              oldLO$alpha, oldLO$f, oldLO$d, debug)
}

mixedInterp <- function(
  bracket, bracketFval, bracketDval,
  Tpos,
  oldLOval, oldLOFval, oldLODval,
  debug) {

  # Mixed Case
  nonTpos <- -Tpos + 3


    gtdT <- bracketDval[Tpos]
    gtdNonT <- bracketDval[nonTpos]
    oldLOgtd <- oldLODval
    if (bracketFval[Tpos] > oldLOFval) {
      alpha_c <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)))
      alpha_q <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, NA)))
      if (abs(alpha_c - oldLOval) < abs(alpha_q - oldLOval)) {
        if (debug) {
          message('Cubic Interpolation')
        }
        res <- alpha_c
      }
      else {
        if (debug) {
          message('Mixed Quad/Cubic Interpolation')
        }
        res <- (alpha_q + alpha_c) / 2
      }
    }
    else if (dot(gtdT, oldLOgtd) < 0) {
      alpha_c <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)))
      alpha_s <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, NA),
        c(oldLOgtd, gtdT)))
      if (abs(alpha_c - bracket[Tpos]) >= abs(alpha_s - bracket[Tpos])) {
        if (debug) {
          message('Cubic Interpolation')
        }
        res <- alpha_c
      }
      else {
        if (debug) {
          message('Quad Interpolation')
        }
        res <- alpha_s
      }
    }
    else if (abs(gtdT) <= abs(oldLOgtd)) {
      alpha_c <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)), min(bracket), max(bracket))
      alpha_s <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(NA, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)), min(bracket), max(bracket))

      if (alpha_c > min(bracket) && alpha_c < max(bracket)) {
        if (abs(alpha_c - bracket[Tpos]) < abs(alpha_s - bracket[Tpos])) {
          if (debug) {
            message('Bounded Cubic Extrapolation')
          }
          res <- alpha_c
        }
        else {
          if (debug) {
            message('Bounded Secant Extrapolation')
          }
          res <- alpha_s
        }
      }
      else {
        if (debug) {
          message('Bounded Secant Extrapolation')
        }
        res <- alpha_s
      }

      if (bracket[Tpos] > oldLOval) {
        res <- min(bracket[Tpos] + 0.66 * (bracket[nonTpos] - bracket[Tpos]),
                   res)
      }
      else {
        res <- max(bracket[Tpos] + 0.66 * (bracket[nonTpos] - bracket[Tpos]),
                   res)
      }
    }
    else {
      res <- polyinterp(point_matrix(
        c(bracket[nonTpos], bracket[Tpos]),
        c(bracketFval[nonTpos], bracketFval[Tpos]),
        c(gtdNonT, gtdT)))
    }
    res
  }

# function [minPos] <- polyinterp(points,doPlot,xminBound,xmaxBound)
#
#   Minimum of interpolating polynomial based on function and derivative
#   values
#
#   It can also be used for extrapolation if {xmin,xmax} are outside
#   the domain of the points.
#
#   Input:
#       points(pointNum,[x f g])
#       xmin: min value that brackets minimum (default: min of points)
#       xmax: max value that brackets maximum (default: max of points)
#
#   set f or g to sqrt(-1) if they are not known
#   the order of the polynomial is the number of known f and g values minus 1
# points position, function and gradient values to interpolate.
# An n x 3 matrix where n is the number of points and each row contains
# x, f, g in columns 1-3 respectively.
# @return minPos
polyinterp <- function(points,
                       xminBound = range(points[, 1])[1],
                       xmaxBound = range(points[, 1])[2],
                       debug = FALSE) {

  # the number of known f and g values minus 1
  order <- sum(!is.na(points[, 2:3])) - 1

  # Code for most common case:
  #   - cubic interpolation of 2 points
  #       w/ function and derivative values for both
  if (nrow(points) == 2 && order == 3) {
    if (debug) {
      message("polyinterp common case")
    }
    # Solution in this case (where x2 is the farthest point):
    #    d1 <- g1 + g2 - 3*(f1-f2)/(x1-x2);
    #    d2 <- sqrt(d1^2 - g1*g2);
    #    minPos <- x2 - (x2 - x1)*((g2 + d2 - d1)/(g2 - g1 + 2*d2));
    #    t_new <- min(max(minPos,x1),x2);
    minPos <- which.min(points[, 1])
    notMinPos <- -minPos + 3

    x1 <- points[minPos, 1]
    x2 <- points[notMinPos, 1]
    f1 <- points[minPos, 2]
    f2 <- points[notMinPos, 2]
    g1 <- points[minPos, 3]
    g2 <- points[notMinPos, 3]

    if (x1 - x2 == 0) {
      return(x1)
    }

    d1 <- g1 + g2 - 3 * (f1 - f2) / (x1 - x2)
    d2sq <- d1 ^ 2 - g1 * g2

    if (is_finite_numeric(d2sq) && d2sq >= 0) {
      d2 <- sqrt(d2sq)

      x <- x2 - (x2 - x1) * ((g2 + d2 - d1) / (g2 - g1 + 2 * d2))
      if (debug) { message("d2 is real ", formatC(d2), " x = ", formatC(x)) }

      minPos <- min(max(x, xminBound), xmaxBound)
    }
    else {
      if (debug) { message("d2 is not real, bisecting") }

      minPos <- (xmaxBound + xminBound) / 2
    }

    return(minPos)
  }

  params <- polyfit(points)
  # If polynomial couldn't be found (due to singular matrix), bisect
  if (is.null(params)) {
    return((xminBound + xmaxBound) / 2)
  }

  # Compute Critical Points
  dParams <- rep(0, order)
  for (i in 1:order) {
    dParams[i] <- params[i + 1] * i
  }

  cp <- unique(c(xminBound, points[, 1], xmaxBound))

  # Remove mad, bad and dangerous to know critical points:
  # Must be finite, non-complex and not an extrapolation
  if (all(is.finite(dParams))) {
    cp <- c(cp,
            Re(Filter(
              function(x) {
                abs(Im(x)) < 1e-8 &&
                  Re(x) >= xminBound &&
                  Re(x) <= xmaxBound },
              polyroot(dParams))))
  }

  # Test Critical Points
  fcp <- polyval(cp, params)
  fminpos <- which.min(fcp)
  if (is.finite(fcp[fminpos])) {
    minpos <- cp[fminpos]
  }
  else {
    # Default to bisection if no critical points valid
    minpos <- (xminBound + xmaxBound) / 2
  }
  minpos
}

# Fits a polynomial to the known function and gradient values. The order of
# the polynomial is the number of known function and gradient values, minus one.
# points - an n x 3 matrix where n is the number of points and each row contains
#          x, f, g in columns 1-3 respectively.
# returns an array containing the coefficients of the polynomial in increasing
# order, e.g. c(1, 2, 3) is the polynomial 1 + 2x + 3x^2
# Returns NULL if the solution is singular
polyfit <- function(points) {
  nPoints <- nrow(points)
  # the number of known f and g values minus 1
  order <- sum(!is.na(points[, 2:3])) - 1

  # Constraints Based on available Function Values
  A <- NULL
  b <- NULL
  for (i in 1:nPoints) {
    if (!is.na(points[i, 2])) {
      constraint <- rep(0, order + 1)
      for (j in rev(0:order)) {
        constraint[order - j + 1] <- points[i, 1] ^ j
      }
      if (is.null(A)) {
        A <- constraint
      }
      else {
        A <- rbind(A, constraint)
      }
      if (is.null(b)) {
        b <- points[i, 2]
      }
      else {
        b <- c(b, points[i, 2])
      }
    }
  }

  # Constraints based on available Derivatives
  for (i in 1:nPoints) {
    if (!is.na(points[i, 3])) {
      constraint <- rep(0, order + 1)
      for (j in 1:order) {
        constraint[j] <- (order - j + 1) * points[i, 1] ^ (order - j)
      }
      if (is.null(A)) {
        A <- constraint
      }
      else {
        A <- rbind(A, constraint)
      }
      if (is.null(b)) {
        b <- points[i, 3]
      }
      else {
        b <- c(b, points[i, 3])
      }
    }
  }
  # Find interpolating polynomial
  params <- try(solve(A, b), silent = TRUE)
  if (class(params) == "numeric") {
    params <- rev(params)
  }
  else {
    params <- NULL
  }
}

# Evaluate 1D polynomial with coefs over the set of points x
# coefs - the coefficients for the terms of the polynomial ordered by
#   increasing degree, i.e. c(1, 2, 3, 4) represents the polynomial
#   4x^3 + 3x^2 + 2x + 1. This is the reverse of the ordering used by the Matlab
#   function, but is consistent with R functions like poly and polyroot
#   Also, the order of the arguments is reversed from the Matlab function
# Returns array of values of the evaluated polynomial
polyval <- function(x, coefs) {
  deg <- length(coefs) - 1
  # Sweep multiplies each column of the poly matrix by the coefficient
  rowSums(sweep(stats::poly(x, degree = deg, raw = TRUE),
                2, coefs[2:length(coefs)], `*`)) + coefs[1]
}

point_matrix <- function(xs, fs, gs) {
  matrix(c(xs, fs, gs), ncol = 3)
}

point_matrix_step <- function(step1, step2) {
  point_matrix(c(step1$alpha, step2$alpha), c(step1$f, step2$f),
               c(step1$d, step2$d))
}

