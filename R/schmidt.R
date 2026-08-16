# Translation of Mark Schmidt's minFunc line search code for satisfying the
# Strong Wolfe conditions (and also the Armijo conditions)
# http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.

# Adapters ----------------------------------

# Supported Schmidt Wolfe profile over the shared bracket-and-zoom engine.
schmidt_core <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  direction,
  method_policy
) {
  run_bracket_zoom(
    evaluator = evaluator,
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    condition_policy = condition_policy,
    direction = direction,
    policy = method_policy
  )
}

# step_down if non-NULL, multiply the step size by this value when backtracking
# Otherwise, use a cubic interpolation based on previous function and derivative
# values
schmidt_armijo_backtrack <- function(
  c1 = 0.05,
  step_down = NULL,
  max_fn = Inf
) {
  function(
    phi,
    step0,
    alpha,
    total_max_fn = Inf,
    total_max_gr = Inf,
    total_max_fg = Inf,
    pm
  ) {
    maxfev <- min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    if (maxfev <= 0) {
      return(list(step = step0, nfn = 0, ngr = 0))
    }

    if (!is.null(step_down)) {
      # fixed-size step reduction by a factor of step_down
      LS_interp <- 0
    } else {
      # cubic interpolation
      LS_interp <- 2
    }

    candidate_is_finite <- function(candidate) {
      isTRUE(is.finite(candidate$alpha)) &&
        step_is_finite(candidate) &&
        (is.null(candidate$par) || all(is.finite(candidate$par)))
    }

    best_decrease <- NULL
    last_candidate <- NULL
    tracked_phi <- function(...) {
      candidate <- phi(...)
      last_candidate <<- candidate
      if (
        candidate_is_finite(candidate) &&
          isTRUE(candidate$f < step0$f) &&
          (is.null(best_decrease) || candidate$f < best_decrease$f)
      ) {
        best_decrease <<- candidate
      }
      candidate
    }

    res <- ArmijoBacktrack(
      step = alpha,
      f = step0$f,
      g = step0$df,
      gtd = step0$d,
      c1 = c1,
      LS_interp = LS_interp,
      LS_multi = 0,
      maxLS = maxfev,
      step_down = step_down,
      funObj = tracked_phi,
      varargin = NULL,
      pnorm_inf = max(abs(pm)),
      progTol = 1e-9,
      debug = FALSE
    )

    # ArmijoBacktrack omits par from its returned step, so check the matching
    # evaluated candidate retained by the adapter.
    returned_trial_is_safe <- !is.null(last_candidate) &&
      isTRUE(last_candidate$alpha == res$step$alpha) &&
      candidate_is_finite(last_candidate) &&
      isTRUE(last_candidate$f < step0$f) &&
      isTRUE(armijo_ok_step(step0, last_candidate, c1))
    if (!isTRUE(res$step$alpha == 0) && !returned_trial_is_safe) {
      if (is.null(best_decrease)) {
        res$step <- step0
      } else {
        res$step <- best_decrease
      }
      res$is_gr_curr <- !is.null(res$step$df)
    }

    res
  }
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
# recent change: LS changed to LS_interp and LS_multi
ArmijoBacktrack <-
  function(
    step,
    f,
    g,
    gtd,
    c1 = 1e-4,
    LS_interp = 2,
    LS_multi = 0,
    maxLS = Inf,
    step_down = 0.5,
    funObj,
    varargin = NULL,
    pnorm_inf,
    progTol = 1e-9,
    debug = TRUE
  ) {
    # Evaluate the Objective and Gradient at the Initial Step
    f_prev <- NA
    t_prev <- NA
    g_prev <- NA
    gtd_prev <- NA

    # Don't calculate gradient if we aren't using gradient values in the
    # interpolation
    calc_gradient <- LS_interp != 0
    fun_obj_res <- funObj(step, calc_gradient = calc_gradient)
    f_new <- fun_obj_res$f
    g_new <- fun_obj_res$df
    gtd_new <- fun_obj_res$d
    funEvals <- 1
    grEvals <- ifelse(calc_gradient, 1, 0)

    while (
      funEvals < maxLS && (f_new > f + c1 * step * gtd || !is.finite(f_new))
    ) {
      temp <- step
      if (LS_interp == 0 || !is.finite(f_new)) {
        # Ignore value of new point
        if (debug) {
          message("Fixed BT")
        }
        step <- step_down * step
      } else if (LS_interp == 1 || any(!is.finite(g_new))) {
        # Use function value at new point, but not its derivative
        if (funEvals < 2 || LS_multi == 0 || !is.finite(f_prev)) {
          # Backtracking w/ quadratic interpolation based on two points
          if (debug) {
            message("Quad BT")
          }
          step <- polyinterp(
            point_matrix(c(0, step), c(f, f_new), c(gtd, NA)),
            0,
            step
          )
        } else {
          # Backtracking w/ cubic interpolation based on three points
          if (debug) {
            message("Cubic BT")
          }
          step <-
            polyinterp(
              point_matrix(
                c(0, step, t_prev),
                c(f, f_new, f_prev),
                c(gtd, NA, NA)
              ),
              0,
              step
            )
        }
      } else {
        # Use function value and derivative at new point
        if (funEvals < 2 || LS_multi == 0 || !is.finite(f_prev)) {
          # Backtracking w/ cubic interpolation w/ derivative
          if (debug) {
            message("Grad-Cubic BT")
          }
          step <- polyinterp(
            point_matrix(c(0, step), c(f, f_new), c(gtd, gtd_new)),
            0,
            step
          )
        } else if (any(!is.finite(g_prev))) {
          # Backtracking w/ quartic interpolation 3 points and derivative
          # of two
          if (debug) {
            message("Grad-Quartic BT")
          }

          step <- polyinterp(
            point_matrix(
              c(0, step, t_prev),
              c(f, f_new, f_prev),
              c(gtd, gtd_new, NA)
            ),
            0,
            step
          )
        } else {
          # Backtracking w/ quintic interpolation of 3 points and derivative
          # of three
          if (debug) {
            message("Grad-Quintic BT")
          }

          step <- polyinterp(
            point_matrix(
              c(0, step, t_prev),
              c(f, f_new, f_prev),
              c(gtd, gtd_new, gtd_prev)
            ),
            0,
            step
          )
        }
      }

      if (!is_finite_numeric(step)) {
        step <- temp * 0.6
      }
      # Adjust if change in step is too small/large
      if (step < temp * 1e-3) {
        if (debug) {
          message("Interpolated Value Too Small, Adjusting")
        }
        step <- temp * 1e-3
      } else if (step > temp * 0.6) {
        if (debug) {
          message("Interpolated Value Too Large, Adjusting")
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

      fun_obj_res <- funObj(step, calc_gradient = calc_gradient)
      f_new <- fun_obj_res$f
      g_new <- fun_obj_res$df
      gtd_new <- fun_obj_res$d
      funEvals <- funEvals + 1
      grEvals <- ifelse(calc_gradient, grEvals + 1, grEvals)

      # Check whether step size has become too small
      if (pnorm_inf * step <= progTol) {
        if (debug) {
          message("Backtracking Line Search Failed")
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
      nfn = funEvals,
      ngr = grEvals,
      is_gr_curr = calc_gradient
    )
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
polyinterp <- function(
  points,
  xminBound = range(points[, 1])[1],
  xmaxBound = range(points[, 1])[2],
  debug = FALSE
) {
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
    d2sq <- d1^2 - g1 * g2

    if (is_finite_numeric(d2sq) && d2sq >= 0) {
      d2 <- sqrt(d2sq)

      x <- x2 - (x2 - x1) * ((g2 + d2 - d1) / (g2 - g1 + 2 * d2))
      if (debug) {
        message("d2 is real ", formatC(d2), " x = ", formatC(x))
      }

      minPos <- min(max(x, xminBound), xmaxBound)
    } else {
      if (debug) {
        message("d2 is not real, bisecting")
      }

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
    cp <- c(
      cp,
      Re(Filter(
        function(x) {
          abs(Im(x)) < 1e-8 &&
            Re(x) >= xminBound &&
            Re(x) <= xmaxBound
        },
        polyroot(dParams)
      ))
    )
  }

  # Test Critical Points
  fcp <- polyval(cp, params)
  fminpos <- which.min(fcp)
  if (is.finite(fcp[fminpos])) {
    minpos <- cp[fminpos]
  } else {
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
        constraint[order - j + 1] <- points[i, 1]^j
      }
      if (is.null(A)) {
        A <- constraint
      } else {
        A <- rbind(A, constraint)
      }
      if (is.null(b)) {
        b <- points[i, 2]
      } else {
        b <- c(b, points[i, 2])
      }
    }
  }

  # Constraints based on available Derivatives
  for (i in 1:nPoints) {
    if (!is.na(points[i, 3])) {
      constraint <- rep(0, order + 1)
      for (j in 1:order) {
        constraint[j] <- (order - j + 1) * points[i, 1]^(order - j)
      }
      if (is.null(A)) {
        A <- constraint
      } else {
        A <- rbind(A, constraint)
      }
      if (is.null(b)) {
        b <- points[i, 3]
      } else {
        b <- c(b, points[i, 3])
      }
    }
  }
  # Find interpolating polynomial
  params <- try(solve(A, b), silent = TRUE)
  if (methods::is(params, "numeric")) {
    params <- rev(params)
  } else {
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
  rowSums(sweep(
    stats::poly(x, degree = deg, raw = TRUE),
    2,
    coefs[2:length(coefs)],
    `*`
  )) +
    coefs[1]
}

point_matrix <- function(xs, fs, gs) {
  matrix(c(xs, fs, gs), ncol = 3)
}
