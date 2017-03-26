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
