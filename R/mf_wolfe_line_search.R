# Translation of Mark Schmidt's minFunc line search
# http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.

schmidt <- function(c1 = c2 / 2, c2 = 0.1, max_fn = Inf) {
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
    res <- WolfeLineSearch(alpha = step0$alpha, f = step0$f, g = step0$df,
                           gtd = step0$d,
                           c1 = 1e-4, c2 = 0.1, LS_interp = 2, LS_multi = 0,
                           maxLS = 25,
                           funObj = phi, varargin = NULL,
                           pnorm_inf = max(abs(pm)),
                           progTol = 1e-9, debug = FALSE)
    list(step = list(
      alpha = res$step,
      f = res$f_new,
      df = res$g_new,
      d = res$d_new),
      nfn = res$funEvals, ngr = res$funEvals)
    }
}

#
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
#   maxLS: maximum number of iterations
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
           pnorm_inf, progTol = 1e-9, debug = FALSE) {

    # Evaluate the Objective and Gradient at the Initial Step
    fun_obj_res <- funObj(alpha)
    f_new <- fun_obj_res$f
    g_new <- fun_obj_res$df
    gtd_new <- fun_obj_res$d
    step_new <- list(alpha = alpha, f = f_new, df = g_new, d = gtd_new)

    funEvals <- 1

    # Bracket an Interval containing a point satisfying the
    # Wolfe criteria

    done <- FALSE

    LSiter <- 0
    t_prev <- 0
    f_prev <- f
    g_prev <- g
    gtd_prev <- gtd

    step0 <- list(alpha = 0, f = f, df = g, d = gtd_prev)
    step_prev <- step0

    # Bracketing Phase
    if (debug) {
      message("% BRACKET PHASE %")
    }
    while (LSiter < maxLS) {

      if (!is.finite(f_new) || !is.finite(g_new)) {
        if (debug) {
          message('Extrapolated into illegal region, switching to Armijo line-search')
        }
        alpha <- (alpha + t_prev) / 2

        # Do Armijo
        armijo_res <- ArmijoBacktrack(alpha, f, g, gtd,
                          c1, LS_interp, LS_multi,
                          funObj, varargin,
                          progTol, debug)

        alpha <- armijo_res$t
        f_new <- armijo_res$f_new
        g_new <- armijo_res$g_new
        armijoFunEvals <- armijo_res$armijoFunEvals

        funEvals <- funEvals + armijoFunEvals
        return(list(
          step = list(alpha = alpha, f = f_new, df = g_new),
          nfn = funEvals
        ))
      }

      if (f_new > f + c1 * alpha * gtd || (LSiter > 1 && f_new >= f_prev)) {
        bracket <- c(t_prev, alpha)
        bracketFval <- c(f_prev, f_new)
        # bracketGval <- c(g_prev, g_new)
        bracketDval <- c(gtd_prev, gtd_new)

        bracket_step <- list(step_prev, step_new)

        if (debug) {
          message('Armijo failed or f_new >= f_prev: bracket is [prev new]')
        }
        break
      }
      else if (abs(gtd_new) <= -c2 * gtd) {
        bracket <- alpha
        bracketFval <- c(f_new)
        # bracketGval <- c(g_new)
        bracketDval <- c(gtd_new)
        done <- TRUE

        bracket_step <- list(step_new)
        if (debug) {
          message('Sufficient curvature found: bracket is [new]')
        }
        break
      }
      else if (gtd_new >= 0) {
        bracket <- c(t_prev, alpha)
        bracketFval <- c(f_prev, f_new)
        # bracketGval <- c(g_prev, g_new)
        bracketDval <- c(gtd_prev, gtd_new)

        bracket_step <- list(step_prev, step_new)
        if (debug) {
          message('gtd_new >= 0: bracket is [prev, new]')
        }
        break
      }

      minStep <- alpha + 0.01 * (alpha - t_prev)
      maxStep <- alpha * 10


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
        alpha_new <- polyinterp(point_matrix(
          c(t_prev, alpha), c(f_prev, f_new), c(gtd_prev, gtd_new)),
          minStep, maxStep)
      }
      else {
        # LS_interp == 3
        alpha_new <- mixedExtrap(t_prev, f_prev, gtd_prev,
                            alpha, f_new,  gtd_new, minStep, maxStep, debug)
      }

      t_prev <- alpha
      alpha <- alpha_new
      f_prev <- f_new
      g_prev <- g_new
      gtd_prev <- gtd_new

      step_prev <- step_new

      fun_obj_res <- funObj(alpha)
      f_new <- fun_obj_res$f
      g_new <- fun_obj_res$df
      gtd_new <- fun_obj_res$d

      step_new <- list(alpha = alpha, f = f_new, df = g_new, d = gtd_new)


      # fun_obj_res <- funObj(x + step * d, varargin)
      # f_new <- fun_obj_res$f_new
      # g_new <- fun_obj_res$g_new

      funEvals <- funEvals + 1
      LSiter <- LSiter + 1
    }

    if (LSiter == maxLS) {
      bracket <- c(0, alpha)
      bracketFval <- c(f, f_new)
      # bracketGval <- c(g, g_new)
      bracketDval <- c(gtd, gtd_new)
      bracket_step <- list(step0, step_new)
    }

    ## Zoom Phase
    # We now either have a point satisfying the criteria, or a bracket
    # surrounding a point satisfying the criteria
    # Refine the bracket until we find a point satisfying the criteria

    if (debug) {
      message("% ZOOM PHASE %")
    }

    insufProgress <- FALSE
    Tpos <- 2
    LOposRemoved <- FALSE
    if (debug) {
      if (length(bracket_step) == 2) {
        message("bracket alpha = [", formatC(bracket_step[[1]]$alpha), ", ",
              formatC(bracket_step[[2]]$alpha), "]")
      }
      else {
        message("bracket alpha = [", formatC(bracket_step[[1]]$alpha), "]")
      }
    }
    while (!done && LSiter < maxLS) {
      # Find High and Low Points in bracket
      LOpos <- which.min(bracketFval)
      f_LO <- bracketFval[LOpos]
      HIpos <- -LOpos + 3 # 1 or 2, whichever wasn't the LOpos

      # Compute new trial value
      if (LS_interp <= 1 || !all(is.finite(c(bracketFval, bracketDval)))) {
        alpha <- mean(bracket)
        if (debug) {
          message('Bisecting: trial step = ', formatC(alpha))
        }

      }
      else if (LS_interp == 2) {
        alpha <- polyinterp(point_matrix(
          c(bracket[1], bracket[2]),
          c(bracketFval[1], bracketFval[2]),
          c(bracketDval[1], bracketDval[2])), debug = debug)
        if (debug) {
          message('Grad-Cubic Interpolation: trial step = ', formatC(alpha),
                   ' a1 = ', formatC(bracket[1]), ' a2 = ', formatC(bracket[2]),
                   ' f1 = ', formatC(bracketFval[1]), ' f2 = ', formatC(bracketFval[2]),
                   ' d1 = ', formatC(bracketDval[1]), ' d2 = ', formatC(bracketDval[2])
                  )
        }
      }
      else {
        # Mixed Case #
        nonTpos <- -Tpos + 3
        if (!LOposRemoved) {
          oldLOval <- bracket[nonTpos]
          oldLOFval <- bracketFval[nonTpos]
          # oldLOGval <- bracketGval[, nonTpos]
          oldLODval <- bracketDval[nonTpos]
        }

        alpha <- mixedInterp(bracket, bracketFval, bracketDval,
                            Tpos,
                            oldLOval, oldLOFval, oldLODval, debug)
        if (debug) {
          message('Mixed Interpolation: trial step = ', formatC(alpha))
        }
      }

      # Test that we are making sufficient progress
      if (min(max(bracket) - alpha, alpha - min(bracket)) / (max(bracket) - min(bracket)) < 0.1) {
        if (debug) {
          message('Interpolation close to boundary')
        }

        if (insufProgress || alpha >= max(bracket) || alpha <= min(bracket)) {
          if (debug) {
            message('Evaluating at 0.1 away from boundary')
          }
          if (abs(alpha - max(bracket)) < abs(alpha - min(bracket))) {
            alpha <- max(bracket) - 0.1 * (max(bracket) - min(bracket))
          }
          else {
            alpha <- min(bracket) + 0.1 * (max(bracket) - min(bracket))
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

      # Evaluate new point
      fun_obj_res <- funObj(alpha)
      f_new <- fun_obj_res$f
      g_new <- fun_obj_res$df
      gtd_new <- fun_obj_res$d

      # fun_obj_res <- funObj(x + step * d, varargin)
      # f_new <- fun_obj_res$f_new
      # g_new <- fun_obj_res$g_new

      funEvals <- funEvals + 1
      LSiter <- LSiter + 1

      armijo <- f_new < f + c1 * alpha * gtd
      if (debug) {
        message('Armijo check f_new ', f_new,
                ' < f + c1 * t * gtd = ',
                formatC(f), ' + ', formatC(c1), ' * ', formatC(alpha), ' * ',
                formatC(gtd), ' = ', formatC(f + c1 * alpha * gtd), ' ? ')
      }
      if (!armijo || f_new >= f_LO) {
        if (debug) {
          if (!armijo) {
            message("Armijo not satisfied")
          }
          if (f_new >= f_LO) {
            message('f_new >= f_LO ', f_new, " ", f_LO)
          }
        }
        if (debug) {
          message("New point becomes new HI")
        }
        # Armijo condition not satisfied or not lower than lowest point
        bracket[HIpos] <- alpha
        bracketFval[HIpos] <- f_new
        # bracketGval[, HIpos] <- g_new
        bracketDval[HIpos] <- gtd_new
        Tpos <- HIpos
      }
      else {
        if (debug) {
          message("Wolfe conditions check |gtd_new| ",
                  formatC(abs(gtd_new)), " <= -c2 * gtd = ",
                  formatC(-c2), " * " , formatC(gtd), " = ",
                  formatC(-c2 * gtd), " ?")
        }
        if (abs(gtd_new) <= -c2 * gtd) {
          if (debug) {
            message("Wolfe conditions satisfied, done!")
          }
          # Wolfe conditions satisfied
          done <- TRUE
        }
        else if (gtd_new * (bracket[HIpos] - bracket[LOpos]) >= 0) {
          if (debug) {
            message("Old HI becomes new LO")
          }
          # Old HI becomes new LO
          bracket[HIpos] <- bracket[LOpos]
          bracketFval[HIpos] <- bracketFval[LOpos]
          # bracketGval[, HIpos] <- bracketGval[, LOpos]
          bracketDval[HIpos] <- bracketDval[LOpos]
          if (LS_interp == 3) {
            if (debug) {
              message('LO Pos is being removed!')
            }
            LOposRemoved <- TRUE
            oldLOval <- bracket[LOpos]
            oldLOFval <- bracketFval[LOpos]
            # oldLOGval <- bracketGval[, LOpos]
            oldLODval <- bracketDval[LOpos]
          }
        }
        if (debug) {
          message("New point becomes new LO")
        }
        # New point becomes new LO
        bracket[LOpos] <- alpha
        bracketFval[LOpos] <- f_new
        # bracketGval[, LOpos] <- g_new
        bracketDval[LOpos] <- gtd_new
        Tpos <- LOpos
      }

      if (debug) {
        message('Bracket ', formatC(bracket[1]), " ",
                formatC(bracket[2]), " width = ",
                formatC(abs(bracket[1] - bracket[2])))
      }


      if (!done && abs(bracket[1] - bracket[2]) * pnorm_inf < progTol) {
        if (debug) {
          message('Line-search bracket has been reduced below progTol')
        }
        break
      }
    } # end of while loop

    ##
    if (LSiter == maxLS) {
      if (debug) {
        message('Line Search Exceeded Maximum Line Search Iterations')
      }
    }

    LOpos <- which.min(bracketFval)
    f_LO <- bracketFval[LOpos]
    alpha <- bracket[LOpos]
    f_new <- bracketFval[LOpos]
    # g_new <- bracketGval[, LOpos]
    d_new <- bracketDval[LOpos]

    list(
      step = list(alpha = alpha, f = f_new, df = g_new, d = d_new),
      nfn = funEvals
    )
  }
