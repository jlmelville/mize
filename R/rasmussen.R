# Rasmussen Line Search
#
# Line Search Factory Function
#
# Line search algorithm originally written by Carl Edward Rasmussen in his
# conjugate gradient routine. It consists of two main parts:
# \enumerate{
#  \item Using cubic extrapolation from an initial starting guess for the step
#    size until either the sufficient decrease condition is not met or the
#    curvature condition is met.
#  \item Interpolation (quadratic or cubic) between that point and the start
#    point of the search until either a step size is found which meets the
#    Strong Wolfe conditions or the maximum number of allowed function
#    evaluations is reached.
# }
#
# The extrapolation and interpolation steps are bounded at each stage to ensure
# they don't represent too large or small a change to the step size.
#
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @param ext Extrapolation constant. Prevents step size extrapolation being
#   too large.
# @param int Interpolation constant. Prevents step size being too small.
# @param max_fn Maximum number of function evaluations allowed per line search.
# @return Line search function.
# @seealso Line search based on Matlab code by
#  \href{http://learning.eng.cam.ac.uk/carl/code/minimize/}{Carl Edward Rasmussen}
#  and also part of the Matlab
#  \href{(http://www.gaussianprocess.org/gpml/code/matlab/doc/)}{GPML} package.
rasmussen <- function(c1 = c2 / 2, c2 = 0.1, int = 0.1, ext = 3.0,
                      max_fn = Inf, xtol = 1e-6, eps = 1e-6, approx_armijo = FALSE,
                      strong_curvature = TRUE, verbose = FALSE) {
  if (c2 < c1) {
    stop("rasmussen line search: c2 < c1")
  }

  if (approx_armijo) {
    armijo_check_fn <- make_approx_armijo_ok_step(eps)
  }
  else {
    armijo_check_fn <- armijo_ok_step
  }

  wolfe_ok_step_fn <- make_wolfe_ok_step_fn(strong_curvature = strong_curvature,
                                            approx_armijo = approx_armijo,
                                            eps = eps)

  function(phi, step0, alpha,
           total_max_fn = Inf, total_max_gr = Inf, total_max_fg = Inf,
           pm = NULL) {
    maxfev <- min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    if (maxfev <= 0) {
      return(list(step = step0, nfn = 0, ngr = 0))
    }

    res <- ras_ls(phi, alpha, step0, c1 = c1, c2 = c2, ext = ext, int = int,
                  max_fn = maxfev, armijo_check_fn = armijo_check_fn,
                  wolfe_ok_step_fn = wolfe_ok_step_fn, verbose = verbose)
    list(step = res$step, nfn = res$nfn, ngr = res$nfn)
  }
}

# Rasmussen Line Search
#
# Line Search Method
#
# @param phi Line function.
# @param alpha Initial guess for step size.
# @param step0 Line search values at starting point of line search.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @param ext Extrapolation constant. Prevents step size extrapolation being
#   too large.
# @param int Interpolation constant. Prevents step size being too small.
# @param max_fn Maximum number of function evaluations allowed.
# @return List containing:
# \itemize{
#   \item step Valid step size or the last step size evaluated.
#   \item nfn Number of function evaluations.
# }
ras_ls <- function(phi, alpha, step0, c1 = 0.1, c2 = 0.1 / 2, ext = 3.0,
                   int = 0.1, max_fn = Inf, xtol = 1e-6,
                   armijo_check_fn = armijo_ok_step,
                   wolfe_ok_step_fn = strong_wolfe_ok_step,
                   verbose = verbose) {
  if (c2 < c1) {
    stop("Rasmussen line search: c2 < c1")
  }
  # extrapolate from initial alpha until either curvature condition is met
  # or the armijo condition is NOT met
  if (verbose) {
    message("Bracketing with initial step size = ", formatC(alpha))
  }
  ex_result <- extrapolate_step_size(phi, alpha, step0, c1, c2, ext, int,
                                     max_fn, armijo_check_fn, verbose = verbose)

  step <- ex_result$step
  nfn <- ex_result$nfn
  max_fn <- max_fn - nfn
  if (max_fn <= 0) {
    return(ex_result)
  }

  if (!ex_result$ok) {
    if (verbose) {
      message("Bracket phase failed, returning best step")
    }
    return(list(step = best_bracket_step(list(step0, step))), nfn = nfn)
  }

  if (verbose) {
    message("Bracket: ", format_bracket(list(step0, step)), " fn = ", nfn)
  }

  # interpolate until the Strong Wolfe conditions are met
  int_result <- interpolate_step_size(phi, step0, step, c1, c2, int, max_fn,
                                      xtol = xtol,
                                      armijo_check_fn = armijo_check_fn,
                                      wolfe_ok_step_fn = wolfe_ok_step_fn,
                                      verbose = verbose)
  if (verbose) {
    message("alpha = ", formatC(int_result$step$alpha))
  }
  int_result$nfn <- int_result$nfn + nfn
  int_result
}

# Increase Step Size
#
# @param phi Line function.
# @param alpha Initial step size.
# @param step0 Line search value at the initial step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @param ext Extrapolation constant. Prevents step size extrapolation being
#   too large.
# @param int Interpolation constant. Prevents step size extrapolation being
#   too small.
# @param max_fn Maximum number of function evaluations allowed.
# @return List containing:
# \itemize{
#   \item step Valid step size or the last step size evaluated.
#   \item nfn Number of function evaluations.
# }
extrapolate_step_size <- function(phi, alpha, step0, c1, c2, ext, int,
                                  max_fn = 20,
                                  armijo_check_fn = armijo_ok_step,
                                  verbose = FALSE) {
  # holds the largest finite-valued step
  finite_step <- step0
  ext_alpha <- alpha
  ok <- FALSE
  nfn <- 0
  while (TRUE) {
    result <- find_finite(phi, ext_alpha, max_fn, min_alpha = 0)
    nfn <- nfn + result$nfn
    max_fn <- max_fn - result$nfn
    if (!result$ok) {
      if (verbose) {
        message("Couldn't find a finite alpha during extrapolation")
      }
      break
    }

    finite_step <- result$step

    if (extrapolation_ok(step0, finite_step, c1, c2, armijo_check_fn)) {
      ok <- TRUE
      break
    }
    if (max_fn <= 0) {
      break
    }

    ext_alpha <- tweaked_extrapolation(step0, finite_step, ext, int)
  }

  list(step = finite_step, nfn = nfn, ok = ok)
}

# Extrapolation Check
#
# Checks that an extrapolated step size is sufficiently large: either by
# passing the curvature condition or failing the sufficient decrease condition.
#
# @param step0 Line search values at starting point of line search.
# @param step Line search value at candiate step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @return TRUE if the extrapolated point is sufficiently large.
extrapolation_ok <- function(step0, step, c1, c2, armijo_check_fn) {
  curvature_ok_step(step0, step, c2) || !armijo_check_fn(step0, step, c1)
}

# Extrapolate and Tweak Step Size
#
# Carries out an extrapolation of the step size, tweaked to not be too small
# or large.
#
# @param step0 Line search values at starting point of line search.
# @param step Line search value at candiate step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @return extrapolated step size.
tweaked_extrapolation <- function(step0, step, ext, int) {
  ext_alpha <- cubic_extrapolate_step(step0, step)
  tweak_extrapolation(ext_alpha, step0$alpha, step$alpha, ext, int)
}

# Interpolate Step Size to Meet Strong Wolfe Condition.
#
# @param phi Line function.
# @param step0 Line search values at starting point of line search.
# @param step Line search value at candiate step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @param int Interpolation constant. Prevents step size being too small.
# @param max_fn Maximum number of function evaluations allowed.
# @return List containing:
# \itemize{
#   \item step Valid step size or the last step size evaluated.
#   \item nfn Number of function evaluations.
# }
interpolate_step_size <- function(phi, step0, step, c1, c2, int, max_fn = 20,
                                  xtol = 1e-6,
                                  armijo_check_fn = armijo_ok_step,
                                  wolfe_ok_step_fn = strong_wolfe_ok_step,
                                  verbose = FALSE) {
  step2 <- step0
  step3 <- step
  nfn <- 0
  if (verbose) {
    message("Interpolating")
  }

  while (!wolfe_ok_step_fn(step0, step3, c1, c2) && nfn < max_fn) {
    if (step3$d > 0 || !armijo_check_fn(step0, step3, c1)) {
      step4 <- step3
    } else {
      step2 <- step3
    }

    if (step4$f > step0$f) {
      step3$alpha <- quadratic_interpolate_step(step2, step4)
    } else {
      step3$alpha <- cubic_interpolate_step(step2, step4)
    }

    if (verbose) {
      message("Bracket: ", format_bracket(list(step2, step4)),
              " alpha: ", formatC(step3$alpha), " f: ", formatC(step3$f),
              " d: ", formatC(step3$d), " nfn: ", nfn, " max_fn: ", max_fn)
    }
    step3$alpha <- tweak_interpolation(step3$alpha, step2$alpha, step4$alpha,
                                       int)
    # Check interpolated step is finite, and bisect if not, as in extrapolation
    # stage
    result <- find_finite(phi, step3$alpha, max_fn - nfn,
                          min_alpha = bracket_min_alpha(list(step2, step4)))
    nfn <- nfn + result$nfn
    if (!result$ok) {
      if (verbose) {
        message("Couldn't find a finite alpha during interpolation, aborting")
      }
      step3 <- best_bracket_step(list(step2, step4))
      break
    }
    step3 <- result$step


    if (bracket_width(list(step2, step4)) < xtol * step3$alpha) {
      if (verbose) {
        message("Bracket width: ", formatC(bracket_width(list(step2, step4))),
                " reduced below tolerance ", formatC(xtol * step3$alpha))
      }
      break
    }

  }
  list(step = step3, nfn = nfn)
}
