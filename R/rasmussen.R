#' Rasmussen Line Search
#'
#' Line Search Factory Function
#'
#' Line search algorithm originally written by Carl Edward Rasmussen in his
#' conjugate gradient routine. It consists of two main parts:
#' \enumerate{
#'  \item Using cubic extrapolation from an initial starting guess for the step
#'    size until either the sufficient decrease condition is not met or the
#'    strong curvature condition is met.
#'  \item Interpolation (quadratic or cubic) between that point and the start
#'    point of the search until either a step size is found which meets the
#'    Strong Wolfe conditions or the maximum number of allowed function
#'    evaluations is reached.
#' }
#'
#' The extrapolation and interpolation steps are bounded at each stage to ensure
#' they don't represent too large or small a change to the step size.
#'
#' @param c1 Constant used in sufficient decrease condition. Should take a value
#'   between 0 and 1.
#' @param c2 Constant used in curvature condition. Should take a value between
#'   c1 and 1.
#' @param ext Extrapolation constant. Prevents step size extrapolation being
#'   too large.
#' @param int Interpolation constant. Prevents step size being too small.
#' @param max_fn Maximum number of function evaluations allowed.
#' @return Line search function.
#' @seealso Line search based on Matlab code by
#'  \href{http://learning.eng.cam.ac.uk/carl/code/minimize/}{Carl Edward Rasmussen}
#'  and also part of the Matlab
#'  \href{(http://www.gaussianprocess.org/gpml/code/matlab/doc/)}{GPML} package.
rasmussen <- function(c1 = c2 / 2, c2 = 0.1, int = 0.1, ext = 3.0,
                      max_fn = 20) {
  if (c2 < c1) {
    stop("rasmussen line search: c2 < c1")
  }
  function(phi, step0, alpha) {
    res <- ras_ls(phi, alpha, step0, c1 = c1, c2 = c2, ext = ext, int = int,
                  max_fn = max_fn)
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
                   int = 0.1, max_fn = Inf) {
  if (c2 < c1) {
    stop("Rasmussen line search: c2 < c1")
  }

  # extrapolate from initial alpha until either curvature condition is met
  # or the armijo condition is NOT met
  ex_result <- extrapolate_step_size(phi, alpha, step0, c1, c2, ext, int,
                                     max_fn)

  step <- ex_result$step
  nfn <- ex_result$nfn
  max_fn <- max_fn - nfn
  if (max_fn <= 0) {
    return(ex_result)
  }

  # interpolate until the Strong Wolfe conditions are met
  int_result <- interpolate_step_size(phi, step0, step, c1, c2, int, max_fn)
  int_result$nfn <- int_result$nfn + nfn
  int_result
}

# Ensure Valid Step Size
#
# Given an initial step size, if either the function value or the directional
# derivative is non-finite (NaN or infinite), reduce the step size until
# finite values are found.
#
# @param phi Line function.
# @param alpha Initial step size.
# @param min_alpha Minimum step size.
# @param max_fn Maximum number of function evaluations allowed.
# @return List containing:
# \itemize{
#   \item step Valid step size or the last step size evaluated.
#   \item nfn Number of function evaluations.
# }
find_finite <- function(phi, alpha, min_alpha = 0, max_fn = 20) {
  nfn <- 0
  while (nfn < max_fn && alpha > min_alpha) {
    step <- phi(alpha)
    nfn <- nfn + 1
    if (is.finite(step$f) && is.finite(step$df)) {
      break
    }
    alpha <- (min_alpha + alpha) / 2
  }
  list(step = step, nfn = nfn)
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
                                  max_fn = 20) {
  step <- list(alpha = alpha)

  nfn <- 0
  while (1) {
    result <- find_finite(phi, step$alpha, max_fn, min_alpha = 0)
    nfn <- nfn + result$nfn
    max_fn <- max_fn - result$nfn
    step <- result$step

    if (extrapolation_ok(step0, step, c1, c2) || max_fn == 0) {
      break
    }
    step$alpha <- tweaked_extrapolation(step0, step, ext, int)
  }

  list(step = step, nfn = nfn)
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
extrapolation_ok <- function(step0, step, c1, c2) {
  curvature_ok_step(step0, step, c2) || !armijo_ok_step(step0, step, c1)
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
interpolate_step_size <- function(phi, step0, step, c1, c2, int, max_fn = 20) {
  step2 <- step0
  step3 <- step
  nfn <- 0

  while (!strong_wolfe_ok_step(step0, step3, c1, c2) && nfn < max_fn) {
    if (step3$d > 0 || !armijo_ok_step(step0, step3, c1)) {
      step4 <- step3
    } else {
      step2 <- step3
    }

    if (step4$f > step0$f) {
      step3$alpha <- quadratic_interpolate_step(step2, step4)
    } else {
      step3$alpha <- cubic_interpolate_step(step2, step4)
    }

    step3$alpha <- tweak_interpolation(step3$alpha, step2$alpha, step4$alpha, int)
    step3 <- phi(step3$alpha)
    nfn <- nfn + 1
  }
  list(step = step3, nfn = nfn)
}



