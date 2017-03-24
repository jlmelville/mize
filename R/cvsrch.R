# More'-Thuente Line Search
#
# Line Search Factory Function
#
# Returns a line search function using a variant of the More-Thuente
#  line search originally implemented in
#  \href{http://www.netlib.org/minpack/}{MINPACK}.
#
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   \code{c1} and 1.
# @param max_fn Maximum number of function evaluations allowed.
# @return Line search function.
# @references
# More, J. J., & Thuente, D. J. (1994).
# Line search algorithms with guaranteed sufficient decrease.
# \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{20}(3),
# 286-307.
# @seealso This code is based on a translation of the original MINPACK code
#  for Matlab by
#  \href{https://www.cs.umd.edu/users/oleary/software/}{Dianne O'Leary}.
more_thuente <- function(c1 = 1e-4, c2 = 0.1, max_fn = Inf) {
  function(phi, step0, alpha,
           total_max_fn = Inf, total_max_gr = Inf, total_max_fg = Inf,
           pm = NULL) {
    maxfev <- min(max_fn, total_max_fn, total_max_gr, floor(total_max_fg / 2))
    if (maxfev <= 0) {
      return(list(step = step0, nfn = 0, ngr = 0))
    }
    res <- cvsrch(phi, step0, alpha = alpha, c1 = c1, c2 = c2,
                  maxfev = maxfev)
    list(step = res$step, nfn = res$nfn, ngr = res$nfn)
  }
}

# More'-Thuente Line Search
#
# This routine is a translation of Dianne O'Leary's Matlab code, which was
# itself a translation of the MINPACK original. Original comments to the Matlab
# code are at the end.
# @param phi Line function.
# @param step0 Line search values at starting point of line search.
# @param alpha Initial guess for step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @param xtol Relative width tolerance: convergence is reached if width falls
#   below xtol * maximum step size.
# @param alpha_min Smallest acceptable value of the step size.
# @param alpha_max Largest acceptable value of the step size.
# @param maxfev Maximum number of function evaluations allowed.
# @param delta Value to force sufficient decrease of interval size on
#   successive iterations. Should be a positive value less than 1.
# @return List containing:
# \itemize{
#   \item \code{step} Best step found and associated line search info.
#   \item \code{info} Return code from convergence check.
#   \item \code{nfn}  Number of function evaluations.
# }
#   Translation of minpack subroutine cvsrch
#   Dianne O'Leary   July 1991
#     **********
#
#     Subroutine cvsrch
#
#     The purpose of cvsrch is to find a step which satisfies
#     a sufficient decrease condition and a curvature condition.
#     The user must provide a subroutine which calculates the
#     function and the gradient.
#
#     At each stage the subroutine updates an interval of
#     uncertainty with endpoints stx and sty. The interval of
#     uncertainty is initially chosen so that it contains a
#     minimizer of the modified function
#
#          f(x+alpha*pv) - f(x) - c1*alpha*(gradf(x)'pv).
#
#     If a step is obtained for which the modified function
#     has a nonpositive function value and nonnegative derivative,
#     then the interval of uncertainty is chosen so that it
#     contains a minimizer of f(x+alpha*pv).
#
#     The algorithm is designed to find a step which satisfies
#     the sufficient decrease condition
#
#           f(x+alpha*pv) <= f(x) + c1*alpha*(gradf(x)'pv),
#
#     and the curvature condition
#
#           abs(gradf(x+alpha*pv)'pv)) <= c2*abs(gradf(x)'pv).
#
#     If c1 is less than c2 and if, for example, the function
#     is bounded below, then there is always a step which satisfies
#     both conditions. If no step can be found which satisfies both
#     conditions, then the algorithm usually stops when rounding
#     errors prevent further progress. In this case alpha only
#     satisfies the sufficient decrease condition.
#
#     The subroutine statement is
#
#        subroutine cvsrch(fcn,n,x,f,g,pv,alpha,c1,c2,xtol,
#                          alpha_min,alpha_max,maxfev,info,nfev)
#     where
#
#	fcn is the name of the user-supplied subroutine which
#         calculates the function and the gradient.  fcn must
#      	  be declared in an external statement in the user
#         calling program, and should be written as follows.
#
#         function [f,g] = fcn(n,x) (Matlab)     (10/2010 change in documentation)
#	  (derived from Fortran subroutine fcn(n,x,f,g) )
#         integer n
#         f
#         x(n),g(n)
#	  ----------
#         Calculate the function at x and
#         return this value in the variable f.
#         Calculate the gradient at x and
#         return this vector in g.
#	  ----------
#	  return
#	  end
#
#       n is a positive integer input variable set to the number
#	  of variables.
#
#	x is an array of length n. On input it must contain the
#	  base point for the line search. On output it contains
#         x + alpha*pv.
#
#	f is a variable. On input it must contain the value of f
#         at x. On output it contains the value of f at x + alpha*pv.
#
#	g is an array of length n. On input it must contain the
#         gradient of f at x. On output it contains the gradient
#         of f at x + alpha*pv.
#
#	pv is an input array of length n which specifies the
#         search direction.
#
#	alpha is a nonnegative variable. On input alpha contains an
#         initial estimate of a satisfactory step. On output
#         alpha contains the final estimate.
#
#       c1 and c2 are nonnegative input variables. Termination
#         occurs when the sufficient decrease condition and the
#         directional derivative condition are satisfied.
#
#	xtol is a nonnegative input variable. Termination occurs
#         when the relative width of the interval of uncertainty
#	  is at most xtol.
#
#	alpha_min and alpha_max are nonnegative input variables which
#	  specify lower and upper bounds for the step.
#
#	maxfev is a positive integer input variable. Termination
#         occurs when the number of calls to fcn is at least
#         maxfev by the end of an iteration.
#
#	info is an integer output variable set as follows:
#
#	  info = 0  Improper input parameters.
#
#	  info = 1  The sufficient decrease condition and the
#                   directional derivative condition hold.
#
#	  info = 2  Relative width of the interval of uncertainty
#		    is at most xtol.
#
#	  info = 3  Number of calls to fcn has reached maxfev.
#
#	  info = 4  The step is at the lower bound alpha_min.
#
#	  info = 5  The step is at the upper bound alpha_max.
#
#	  info = 6  Rounding errors prevent further progress.
#                   There may not be a step which satisfies the
#                   sufficient decrease and curvature conditions.
#                   Tolerances may be too small.
#
#       nfev is an integer output variable set to the number of
#         calls to fcn.
#
#
#     Subprograms called
#
#	user-supplied......fcn
#
#	MINPACK-supplied...cstep
#
#	FORTRAN-supplied...abs,max,min
#
#     Argonne National Laboratory. MINPACK Project. June 1983
#     Jorge J. More', David J. Thuente
#
#     **********
cvsrch <- function(phi, step0, alpha = 1,
                    c1 = 1e-4, c2 = 0.1, xtol = .Machine$double.eps,
                    alpha_min = 0, alpha_max = Inf,
                    maxfev = Inf, delta = 0.66) {

  xtrapf <- 4
  infoc <- 1

  # Check the input parameters for errors.
  params_ok <- TRUE
  problems <- c()
  if (alpha <= 0.0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("alpha <= 0.0: ", formatC(alpha)))
  }
  if (c1 < 0.0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("c1 < 0.0: ", formatC(c1)))
  }
  if (c2 < 0.0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("c2 < 0.0: ", formatC(c2)))
  }
  if (xtol < 0.0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("xtol < 0.0: ", formatC(xtol)))
  }
  if (alpha_min < 0.0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("alpha_min < 0.0: ", formatC(alpha_min)))
  }
  if (alpha_max < alpha_min) {
    params_ok <- FALSE
    problems <- c(problems, paste0("alpha_max ", formatC(alpha_max)
                         , " < alpha_min ", formatC(alpha_min)))
  }
  if (maxfev <= 0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("maxfev <= 0: ", formatC(maxfev)))
  }
  if (!params_ok) {
    problems <- paste(problems, collapse = "; ")
    stop("Parameter errors detected: ", problems)
  }

  # Check that pv is a descent direction.
  d0 <- step0$d
  if (step0$d >= 0.0) {
    stop("Not a descent direction")
  }
  dgtest <- c1 * step0$d

  # Initialize local variables.
  brackt <- FALSE
  stage1 <- TRUE
  nfev <- 0

  width <- alpha_max - alpha_min
  width_old <- 2 * width

  # The variables stx, fx, dgx contain the values of the step,
  # function, and directional derivative at the best step.
  # The variables sty, fy, dgy contain the value of the step,
  # function, and derivative at the other endpoint of
  # the interval of uncertainty.
  # The variables alpha, f, dg contain the values of the step,
  # function, and derivative at the current step.
  stepx <- step0
  stepy <- step0
  step <- list(alpha = alpha)

  #     Start of iteration.
  iter <- 0
  while (1) {
    iter <- iter + 1
    # Set the minimum and maximum steps to correspond
    # to the present interval of uncertainty.
    if (brackt) {
      stmin <- min(stepx$alpha, stepy$alpha)
      stmax <- max(stepx$alpha, stepy$alpha)
    } else {
      stmin <- stepx$alpha
      stmax <- step$alpha + xtrapf * (step$alpha - stepx$alpha)
    }

    # Force the step to be within the bounds alpha_max and alpha_min.
    step$alpha <- max(step$alpha, alpha_min)
    step$alpha <- min(step$alpha, alpha_max)

    # Evaluate the function and gradient at alpha
    # and compute the directional derivative.
    step <- phi(step$alpha)
    nfev <- nfev + 1
    # Test for convergence.
    info <- check_convergence(step0, step, brackt, infoc, stmin, stmax,
                              alpha_min, alpha_max, c1, c2, dgtest, nfev,
                              maxfev, xtol)
    # Check for termination.
    if (info != 0) {
      # If an unusual termination is to occur, then set step to the best step
      # found
      if (info == 2 || info == 3 || info == 6) {
        step <- stepx
      }
      return(list(step = step, info = info, nfn = nfev))
    }

    # In the first stage we seek a step for which the modified
    # function has a nonpositive value and nonnegative derivative.
    if (stage1 && armijo_ok_step(step0, step, c1) && step$df >= min(c1, c2) * d0) {
      stage1 <- FALSE
    }

    # A modified function is used to predict the step only if
    # we have not obtained a step for which the modified
    # function has a nonpositive function value and nonnegative
    # derivative, and if a lower function value has been
    # obtained but the decrease is not sufficient.
    if (stage1 && step$f <= stepx$f && !armijo_ok_step(step0, step, c1)) {
      # Define the modified function and derivative values.

      stepxm <- modify_step(stepx, dgtest)
      stepym <- modify_step(stepy, dgtest)
      stepm <- modify_step(step, dgtest)

      step_result <- cstep(stepxm, stepym, stepm, brackt, stmin, stmax)

      brackt <- step_result$brackt
      infoc <- step_result$info
      stepxm <- step_result$stepx
      stepym <- step_result$stepy
      stepm <- step_result$step

      # Reset the function and gradient values for f.
      stepx <- unmodify_step(stepxm, dgtest)
      stepy <- unmodify_step(stepym, dgtest)
      step$alpha <- stepm$alpha
    } else {
      # Call cstep to update the interval of uncertainty
      # and to compute the new step.
      step_result <- cstep(stepx, stepy, step, brackt, stmin, stmax)
      brackt <- step_result$brackt
      infoc <- step_result$info
      stepx <- step_result$stepx
      stepy <- step_result$stepy
      step <- step_result$step
    }

    # Force a sufficient decrease in the size of the interval of uncertainty.
    if (brackt) {
      # if the length of I does not decrease by a factor of delta < 1
      # then use a bisection step for the next trial alpha
      width_new <- abs(stepy$alpha - stepx$alpha)
      if (width_new >= delta * width_old) {
        step$alpha <- stepx$alpha + 0.5 * (stepy$alpha - stepx$alpha)
      }
      width_old <- width
      width <- width_new
    }
  }
}


# Modify Line Search Values
#
# Modifies a line search function and directional derivative value.
# Used by MINPACK version of More'-Thuente line search algorithm.
#
# @param step Line search information.
# @param dgtest Product of the initial line search directional derivative and
#   the sufficent decrease condition constant.
# @return Modified step size.
modify_step <- function(step, dgtest) {
  stepm <- step
  stepm$f <- step$f - step$alpha * dgtest
  stepm$d <- step$d - dgtest
  stepm
}

# Un-modify Line Search Values
#
# Un-modifies a line search function and directional derivative value that was
# modified by the modify_step function. Used by MINPACK version of More'-Thuente
# line search algorithm.
#
# @param stepm Modified line search information.
# @param dgtest Product of the initial line search directional derivative and
#   the sufficent decrease condition constant.
# @return Unmodified step size.
unmodify_step <- function(stepm, dgtest) {
  stepm$f <- stepm$f + stepm$alpha * dgtest
  stepm$d <- stepm$d + dgtest
  stepm
}

# Check Convergence of More'-Thuente Line Search
#
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param brackt TRUE if the step has been bracketed.
# @param infoc Return code of the last step size update.
# @param stmin Smallest value of the step size interval.
# @param stmax Largest value of the step size interval.
# @param alpha_min Smallest acceptable value of the step size.
# @param alpha_max Largest acceptable value of the step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @param dgtest Product of the initial line search directional derivative and
#   the sufficent decrease condition constant.
# @param nfev Current number of function evaluations.
# @param maxfev Maximum number of function evaluations allowed.
# @param xtol Relative width tolerance: convergence is reached if width falls
#   below xtol * stmax.
# @return Integer code indicating convergence state:
#  \itemize{
#   \item \code{0} No convergence.
#   \item \code{1} The sufficient decrease condition and the directional
#     derivative condition hold.
#	  \item \code{2} Relative width of the interval of uncertainty
#		    is at most xtol.
#	  \item \code{3} Number of calls to fcn has reached maxfev.
#	  \item \code{4} The step is at the lower bound alpha_min.
#	  \item \code{5} The step is at the upper bound alpha_max.
#	  \item \code{6} Rounding errors prevent further progress.
# }
check_convergence <- function(step0, step, brackt, infoc, stmin, stmax,
                              alpha_min, alpha_max, c1, c2, dgtest, nfev,
                              maxfev, xtol) {
  info <- 0
  if (!is.finite(step$f) || any(is.nan(step$df))) {
    return(6)
  }
  if ((brackt && (step$alpha <= stmin || step$alpha >= stmax)) || infoc == 0) {
    # rounding errors prevent further progress
    info <- 6
  }
  if (step$alpha == alpha_max && armijo_ok_step(step0, step, c1) && step$d <= dgtest) {
    # reached alpha_max
    info <- 5
  }
  if (step$alpha == alpha_min && (!armijo_ok_step(step0, step, c1) || step$d >= dgtest)) {
    # reached alpha_min
    info <- 4
  }
  if (nfev >= maxfev) {
    # maximum number of function evaluations reached
    info <- 3
  }
  if (brackt && stmax - stmin <= xtol * stmax) {
    # interval width is below xtol
    info <- 2
  }
  if (strong_wolfe_ok_step(step0, step, c1, c2)) {
    # success!
    info <- 1
  }
  info
}
