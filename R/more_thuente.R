# More'-Thuente Line Search
#
# Combination of the cvsrch and cstep matlab files.
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
#  \href{http://www.cs.umd.edu/users/oleary/software/}{Dianne O'Leary}.
more_thuente <- function(c1 = 1e-4, c2 = 0.1, max_fn = Inf, eps = 1e-6,
                         approx_armijo = FALSE,
                         strong_curvature = TRUE,
                         verbose = FALSE) {
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
    res <- cvsrch(phi, step0, alpha = alpha, c1 = c1, c2 = c2,
                  maxfev = maxfev,
                  armijo_check_fn = armijo_check_fn,
                  wolfe_ok_step_fn = wolfe_ok_step_fn, verbose = verbose)
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
#	  ---
#         Calculate the function at x and
#         return this value in the variable f.
#         Calculate the gradient at x and
#         return this vector in g.
#	  ---
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
                   maxfev = Inf, delta = 0.66,
                   armijo_check_fn = armijo_ok_step,
                   wolfe_ok_step_fn = strong_wolfe_ok_step,
                   verbose = FALSE) {

  # increase width by this amount during zoom phase
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
  if (maxfev < 0) {
    params_ok <- FALSE
    problems <- c(problems, paste0("maxfev < 0: ", formatC(maxfev)))
  }
  if (!params_ok) {
    problems <- paste(problems, collapse = "; ")
    stop("Parameter errors detected: ", problems)
  }

  if (maxfev == 0) {
    return(list(step = step0, nfn = 0, info = 3))
  }

  # Check that pv is a descent direction: if not, return a zero step.
  if (step0$d >= 0.0) {
    return(list(step = step0, info = 6, nfn = 0))
  }
  dgtest <- c1 * step0$d

  # Initialize local variables.
  bracketed <- FALSE
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

    if (verbose) {
    message("Bracket: [", formatC(stmin), ", ", formatC(stmax),
            "] alpha = ", formatC(step$alpha))
    }

    # Evaluate the function and gradient at alpha
    # and compute the directional derivative.
    # Additional check: bisect (if needed) until a finite value is found
    # (most important for first iteration)
    ffres <- find_finite(phi, step$alpha, maxfev - nfev, min_alpha = stmin)
    nfev <- nfev + ffres$nfn
    if (!ffres$ok) {
      if (verbose) {
        message("Unable to create finite alpha")
      }

      return(list(step = step0, nfn = nfev, info = 7))
    }
    step <- ffres$step

    # Test for convergence.
    info <- check_convergence(step0, step, brackt, infoc, stmin, stmax,
                              alpha_min, alpha_max, c1, c2, nfev,
                              maxfev, xtol,
                              armijo_check_fn = armijo_check_fn,
                              wolfe_ok_step_fn = wolfe_ok_step_fn,
                              verbose = verbose)

    # Check for termination.
    if (info != 0) {
      # If an unusual termination is to occur, then set step to the best step
      # found
      if (info == 2 || info == 3 || info == 6) {
        step <- stepx
      }
      if (verbose) {
        message("alpha = ", formatC(step$alpha))
      }
      return(list(step = step, info = info, nfn = nfev))
    }

    # In the first stage we seek a step for which the modified
    # function has a nonpositive value and nonnegative derivative.

    # In the original MINPACK the following test is:
    # if (stage1 .and. f .le. ftest1 .and.
    #    *       dg .ge. min(ftol,gtol)*dginit) stage1
    # which translates to: step$f <= f0 + alpha * c1 * d0 &&
    #            step$df >= min(c1, c2) * alpha * d0
    # The second test is the armijo condition and the third is the
    # curvature condition but using the smaller of c1 and
    # c2. This is nearly the standard Wolfe conditions, but because c1 is
    # always <= c2 for a convergent line search, this means
    # we would always use c1 for the curvature condition.
    # I have translated this faithfully, but it seems odd. Using c2 has no
    # effect on the test function from the More'-Thuente paper
    if (stage1 && wolfe_ok_step(step0, step, c1, min(c1, c2))) {
      stage1 <- FALSE
    }

    # A modified function is used to predict the step only if
    # we have not obtained a step for which the modified
    # function has a nonpositive function value and nonnegative
    # derivative, and if a lower function value has been
    # obtained but the decrease is not sufficient.
    if (stage1 && step$f <= stepx$f && !armijo_check_fn(step0, step, c1)) {
      # Define the modified function and derivative values.
      stepxm <- modify_step(stepx, dgtest)
      stepym <- modify_step(stepy, dgtest)
      stepm <- modify_step(step, dgtest)

      step_result <- cstep(stepxm, stepym, stepm, brackt, stmin, stmax,
                           verbose = verbose)

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
      step_result <- cstep(stepx, stepy, step, brackt, stmin, stmax,
                           verbose = verbose)
      brackt <- step_result$brackt
      infoc <- step_result$info
      stepx <- step_result$stepx
      stepy <- step_result$stepy
      step <- step_result$step
    }

    if (!bracketed && brackt) {
      bracketed <- TRUE
      if (verbose) {
        message("Bracketed")
      }
    }

    # Force a sufficient decrease in the size of the interval of uncertainty.
    if (brackt) {
      # if the length of I does not decrease by a factor of delta < 1
      # then use a bisection step for the next trial alpha
      width_new <- abs(stepy$alpha - stepx$alpha)
      if (width_new >= delta * width_old) {
        if (verbose) {
          message("Interval did not decrease sufficiently: bisecting")
        }
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
# NB dgtest was originally used in testing for min/max alpha test (code 4 and 5)
# but has been replaced with a call to the curvature test using c1 instead of c2
# so dgtest is no longer used in the body of the function.
check_convergence <- function(step0, step, brackt, infoc, stmin, stmax,
                              alpha_min, alpha_max, c1, c2, nfev,
                              maxfev, xtol, armijo_check_fn = armijo_ok_step,
                              wolfe_ok_step_fn = strong_wolfe_ok_step,
                              verbose = FALSE) {
  info <- 0
  if ((brackt && (step$alpha <= stmin || step$alpha >= stmax)) || infoc == 0) {
    if (verbose) {
      message("MT: Rounding errors prevent further progress: stmin = ",
            formatC(stmin), " stmax = ", formatC(stmax))
    }
    # rounding errors prevent further progress
    info <- 6
  }
  # use of c1 in curvature check is on purpose (it's in the MINPACK code)
  if (step$alpha == alpha_max && armijo_check_fn(step0, step, c1) &&
      !curvature_ok_step(step0, step, c1)) {
    # reached alpha_max
    info <- 5
    if (verbose) {
      message("MT: Reached alpha max")
    }

  }
  # use of c1 in curvature check here is also in MINPACK code
  if (step$alpha == alpha_min && (!armijo_check_fn(step0, step, c1) ||
                                  curvature_ok_step(step0, step, c1))) {
    # reached alpha_min
    info <- 4
    if (verbose) {
      message("MT: Reached alpha min")
    }

  }
  if (nfev >= maxfev) {
    # maximum number of function evaluations reached
    info <- 3
    if (verbose) {
      message("MT: exceeded fev")
    }
  }
  if (brackt && stmax - stmin <= xtol * stmax) {
    # interval width is below xtol
    info <- 2
    if (verbose) {
      message("MT: interval width is <= xtol: ", formatC(xtol * stmax))
    }

  }
  if (wolfe_ok_step_fn(step0, step, c1, c2)) {
    # success!
    info <- 1
    if (verbose) {
      message("Success!")
    }
  }
  info
}


# Part of the More'-Thuente line search.
#
# Updates the interval of uncertainty of the current step size and updates the
# current best step size.
#
# This routine is a translation of Dianne O'Leary's Matlab code, which was
# itself a translation of the MINPACK original. Original comments to the Matlab
# code are at the end.
#
# @param stepx One side of the updated step interval, and the associated
#     line search values.
# @param stepy Other side of the updated step interval, and the
#     associated line search values.
# @param step Optimal step size and associated line search
#     value.
# @param brackt TRUE if the interval has been bracketed.
# @param stpmin Minimum allowed interval length.
# @param stpmax Maximum allowed interval length.
# @return List containing:
# \itemize{
#   \item \code{stepx} One side of the updated step interval, and the associated
#     line search values.
#   \item \code{stepy} Other side of the updated step interval, and the
#     associated line search values.
#   \item \code{step} Updated optimal step size and associated line search
#     value.
#   \item \code{brackt} TRUE if the interval has been bracketed.
#   \item \code{info} Integer return code.
# }
# The possible integer return codes refer to the cases 1-4 enumerated in the
# original More'-Thuente paper that correspond to different line search values
# at the ends of the interval and the current best step size (and therefore
# the type of cubic or quadratic interpolation). An integer value of 0 indicates
# that the input parameters are invalid.
#
#
#   Translation of minpack subroutine cstep
#   Dianne O'Leary   July 1991
#     **********
#
#     Subroutine cstep
#
#     The purpose of cstep is to compute a safeguarded step for
#     a linesearch and to update an interval of uncertainty for
#     a minimizer of the function.
#
#     The parameter stx contains the step with the least function
#     value. The parameter stp contains the current step. It is
#     assumed that the derivative at stx is negative in the
#     direction of the step. If brackt is set true then a
#     minimizer has been bracketed in an interval of uncertainty
#     with end points stx and sty.
#
#     The subroutine statement is
#
#       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
#                        stpmin,stpmax,info)
#
#     where
#
#       stx, fx, and dx are variables which specify the step,
#         the function, and the derivative at the best step obtained
#         so far. The derivative must be negative in the direction
#         of the step, that is, dx and stp-stx must have opposite
#         signs. On output these parameters are updated appropriately.
#
#       sty, fy, and dy are variables which specify the step,
#         the function, and the derivative at the other endpoint of
#         the interval of uncertainty. On output these parameters are
#         updated appropriately.
#
#       stp, fp, and dp are variables which specify the step,
#         the function, and the derivative at the current step.
#         If brackt is set true then on input stp must be
#         between stx and sty. On output stp is set to the new step.
#
#       brackt is a logical variable which specifies if a minimizer
#         has been bracketed. If the minimizer has not been bracketed
#         then on input brackt must be set false. If the minimizer
#         is bracketed then on output brackt is set true.
#
#       stpmin and stpmax are input variables which specify lower
#         and upper bounds for the step.
#
#       info is an integer output variable set as follows:
#         If info <- 1,2,3,4,5, then the step has been computed
#         according to one of the five cases below. Otherwise
#         info <- 0, and this indicates improper input parameters.
#
#     Subprograms called
#
#       FORTRAN-supplied ... abs,max,min,sqrt
#                        ... dble
#
#     Argonne National Laboratory. MINPACK Project. June 1983
#     Jorge J. More', David J. Thuente
#
#     **********
cstep <-  function(stepx, stepy, step, brackt, stpmin, stpmax,
                   verbose = FALSE) {

  stx <- stepx$alpha
  fx <- stepx$f
  dx <- stepx$d
  dfx <- stepx$df

  sty <- stepy$alpha
  fy <- stepy$f
  dy <- stepy$d
  dfy <- stepy$df

  stp <- step$alpha
  fp <- step$f
  dp <- step$d
  dfp <- step$df

  delta <- 0.66
  info <- 0
  # Check the input parameters for errors.
  if ((brackt && (stp <= min(stx, sty) ||
                  stp >= max(stx, sty))) ||
      dx * (stp - stx) >= 0.0 || stpmax < stpmin) {
    list(
      stepx = stepx,
      stepy = stepy,
      step = step,
      brackt = brackt, info = info)
  }
  # Determine if the derivatives have opposite sign.
  sgnd <- dp * (dx / abs(dx))

  # First case. Trial function value is larger, so choose a step which is
  # closer to stx.
  # The minimum is bracketed.
  # If the cubic step is closer to stx than the quadratic step, the cubic step
  # is taken else the average of the cubic and quadratic steps is taken.
  if (fp > fx) {
    info <- 1
    bound <- TRUE

    stpc <- cubic_interpolate(stx, fx, dx, stp, fp, dp, ignoreWarnings = TRUE)
    stpq <- quadratic_interpolate(stx, fx, dx, stp, fp)

    if (is.nan(stpc)) {
      stpf <- stpq
    }
    else {
      if (abs(stpc - stx) < abs(stpq - stx)) {
        stpf <- stpc
      } else {
        stpf <- stpc + (stpq - stpc) / 2
      }
    }
    brackt <- TRUE

    # Second case. A lower function value and derivatives of
    # opposite sign. The minimum is bracketed. If the cubic
    # step is closer to stx than the quadratic (secant) step,
    # the cubic step is taken, else the quadratic step is taken.
  } else if (sgnd < 0.0) {
    info <- 2
    bound <- FALSE

    stpc <- cubic_interpolate(stx, fx, dx, stp, fp, dp, ignoreWarnings = TRUE)
    stpq <- quadratic_interpolateg(stp, dp, stx, dx)
    if (is.nan(stpc)) {
      stpf <- stpq
    }
    else {
      if (abs(stpc - stp) > abs(stpq - stp)) {
        stpf <- stpc
      } else {
        stpf <- stpq
      }
    }

    brackt <- TRUE

    # Third case. A lower function value, derivatives of the
    # same sign, and the magnitude of the derivative decreases.
    # The next trial step exists outside the interval so is an extrapolation.
    # The cubic may not have a minimizer. If it does, it may be in the
    # wrong direction, e.g. stc < stx < stp
    # The cubic step is only used if the cubic tends to infinity
    # in the direction of the step and if the minimum of the cubic
    # is beyond stp. Otherwise the cubic step is defined to be
    # either stpmin or stpmax. The quadratic (secant) step is also
    # computed and if the minimum is bracketed then the the step
    # closest to stx is taken, else the step farthest away is taken.
  } else if (abs(dp) < abs(dx)) {
    info <- 3
    bound <- TRUE
    theta <- 3 * (fx - fp) / (stp - stx) + dx + dp
    s <- norm(rbind(theta, dx, dp), "i")
    # The case gamma = 0 only arises if the cubic does not tend
    # to infinity in the direction of the step.
    gamma <- s * sqrt(max(0.,(theta / s) ^ 2 - (dx / s) * (dp / s)))
    if (stp > stx) {
      gamma <- -gamma
    }
    p <- (gamma - dp) + theta
    q <- (gamma + (dx - dp)) + gamma
    r <- p / q

    if (r < 0.0 && gamma != 0.0) {
      stpc <- stp + r * (stx - stp)
    } else if (stp > stx) {
      stpc <- stpmax
    } else {
      stpc <- stpmin
    }

    stpq <- quadratic_interpolateg(stp, dp, stx, dx)

    if (brackt) {
      if (abs(stp - stpc) < abs(stp - stpq)) {
        stpf <- stpc
      } else {
        stpf <- stpq
      }
    } else {
      if (abs(stp - stpc) > abs(stp - stpq)) {
        stpf <- stpc
      } else {
        stpf <- stpq
      }
    }
    # Fourth case. A lower function value, derivatives of the
    # same sign, and the magnitude of the derivative does
    # not decrease. If the minimum is not bracketed, the step
    # is either stpmin or stpmax, else the cubic step is taken.
  } else {
    info <- 4
    bound <- FALSE
    if (brackt) {
      stpc <- cubic_interpolate(sty, fy, dy, stp, fp, dp, ignoreWarnings = TRUE)
      if (is.nan(stpc)) {
        stpc <- (sty + stp) / 2
      }
      stpf <- stpc
    } else if (stp > stx) {
      stpf <- stpmax
    } else {
      stpf <- stpmin
    }
  }

  # Update the interval of uncertainty. This update does not
  # depend on the new step or the case analysis above.
  if (fp > fx) {
    sty <- stp
    fy <- fp
    dy <- dp
    dfy <- dfp
  } else {
    if (sgnd < 0.0) {
      sty <- stx
      fy <- fx
      dy <- dx
      dfy <- dfx
    }
    stx <- stp
    fx <- fp
    dx <- dp
    dfx <- dfp
  }

  # Compute the new step and safeguard it.
  stpf <- min(stpmax, stpf)
  stpf <- max(stpmin, stpf)
  stp <- stpf
  if (brackt && bound) {
    # if the new step is too close to an end point
    # replace with a (weighted) bisection (delta = 0.66 in the paper)
    if (verbose) {
      message("Step too close to end point, weighted bisection")
    }
    stb <- stx + delta * (sty - stx)
    if (sty > stx) {
      stp <- min(stb, stp)
    } else {
      stp <- max(stb, stp)
    }
  }
  list(
    stepx = list(alpha = stx, f = fx, d = dx, df = dfx),
    stepy = list(alpha = sty, f = fy, d = dy, df = dfy),
    step = list(alpha = stp, f = fp, d = dp, df = dfp),
    brackt = brackt, info = info)
}
