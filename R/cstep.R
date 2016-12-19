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
cstep <-  function(stepx, stepy, step, brackt, stpmin, stpmax) {

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
        stpf <- (sty + stp) / 2
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
