# Interpolation and extrapolation functions.

# Estimate Minimum By Cubic Extrapolation
#
# Carries out cubic extrapolation based on the x, f(x), and f'(x) values
# at two points to find minimum value of x.
#
# @param x1 x value at first point.
# @param f1 f(x) value at first point.
# @param g1 f'(x) value at first point.
# @param x2 x value at second point.
# @param f2 f(x) value at second point.
# @param g2 f'(x) value at second point.
# @param ignoreWarnings If TRUE, don't warn if the extrapolation creates a
#   non-finite value.
# @return Cubic extrapolated estimate of minimum value of x.
cubic_extrapolate <- function(x1, f1, g1, x2, f2, g2, ignoreWarnings = FALSE) {
  A <- 6 * (f1 - f2) + 3 * (g2 + g1) * (x2 - x1)
  B <- 3 * (f2 - f1) - (2 * g1 + g2) * (x2 - x1)
  if (ignoreWarnings) {
    suppressWarnings(
      x1 - g1 * (x2 - x1) ^ 2 / (B + sqrt(B * B - A * g1 * (x2 - x1)))
    )
  }
}

# Estimate Step Size Minimum By Cubic Extrapolation
#
# Estimates step size corresponding to minimum of the line function using
# cubic extrapolation from two line search evaluations with both function
# and directional derivatives calculated.
#
# @param step1 Line search information for first step value.
# @param step2 Line search information for second step value.
# @return Cubic extrapolated estimate of step size which minimizes the line
#   function.
cubic_extrapolate_step <- function(step1, step2) {
  cubic_extrapolate(step1$alpha, step1$f, step1$d, step2$alpha, step2$f,
                    step2$d, ignoreWarnings = TRUE)
}

# Estimate Minimum By Cubic Interpolation
#
# Carries out cubic interpolation based on the x, f(x), and f'(x) values
# at two points to find minimum value of x.
#
# @param x1 x value at first point.
# @param f1 f(x) value at first point.
# @param g1 f'(x) value at first point.
# @param x2 x value at second point.
# @param f2 f(x) value at second point.
# @param g2 f'(x) value at second point.
# @param ignoreWarnings If TRUE, don't warn if the interpolation creates a
#   non-finite value.
# @return Cubic interpolated estimate of minimum value of x.
cubic_interpolate <- function(x1, f1, g1, x2, f2, g2, ignoreWarnings = FALSE) {
  # nwc(x1, f1, g1, x2, f2, g2)
  # A <- 6 * (f1 - f2) / (x2 - x1) + 3 * (g2 + g1)
  # B <- 3 * (f2 - f1) - (2 * g1 + g2) * (x2 - x1)
  # # num. error possible, ok!
  # suppressWarnings(
  #   x1 + (sqrt(B * B - A * g1 * (x2 - x1) ^ 2) -  B) / A
  # )
#  A <- 6 * (f1 - f2) + 3 * (g2 + g1) * (x2 - x1)
#  B <- 3 * (f2 - f1) - (2 * g1 + g2) * (x2 - x1)
#  x1 - g1 * (x2 - x1) ^ 2 / (B + sqrt(B * B - A * g1 * (x2 - x1)))
  d1 <- g1 + g2 - 3 * ((f1 - f2) / (x1 - x2))

  if (ignoreWarnings) {
    suppressWarnings(
      d2 <- sign(x2 - x1) * sqrt(d1 * d1 - g1 * g2)
    )
  }
  else {
    d2 <- sign(x2 - x1) * sqrt(d1 * d1 - g1 * g2)
  }
  x2 - (x2 - x1) * ((g2 + d2 - d1) / (g2 - g1 + 2 * d2))
}

# Estimate Step Size Minimum By Cubic Interpolation
#
# Estimates step size corresponding to minimum of the line function using
# cubic interpolation from two line search evaluations with both function
# and directional derivatives calculated.
#
# @param step1 Line search information for first step value.
# @param step2 Line search information for second step value.
# @return Cubic interpolated estimate of step size which minimizes the line
#   function.
cubic_interpolate_step <- function(step1, step2) {
  cubic_interpolate(step1$alpha, step1$f, step1$d,
                    step2$alpha, step2$f, step2$d)
}

# Estimate Minimum By Quadratic Interpolation With One Gradient
#
# Carries out quadratic interpolation based on the x and f(x) values at two
# points, and the f'(x) value at the first point, to find minimum value of x.
#
# @param x1 x value at first point.
# @param f1 f(x) value at first point.
# @param g1 f'(x) value at first point.
# @param x2 x value at second point.
# @param f2 f(x) value at second point.
# @return Quadratic interpolated estimate of minimum value of x.
quadratic_interpolate <- function(x1, f1, g1, x2, f2) {
  x1 - (0.5 * g1 * (x2 - x1) ^ 2) / (f2 - f1 - g1 * (x2 - x1))
}


# Estimate Step Size Minimum By Quadratic Interpolation
#
# Estimates step size corresponding to minimum of the line function using
# quadratic interpolation from two line search evaluations. The function must
# have been evaluated at both points, but only the directional derivative at
# the first point is used.
#
# @param step1 Line search information for first step value.
# @param step2 Line search information for second step value.
# @return Quadratic interpolated estimate of step size which minimizes the line
#   function.
quadratic_interpolate_step <- function(step1, step2) {
  quadratic_interpolate(step1$alpha, step1$f, step1$d,
                        step2$alpha, step2$f)
}

# Estimate Minimum By Quadratic Interpolation With Two Gradients
#
# Carries out quadratic interpolation based on the x and f'(x) values at two
# points. Note that this does not use the function values at either of the
# points.
#
# @param x1 x value at first point.
# @param g1 f'(x) value at first point.
# @param x2 x value at second point.
# @param g2 f'(x) value at second point.
# @return Quadratic interpolated estimate of minimum value of x.
quadratic_interpolateg <- function(x1, g1, x2, g2) {
  x2 + (x1 - x2) * g2 / (g2 - g1)
}

# Tweak Extrapolated Point
#
# Prevents the extrapolated point from being too far away from or to close to
# the points used in the extrapolation.
#
# @param xnew 1D position of the new point.
# @param x1 1D position of the first points used in the extrapolation.
# @param x2 1D position of the second point used in the extrapolation.
# @param ext Maximum multiple of \code{x2} that \code{xnew} is allowed to be
#  extrapolated to.
# @param int Given the distance between \code{x1} and \code{x2}, specified what
#  multiple of that distance is the minimum allowed distance for \code{xnew}
#  from \code{x2}.
# @return A value of \code{xnew} that obeys the minimum and maximum distance
#  constraints from \code{x2}.
tweak_extrapolation <- function(xnew, x1, x2, ext, int) {
  # num prob | wrong sign?
  if (!is.double(xnew) || is.nan(xnew) || is.infinite(xnew) || xnew < 0) {
    # extrapolate maximum amount
    xnew <- x2 * ext
  } else if (xnew > (x2 * ext)) {
    # new point beyond extrapolation limit?
    # extrapolate maximum amount
    xnew <- x2 * ext
  } else if (xnew < (x2 + int * (x2 - x1))) {
    # new point too close to previous point?
    xnew <- x2 + int * (x2 - x1)
  }
  xnew
}

# Tweak Interpolated Point
#
# Prevents interpolated point from getting too close to either of the
# points used for the interpolation. If the point is not a number or infinite,
# then it is set to the bisection of the position of the two interpolating
# points before the check for a too-close approach is carried out.
#
# @param xnew Position of the interpolated point.
# @param x1 Position of the first point used for interpolation.
# @param x2 Position of the second point used for interpolation.
# @param int Given the distance between \code{x1} and \code{x2}, specifies what
#  multiple of that distance is the minimum allowed distance for \code{xnew}
#  from \code{x1} or \code{x2}.
# @return Tweaked position of \code{xnew} such that it is not too close to
#  either \code{x1} or \code{x2}.
tweak_interpolation <- function(xnew, x1, x2, int) {
  if (is.nan(xnew) || is.infinite(xnew)) {
    # if we had a numerical problem then bisect
    xnew <- (x1 + x2) / 2
  }
  # don't accept too close
  max(min(xnew, x2 - int * (x2 - x1)), x1 + int * (x2 - x1))
}
