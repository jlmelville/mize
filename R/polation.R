# Interpolation and extrapolation functions.

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
# @param suppress_warnings If TRUE, don't warn if the interpolation creates a
#   non-finite value.
# @return Cubic interpolated estimate of minimum value of x.
cubic_interpolate <- function(
  x1,
  f1,
  g1,
  x2,
  f2,
  g2,
  suppress_warnings = FALSE
) {
  cubic_shape <- g1 + g2 - 3 * ((f1 - f2) / (x1 - x2))

  if (suppress_warnings) {
    suppressWarnings(
      cubic_root <- sign(x2 - x1) * sqrt(cubic_shape^2 - g1 * g2)
    )
  } else {
    cubic_root <- sign(x2 - x1) * sqrt(cubic_shape^2 - g1 * g2)
  }
  x2 -
    (x2 - x1) *
      ((g2 + cubic_root - cubic_shape) / (g2 - g1 + 2 * cubic_root))
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
  x1 - (0.5 * g1 * (x2 - x1)^2) / (f2 - f1 - g1 * (x2 - x1))
}

propose_quadratic_alpha <- function(first_step, second_step) {
  quadratic_interpolate(
    first_step$alpha,
    first_step$f,
    first_step$d,
    second_step$alpha,
    second_step$f
  )
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
