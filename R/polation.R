# Interpolation and extrapolation functions.

# Estimate Minimum By Cubic Interpolation
#
# Carries out cubic interpolation based on the alpha, objective, and slope
# values at two line points to estimate the minimizing alpha.
#
# @param first_alpha Alpha at the first point.
# @param first_value Objective value at the first point.
# @param first_slope Directional derivative at the first point.
# @param second_alpha Alpha at the second point.
# @param second_value Objective value at the second point.
# @param second_slope Directional derivative at the second point.
# @return Cubic interpolated estimate of the minimizing alpha, or `NA_real_` when
#   no finite proposal can be computed.
cubic_interpolate <- function(
  first_alpha,
  first_value,
  first_slope,
  second_alpha,
  second_value,
  second_slope
) {
  alpha_difference <- first_alpha - second_alpha
  if (!isTRUE(is.finite(alpha_difference)) || alpha_difference == 0) {
    return(NA_real_)
  }

  cubic_shape <- first_slope +
    second_slope -
    3 * ((first_value - second_value) / alpha_difference)
  discriminant <- cubic_shape^2 - first_slope * second_slope
  if (!isTRUE(is.finite(discriminant)) || discriminant < 0) {
    return(NA_real_)
  }

  discriminant_root <- sign(second_alpha - first_alpha) * sqrt(discriminant)
  denominator <- second_slope - first_slope + 2 * discriminant_root
  if (!isTRUE(is.finite(denominator)) || denominator == 0) {
    return(NA_real_)
  }

  proposed_alpha <- second_alpha -
    (second_alpha - first_alpha) *
      ((second_slope + discriminant_root - cubic_shape) / denominator)
  if (!isTRUE(is.finite(proposed_alpha))) {
    return(NA_real_)
  }
  proposed_alpha
}

# Estimate Minimum By Quadratic Interpolation With One Gradient
#
# Carries out quadratic interpolation based on the alpha and objective values
# at two points, and the slope at the first point.
#
# @param first_alpha Alpha at the first point.
# @param first_value Objective value at the first point.
# @param first_slope Directional derivative at the first point.
# @param second_alpha Alpha at the second point.
# @param second_value Objective value at the second point.
# @return Quadratic interpolated estimate of the minimizing alpha.
quadratic_interpolate <- function(
  first_alpha,
  first_value,
  first_slope,
  second_alpha,
  second_value
) {
  first_alpha -
    (0.5 * first_slope * (second_alpha - first_alpha)^2) /
      (second_value - first_value - first_slope * (second_alpha - first_alpha))
}

propose_quadratic_alpha <- function(first_point, second_point) {
  quadratic_interpolate(
    first_point$alpha,
    first_point$value,
    first_point$slope,
    second_point$alpha,
    second_point$value
  )
}

# Estimate Minimum By Slope Secant
#
# Finds where the secant through two directional derivatives reaches zero.
#
# @param first_alpha Alpha at the first point.
# @param first_slope Directional derivative at the first point.
# @param second_alpha Alpha at the second point.
# @param second_slope Directional derivative at the second point.
# @return Secant estimate of the alpha with zero directional derivative.
propose_slope_secant_alpha <- function(
  first_alpha,
  first_slope,
  second_alpha,
  second_slope
) {
  second_alpha +
    (first_alpha - second_alpha) * second_slope / (second_slope - first_slope)
}
