# function [minPos] <- polyinterp(points,doPlot,xminBound,xmaxBound)
#
#   Minimum of interpolating polynomial based on function and derivative
#   values
#
#   It can also be used for extrapolation if {xmin,xmax} are outside
#   the domain of the points.
#
#   Input:
#       points(pointNum,[x f g])
#       xmin: min value that brackets minimum (default: min of points)
#       xmax: max value that brackets maximum (default: max of points)
#
#   set f or g to sqrt(-1) if they are not known
#   the order of the polynomial is the number of known f and g values minus 1
# points position, function and gradient values to interpolate.
# An n x 3 matrix where n is the number of points and each row contains
# x, f, g in columns 1-3 respectively.
# @return minPos
polyinterp <- function(points,
                       xminBound = range(points[, 1])[1],
                       xmaxBound = range(points[, 1])[2],
                       debug = FALSE) {

  # the number of known f and g values minus 1
  order <- sum(!is.na(points[, 2:3])) - 1

  # Code for most common case:
  #   - cubic interpolation of 2 points
  #       w/ function and derivative values for both
  if (nrow(points) == 2 && order == 3) {
    if (debug) {
      message("polyinterp common case")
    }
    # Solution in this case (where x2 is the farthest point):
    #    d1 <- g1 + g2 - 3*(f1-f2)/(x1-x2);
    #    d2 <- sqrt(d1^2 - g1*g2);
    #    minPos <- x2 - (x2 - x1)*((g2 + d2 - d1)/(g2 - g1 + 2*d2));
    #    t_new <- min(max(minPos,x1),x2);
    minPos <- which.min(points[, 1])
    notMinPos <- -minPos + 3

    x1 <- points[minPos, 1]
    x2 <- points[notMinPos, 1]
    f1 <- points[minPos, 2]
    f2 <- points[notMinPos, 2]
    g1 <- points[minPos, 3]
    g2 <- points[notMinPos, 3]

    d1 <- g1 + g2 - 3 * (f1 - f2) / (x1 - x2)
    # d1 <- points[minPos, 3] + points[notMinPos, 3] - 3 *
    #   (points[minPos, 2] - points[notMinPos, 2]) /
    #   (points[minPos, 1] - points[notMinPos, 1])

    #d2sq <- d1 ^ 2 - points[minPos, 3] * points[notMinPos, 3]
    d2sq <- d1 ^ 2 - g1 * g2

    if (d2sq >= 0) {
      d2 <- sqrt(d2sq)

      # x <- points[notMinPos, 1] -
      #   (points[notMinPos, 1] - points[minPos, 1]) *
      #   ((points[notMinPos, 3] + d2 - d1) /
      #      (points[notMinPos, 3] - points[minPos, 3] + 2 * d2))
      x <- x2 - (x2 - x1) * ((g2 + d2 - d1) / (g2 - g1 + 2 * d2))
      if (debug) { message("d2 is real ", formatC(d2), " x = ", formatC(x)) }

      minPos <- min(max(x, xminBound), xmaxBound)
    }
    else {
      if (debug) { message("d2 is not real, bisecting") }

      minPos <- (xmaxBound + xminBound) / 2
    }
    return(minPos)
  }

  params <- polyfit(points)
  # If polynomial couldn't be found (due to singular matrix), bisect
  if (is.null(params)) {
    return((xminBound + xmaxBound) / 2)
  }

  # Compute Critical Points
  dParams <- rep(0, order)
  for (i in 1:order) {
     dParams[i] <- params[i + 1] * i
  }

  cp <- unique(c(xminBound, points[, 1], xmaxBound))

  # Remove mad, bad and dangerous to know critical points:
  # Must be finite, non-complex and not an extrapolation
  if (all(is.finite(dParams))) {
    cp <- c(cp,
            Re(Filter(
              function(x) {
                abs(Im(x)) < 1e-8 &&
                Re(x) >= xminBound &&
                Re(x) <= xmaxBound },
              polyroot(dParams))))
  }

  # Test Critical Points
  fcp <- polyval(cp, params)
  fminpos <- which.min(fcp)
  if (is.finite(fcp[fminpos])) {
    minpos <- cp[fminpos]
  }
  else {
    # Default to bisection if no critical points valid
    minpos <- (xminBound + xmaxBound) / 2
  }
  minpos
}

# Fits a polynomial to the known function and gradient values. The order of
# the polynomial is the number of known function and gradient values, minus one.
# points - an n x 3 matrix where n is the number of points and each row contains
#          x, f, g in columns 1-3 respectively.
# returns an array containing the coefficients of the polynomial in increasing
# order, e.g. c(1, 2, 3) is the polynomial 1 + 2x + 3x^2
# Returns NULL if the solution is singular
polyfit <- function(points) {
  nPoints <- nrow(points)
  # the number of known f and g values minus 1
  order <- sum(!is.na(points[, 2:3])) - 1

  # Constraints Based on available Function Values
  A <- NULL
  b <- NULL
  for (i in 1:nPoints) {
    if (!is.na(points[i, 2])) {
      constraint <- rep(0, order + 1)
      for (j in rev(0:order)) {
        constraint[order - j + 1] <- points[i, 1] ^ j
      }
      if (is.null(A)) {
        A <- constraint
      }
      else {
        A <- rbind(A, constraint)
      }
      if (is.null(b)) {
        b <- points[i, 2]
      }
      else {
        b <- c(b, points[i, 2])
      }
    }
  }

  # Constraints based on available Derivatives
  for (i in 1:nPoints) {
    if (!is.na(points[i, 3])) {
      constraint <- rep(0, order + 1)
      for (j in 1:order) {
        constraint[j] <- (order - j + 1) * points[i, 1] ^ (order - j)
      }
      if (is.null(A)) {
        A <- constraint
      }
      else {
        A <- rbind(A, constraint)
      }
      if (is.null(b)) {
        b <- points[i, 3]
      }
      else {
        b <- c(b, points[i, 3])
      }
    }
  }
  # Find interpolating polynomial
  params <- try(solve(A, b), silent = TRUE)
  if (class(params) == "numeric") {
    params <- rev(params)
  }
  else {
    params <- NULL
  }
}

# Evaluate 1D polynomial with coefs over the set of points x
# coefs - the coefficients for the terms of the polynomial ordered by
#   increasing degree, i.e. c(1, 2, 3, 4) represents the polynomial
#   4x^3 + 3x^2 + 2x + 1. This is the reverse of the ordering used by the Matlab
#   function, but is consistent with R functions like poly and polyroot
#   Also, the order of the arguments is reversed from the Matlab function
# Returns array of values of the evaluated polynomial
polyval <- function(x, coefs) {
  deg <- length(coefs) - 1
  # Sweep multiplies each column of the poly matrix by the coefficient
  rowSums(sweep(stats::poly(x, degree = deg, raw = TRUE),
                2, coefs[2:length(coefs)], `*`)) + coefs[1]
}

point_matrix <- function(xs, fs, gs) {
  matrix(c(xs, fs, gs), ncol = 3)
}
