mixedExtrap <- function(x0, f0, g0, x1, f1, g1, minStep, maxStep, debug) {
  alpha_c <- polyinterp(point_matrix(c(x0, x1), c(f0, f1), c(g0, g1)),
                        minStep, maxStep, debug = TRUE)
  alpha_s <- polyinterp(point_matrix(c(x0, x1), c(f0, NA), c(g0, g1)),
                        minStep, maxStep, debug = TRUE)
  if (debug) {
    message("cubic ext = ", formatC(alpha_c), " secant ext = ", formatC(alpha_s),
            " minStep = ", formatC(minStep),
            " alpha_c > minStep ? ", alpha_c > minStep,
            " |ac - x1| = ", formatC(abs(alpha_c - x1)),
            " |as - x1| = ", formatC( abs(alpha_s - x1))
            )
  }
  if (alpha_c > minStep && abs(alpha_c - x1) < abs(alpha_s - x1)) {
    if (debug) {
      message('Cubic Extrapolation ', formatC(alpha_c))
    }
    res <- alpha_c
  }
  else {
    if (debug) {
      message('Secant Extrapolation ', formatC(alpha_s))
    }
    res <- alpha_s
  }
  res
}
