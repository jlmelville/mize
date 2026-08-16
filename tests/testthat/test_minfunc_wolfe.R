
# fmt: skip file

# Reference values for the supported Schmidt gradient/cubic Wolfe profile.

mfls <- function(fg, x, alpha, c1, c2, pv = -fg$gr(x) / abs(fg$gr(x)),
                 eps = 1e-6, approx_armijo = FALSE,
                 strong_curvature = TRUE) {
  step0 <- make_step0(fg, x, pv)

  res <- new_wolfe_line_search(
    schmidt_core,
    armijo_constant = c1,
    curvature_constant = c2,
    max_evaluations = 10000,
    approximation_tolerance = eps,
    approximate_armijo = approx_armijo,
    strong_curvature = strong_curvature,
    method_policy = new_schmidt_bracket_zoom_policy()
  )(
    phi = make_phi_alpha(x, fg, pv, calc_gradient_default = TRUE),
    step0 = step0,
    alpha = alpha,
    pm = pv
  )

  res$step$par <- x + res$step$alpha * pv
  res
}

## Test data from the More'-Thuente paper.

# Test Grad/Cubic interpolation (default)

# Table 1
test_that("Table 1 Grad/Cubic Interpolation", {
  res11 <- mfls(fg = f1, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.1)
  expect_step(res11, x = 1.4218, f = -0.35355, df = 0.0013375, nfev = 6)
  res12 <- mfls(fg = f1, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.1)
  expect_step(res12, x = 1.4218, f = -0.35355, df = 0.0013220, nfev = 4)
  res13 <- mfls(fg = f1, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.1)
  expect_step(res13, x = 10, f = -0.098039, df = 0.0094195, nfev = 1)
  # progtol
  res14 <- mfls(fg = f1, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.1)
  expect_step(res14, x = 333.34, f = -0.0029999, df = 8.9994e-006, nfev = 37)
})

# Table 2
test_that("Table 2 Grad/Cubic Interpolation", {
  res21 <- mfls(fg = f2, x = 0, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(res21, x = 1.5960, f = -2.6214, df = 2.1828e-014, nfev = 13)
  res22 <- mfls(fg = f2, x = 0, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(res22, x = 1.5960, f = -2.6214, df = -9.4587e-013, nfev = 10)
  res23 <- mfls(fg = f2, x = 0, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(res23, x = 1.5960, f = -2.6214, df = 7.2760e-014, nfev = 9)
  res24 <- mfls(fg = f2, x = 0, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(res24, x = 1.5960, f = -2.6214, df = 6.8758e-012, nfev = 14)
})

# Table 3
test_that("Table 3 Grad/Cubic Interpolation", {
  res31 <- mfls(fg = f3, x = 0, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(res31, x = 1.0, f = -0.011160, df = 8.6906e-006, nfev = 14)
  res32 <- mfls(fg = f3, x = 0, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(res32, x = 1.0, f = -0.011160, df = -1.9129e-005, nfev = 11)
  res33 <- mfls(fg = f3, x = 0, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(res33, x = 1.0, f = -0.011160, df = -2.3751e-006, nfev = 9)
  res34 <- mfls(fg = f3, x = 0, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(res34, x = 1.0, f = -0.011160, df = -4.1259e-006, nfev = 12)
})

# Table 4
test_that("Table 4 Grad/Cubic Interpolation", {
  res41 <- mfls(fg = f4, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(res41, x = 0.030526, f = 0.99902, df = -5.3509e-004, nfev = 4)
  res42 <- mfls(fg = f4, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(res42, x = 0.1, f = 0.99901, df = -4.9330e-005, nfev = 1)
  res43 <- mfls(fg = f4, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(res43, x = 0.34963, f = 0.99900, df = -2.9053e-006, nfev = 3)
  res44 <- mfls(fg = f4, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(res44, x = 0.90475, f = 0.99901, df = 5.4438e-005, nfev = 4)
})

test_that("Table 5 Grad/Cubic Interpolation", {
  res51 <- mfls(fg = f5, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(res51, x = 0.074271, f = 0.99138, df = 1.7174e-005, nfev = 7)
  res52 <- mfls(fg = f5, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(res52, x = 0.071748, f = 0.99139, df = -6.1309e-004, nfev = 3)
  res53 <- mfls(fg = f5, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(res53, x = 0.075308, f = 0.99138, df = 2.5835e-004, nfev = 6)
  res54 <- mfls(fg = f5, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(res54, x = 0.073117, f = 0.99138, df = -2.6323e-004, nfev = 8)
})

test_that("Table 6 Grad/Cubic Interpolation", {
  res61 <- mfls(fg = f6, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(res61, x = 0.92579, f = 0.99138, df = -3.5017e-006, nfev = 8)
  res62 <- mfls(fg = f6, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(res62, x = 0.92706, f = 0.99139, df = 3.0799e-004, nfev = 11)
  res63 <- mfls(fg = f6, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(res63, x = 0.92593, f = 0.99138, df = 3.1671e-005, nfev = 7)
  res64 <- mfls(fg = f6, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(res64, x = 0.92918, f = 0.99139, df = 8.6261e-004, nfev = 10)
})

test_that("MT Function modification Grad/Cubic Interpolation", {
  # progtol
  res4m <- mfls(fg = f4, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(res4m, x = 0.5, f = 0.999, df = 1.3132e-011, alpha = 0.5, nfev = 12)
  # progtol
  res5m <- mfls(fg = f5, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(res5m, x = 0.50112, f = 0.99464, df = 0.0087536, alpha = 0.49888, nfev = 14)
  # progtol
  res6m <- mfls(fg = f6, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(res6m, x = 0.82848, f = 0.99188, df = -0.0072578, alpha = 0.17152, nfev = 14)
})
