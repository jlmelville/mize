context("Schmidt (MinFunc) Line Search")

# Test values produced by running the original minFunc WolfeLineSearch function
# using Octave 4.2.0

mfls <- function(fg, x, pv, alpha, c1, c2, LS_interp = 2, debug = FALSE) {
  step0 <- make_step0(fg, x, pv)

  res <- WolfeLineSearch(alpha = alpha, f = step0$f, g = step0$df,
                         gtd = step0$d,
                         c1 = c1, c2 = c2, LS_interp = LS_interp, LS_multi = 0,
                         maxLS = 10000,
                         funObj = make_phi_alpha(x, fg, pv, calc_gradient_default = TRUE),
                         varargin = NULL,
                         pnorm_inf = max(abs(pv)),
                         progTol = 1e-9, debug = debug)

  res$step$par <- x + res$step$alpha * pv
  res
}

## Test data from the More'-Thuente paper.

# Test Grad/Cubic interpolation (default)

# Table 1
test_that("Table 1 Grad/Cubic Interpolation", {
  pv1 <- -f1$gr(0)/abs(f1$gr(0))
  res11 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e-3, c1 = 0.001, c2 = 0.1, LS_interp = 2)
  expect_step(res11, x = 1.4218, f = -0.35355, df = 0.0013375, alpha = 1.4218, nfev = 6)
  res12 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e-1, c1 = 0.001, c2 = 0.1, LS_interp = 2)
  expect_step(res12, x = 1.4218, f = -0.35355, df = 0.0013220, alpha = 1.4218, nfev = 4)
  res13 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e1, c1 = 0.001, c2 = 0.1, LS_interp = 2)
  expect_step(res13, x = 10, f = -0.098039, df =  0.0094195, alpha = 10, nfev = 1)
  # progtol
  res14 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e3, c1 = 0.001, c2 = 0.1, LS_interp = 2)
  expect_step(res14, x = 333.34, f = -0.0029999, df = 8.9994e-006, alpha = 333.34, nfev = 37)
})

# Table 2
test_that("Table 2 Grad/Cubic Interpolation", {
  pv2 <- -f2$gr(0)/abs(f2$gr(0))
  res21 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e-3, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res21, x = 1.5960, f = -2.6214, df = 2.1828e-014, alpha = 1.5960, nfev = 13)
  res22 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e-1, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res22, x = 1.5960, f = -2.6214, df = -9.4587e-013, alpha = 1.5960, nfev = 10)
  res23 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e1, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res23, x = 1.5960, f = -2.6214, df = 7.2760e-014, alpha = 1.5960, nfev = 9)
  res24 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e3, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res24, x = 1.5960, f = -2.6214, df = 6.8758e-012, alpha = 1.5960, nfev = 14)
})

# Table 3
test_that("Table 3 Grad/Cubic Interpolation", {
  pv3 <- -f3$gr(0)/abs(f3$gr(0))
  res31 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e-3, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res31, x = 1.0, f = -0.011160, df = 8.6906e-006, alpha = 1.0, nfev = 14)
  res32 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e-1, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res32, x = 1.0, f = -0.011160, df = -1.9129e-005, alpha = 1.0, nfev = 11)
  res33 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e1, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res33, x = 1.0, f = -0.011160, df = -2.3751e-006, alpha = 1.0, nfev = 9)
  res34 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e3, c1 = 0.1, c2 = 0.1, LS_interp = 2)
  expect_step(res34, x = 1.0, f = -0.011160, df =  -4.1259e-006, alpha = 1.0, nfev = 12)
})

# Table 4
test_that("Table 4 Grad/Cubic Interpolation", {
  pv4 <- -f4$gr(0)/abs(f4$gr(0))
  res41 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res41, x = 0.030526, f = 0.99902, df = -5.3509e-004, alpha = 0.030526, nfev = 4)
  res42 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res42, x = 0.1, f = 0.99901, df = -4.9330e-005, alpha = 0.1, nfev = 1)
  res43 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res43, x = 0.34963, f = 0.99900, df = -2.9053e-006, alpha = 0.34963, nfev = 3)
  res44 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res44, x = 0.90475, f = 0.99901, df = 5.4438e-005, alpha = 0.90475, nfev = 4)
})

test_that("Table 5 Grad/Cubic Interpolation", {
  pv5 <- -f5$gr(0)/abs(f5$gr(0))
  res51 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res51, x = 0.074271, f = 0.99138, df = 1.7174e-005, alpha = 0.074271, nfev = 7)
  res52 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res52, x = 0.071748, f = 0.99139, df = -6.1309e-004, alpha = 0.071748, nfev = 3)
  res53 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res53, x = 0.075308, f = 0.99138, df = 2.5835e-004, alpha = 0.075308, nfev = 6)
  res54 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res54, x = 0.073117, f = 0.99138, df = -2.6323e-004, alpha = 0.073117, nfev = 8)
})

test_that("Table 6 Grad/Cubic Interpolation", {
  pv6 <- -f6$gr(0)/abs(f6$gr(0))
  res61 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res61, x = 0.92579, f = 0.99138, df = -3.5017e-006, alpha = 0.92579, nfev = 8)
  res62 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res62, x = 0.92706, f = 0.99139, df = 3.0799e-004, alpha = 0.92706, nfev = 11)
  res63 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res63, x = 0.92593, f = 0.99138, df = 3.1671e-005, alpha = 0.92593, nfev = 7)
  res64 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 2)
  expect_step(res64, x = 0.92918, f = 0.99139, df = 8.6261e-004, alpha = 0.92918, nfev = 10)
})

test_that("MT Function modification Grad/Cubic Interpolation", {
  # progtol
  res4m <- mfls(fg = f4, x = 1, pv = -f4$gr(1)/abs(f4$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 2)
  expect_step(res4m, x = 0.5, f = 0.999, df = 1.3132e-011, alpha = 0.5, nfev = 12)
  # progtol
  res5m <- mfls(fg = f5, x = 1, pv = -f5$gr(1)/abs(f5$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 2)
  expect_step(res5m, x = 0.50112, f = 0.99464, df = 0.0087536, alpha = 0.49888, nfev = 14)
  # progtol
  res6m <- mfls(fg = f6, x = 1, pv = -f6$gr(1)/abs(f6$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 2)
  expect_step(res6m, x = 0.82848, f = 0.99188, df = -0.0072578, alpha = 0.17152, nfev = 14)
})


# Test Mixed Quadratic/Cubic interpolation

# Table 1
test_that("Table 1 Quadratic/Cubic interpolation", {
  pv1 <- -f1$gr(0)/abs(f1$gr(0))
  res11 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e-3, c1 = 0.001, c2 = 0.1, LS_interp = 3)
  expect_step(res11, x = 1.3014, f = -0.35234, df = -0.022453, alpha = 1.3014, nfev = 6)
  res12 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e-1, c1 = 0.001, c2 = 0.1, LS_interp = 3)
  expect_step(res12, x = 1.3012, f = -0.35233, df = -0.022500, alpha = 1.3012, nfev = 4)
  res13 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e1, c1 = 0.001, c2 = 0.1, LS_interp = 3)
  expect_step(res13, x = 10, f = -0.098039, df =  0.0094195, alpha = 10, nfev = 1)
  # progtol (NB differs from Octave implementation due to singular matrices, nfev = 21)
  res14 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e3, c1 = 0.001, c2 = 0.1, LS_interp = 3)
  expect_step(res14, x = 333.34, f = -0.0029999, df = 8.9994e-006, alpha = 333.34, nfev = 28)
})

test_that("Table 2 Quadratic/Cubic interpolation", {
  pv2 <- -f2$gr(0)/abs(f2$gr(0))
  res21 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e-3, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res21, x = 1.5960, f = -2.6214, df = 3.4393e-011, alpha = 1.5960, nfev = 13)
  res22 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e-1, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res22, x = 1.5960, f = -2.6214, df = 3.2167e-011, alpha = 1.5960, nfev = 11)
  res23 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e1, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res23, x = 1.5960, f = -2.6214, df = -1.2369e-013, alpha = 1.5960, nfev = 15)
  res24 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e3, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res24, x = 1.5960, f = -2.6214, df = 7.2687e-012, alpha = 1.5960, nfev = 16)
})

test_that("Table 3 Quadratic/Cubic interpolation", {
  pv3 <- -f3$gr(0)/abs(f3$gr(0))
  res31 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e-3, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res31, x = 1.0, f = -0.011160, df = -3.9691e-004, alpha = 1.0, nfev = 18)
  res32 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e-1, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res32, x = 1.0, f = -0.011160, df = -4.1227e-15, alpha = 1.0, nfev = 2)
  res33 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e1, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res33, x = 1.0, f = -0.011160, df = -1.2061e-006, alpha = 1.0, nfev = 14)
  res34 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e3, c1 = 0.1, c2 = 0.1, LS_interp = 3)
  expect_step(res34, x = 1.0, f = -0.011160, df = -4.8577e-007, alpha = 1.0, nfev = 13)
})

test_that("Table 4 Quadratic/Cubic interpolation", {
  pv4 <- -f4$gr(0)/abs(f4$gr(0))
  res41 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res41, x = 0.028791, f = 0.99902, df = -6.0151e-004, alpha = 0.028791, nfev = 12)
  res42 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res42, x = 0.1, f = 0.99901, df = -4.9330e-005, alpha = 0.1, nfev = 1)
  res43 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res43, x = 0.34963, f = 0.99900, df = -2.9053e-006, alpha = 0.34963, nfev = 3)
  res44 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res44, x = 0.90475, f = 0.99901, df = 5.4438e-005, alpha = 0.90475, nfev = 4)
})

test_that("Table 5 Quadratic/Cubic interpolation", {
  pv5 <- -f5$gr(0)/abs(f5$gr(0))
  res51 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res51, x = 0.071381, f = 0.99139, df = -7.1011e-004, alpha = 0.071381, nfev = 10)
  res52 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res52, x = 0.077510, f = 0.99139, df = 7.3935e-004, alpha = 0.077510, nfev = 3)
  res53 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res53, x = 0.072834, f = 0.99139, df = -3.3384e-004, alpha = 0.072834, nfev = 7)
  res54 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res54, x = 0.073825, f = 0.99138, df =  -8.9559e-005, alpha = 0.073825, nfev = 13)
})

test_that("Table 6 Quadratic/Cubic interpolation", {
  pv6 <- -f6$gr(0)/abs(f6$gr(0))
  res61 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res61, x = 0.92537, f = 0.99138, df = -1.0200e-004, alpha = 0.92537, nfev = 19)
  res62 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res62, x = 0.92213, f = 0.99139, df = -8.1497e-004, alpha = 0.92213, nfev = 7)
  res63 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res63, x = 0.92555, f = 0.99138, df = -5.9882e-005, alpha = 0.92555, nfev = 7)
  res64 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 3)
  expect_step(res64, x = 0.92174, f = 0.99139, df = -8.9474e-004, alpha = 0.92174, nfev = 5)
})

test_that("MT Function modification Quadratic/Cubic interpolation", {
  # progtol
  res4m <- mfls(fg = f4, x = 1, pv = -f4$gr(1)/abs(f4$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 3)
  expect_step(res4m, x = 0.5, f = 0.999, df = 6.1914e-012, alpha = 0.5, nfev = 11)
  # progtol
  res5m <- mfls(fg = f5, x = 1, pv = -f5$gr(1)/abs(f5$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 3)
  expect_step(res5m, x = 0.50112, f = 0.99464, df = 0.0087536, alpha = 0.49888, nfev = 13)
  # progtol (differs from Octave version due to singular matrix in polyinterp nfev = 13)
  res6m <- mfls(fg = f6, x = 1, pv = -f6$gr(1)/abs(f6$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 3)
  expect_step(res6m, x = 0.82958, f = 0.99187, df = -0.0072358, alpha = 0.17042, nfev = 14)
})


# Interpolation via fixed step multiplication/bisection

# Table 1
test_that("Table 1 Fixed Step Multiplication/Bisection", {
  pv1 <- -f1$gr(0)/abs(f1$gr(0))
  res11 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e-3, c1 = 0.001, c2 = 0.1, LS_interp = 1)
  expect_step(res11, x = 1.5625, f = -0.35180, df = 0.022377, alpha = 1.5625, nfev = 9)
  res12 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e-1, c1 = 0.001, c2 = 0.1, LS_interp = 1)
  expect_step(res12, x = 1.5626, f = -0.35180, df = 0.022377, alpha = 1.5625, nfev = 7)
  res13 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e1, c1 = 0.001, c2 = 0.1, LS_interp = 1)
  expect_step(res13, x = 10, f = -0.098039, df =  0.0094195, alpha = 10, nfev = 1)
  # progtol
  res14 <- mfls(fg = f1, x = 0, pv = pv1, alpha = 1e3, c1 = 0.001, c2 = 0.1, LS_interp = 1)
  expect_step(res14, x = 500, f = -0.002, df = 3.9999e-6, alpha = 500, nfev = 41)
})

# Table 2
test_that("Table 2 Fixed Step Multiplication/Bisection", {
  pv2 <- -f2$gr(0)/abs(f2$gr(0))
  # Differs from Octave implementation in nfev (39 for Octave, reaching progtol)
  # Difference is due to floating point equality comparison in Armijo check
  res21 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e-3, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res21, x = 1.5960, f = -2.6214, df = 5.5237e-008, alpha = 1.5960, nfev = 36)
  # Differs from Octave implementation in nfev (37 for Octave, reaching progtol)
  # Difference is due to floating point equality comparison in Armijo check
  res22 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e-1, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res22, x = 1.5960, f = -2.6214, df = 5.5237e-008, alpha = 1.5960, nfev = 34)
  # progtol
  res23 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e1, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res23, x = 1.5960, f = -2.6214, df = 5.5237e-008, alpha = 1.5960, nfev = 35)
  res24 <- mfls(fg = f2, x = 0, pv = pv2, alpha = 1e3, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res24, x = 1.5960, f = -2.6214, df = 3.8548e-008, alpha = 1.5960, nfev = 37)
})

# Table 3
test_that("Table 3 Fixed Step Multiplication/Bisection", {
  pv3 <- -f3$gr(0)/abs(f3$gr(0))
  res31 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e-3, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res31, x = 1.0, f = -0.011160, df = -4.1227e-15, alpha = 1.0, nfev = 4)
  res32 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e-1, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res32, x = 1.0, f = -0.011160, df = -4.1227e-15, alpha = 1.0, nfev = 2)
  res33 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e1, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res33, x = 1.0, f = -0.011160, df = 6.1283e-4, alpha = 1.0, nfev = 20)
  res34 <- mfls(fg = f3, x = 0, pv = pv3, alpha = 1e3, c1 = 0.1, c2 = 0.1, LS_interp = 1)
  expect_step(res34, x = 1.0, f = -0.011160, df = 3.2556e-4, alpha = 1.0, nfev = 27)
})

# Table 4
test_that("Table 4 Fixed Step Multiplication/Bisection", {
  pv4 <- -f4$gr(0)/abs(f4$gr(0))
  res41 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res41, x = 0.1, f = 0.99901, df = -4.9330e-005, alpha = 0.1, nfev = 3)
  res42 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res42, x = 0.1, f = 0.99901, df = -4.9330e-005, alpha = 0.1, nfev = 1)
  res43 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res43, x = 0.625, f = 0.99900, df = 2.2733e-6, alpha = 0.625, nfev = 5)
  res44 <- mfls(fg = f4, x = 0, pv = pv4, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res44, x = 0.97656, f = 0.99902, df = 9.0755e-4, alpha = 0.97656, nfev = 11)
})

# Table 5
test_that("Table 5 Fixed Step Multiplication/Bisection", {
  pv5 <- -f5$gr(0)/abs(f5$gr(0))
  res51 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res51, x = 0.0775, f = 0.99139, df = 7.3717e-4, alpha = 0.0775, nfev = 5)
  res52 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res52, x = 0.0750, f = 0.99138, df = 1.8775e-4, alpha = 0.0750, nfev = 3)
  res53 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res53, x = 0.078125, f = 0.99139, df = 8.6648e-4, alpha = 0.078125, nfev = 8)
  res54 <- mfls(fg = f5, x = 0, pv = pv5, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res54, x = 0.076294, f = 0.99139, df = 4.7874e-4, alpha = 0.076294, nfev = 17)
})

# Table 6
test_that("Table 6 Fixed Step Multiplication/Bisection", {
  pv6 <- -f6$gr(0)/abs(f6$gr(0))
  res61 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e-3, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res61, x = 0.92266, f = 0.99139, df = -7.0436e-4, alpha = 0.92266, nfev = 11)
  res62 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e-1, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res62, x = 0.92266, f = 0.99139, df = -7.0436e-4, alpha = 0.92266, nfev = 9)
  res63 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e1, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res63, x = 0.92773, f = 0.99139, df = 4.7843e-004, alpha = 0.92773, nfev = 11)
  res64 <- mfls(fg = f6, x = 0, pv = pv6, alpha = 1e3, c1 = 0.001, c2 = 0.001, LS_interp = 1)
  expect_step(res64, x = 0.92316, f = 0.99139, df = -5.9793e-004, alpha = 0.92316, nfev = 18)
})

test_that("Function modification Fixed Step Multiplication/Bisection", {
  # progtol
  res4m <- mfls(fg = f4, x = 1, pv = -f4$gr(1)/abs(f4$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 1)
  expect_step(res4m, x = 0.5, f = 0.999, df = 0, alpha = 0.5, nfev = 31)
  # progtol
  res5m <- mfls(fg = f5, x = 1, pv = -f5$gr(1)/abs(f5$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 1)
  expect_step(res5m, x = 0.5, f = 0.99463, df = 0.0087527, alpha = 0.5, nfev = 31)
  # progtol
  res6m <- mfls(fg = f6, x = 1, pv = -f6$gr(1)/abs(f6$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9, LS_interp = 1)
  expect_step(res6m, x = 0.75, f = 0.99249, df = -0.0081531, alpha = 0.25, nfev = 31)
})
