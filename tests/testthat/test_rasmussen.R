context("Rasmussen Line Search")

expect_step <- function(actual, x, f, df, alpha, nfev, tolerance = 1e-4) {
  expect_equal(actual$step$par, x, tolerance = tolerance)
  expect_equal(actual$step$f, f, tolerance = tolerance)
  expect_equal(actual$step$df, df, tolerance = tolerance)
  expect_equal(actual$step$alpha, alpha, tolerance = tolerance)
  expect_equal(actual$nfn, nfev)
}

rls <- function(fn, gr, x, pv, alpha, c1, c2) {
  ras_ls(phi = make_phi(fn, gr, x, pv),
                   alpha,
                   step0 = make_step0(fn, gr, x, pv),
                   max_fn = 10000,
                  c1 = c1, c2 = c2)
}

## Test data from the More'-Thuente paper.

# Table 1
test_that("Table 1", {
  pv1 <- -f1$gr(0)/abs(f1$gr(0))
  res11 <- rls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e-3, c1 = 0.001, c2 = 0.1)
  expect_step(res11, x = 1.2132, f = -0.34944, df = -0.043803, alpha = 1.2132, nfev = 9)
  res12 <- rls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e-1, c1 = 0.001, c2 = 0.1)
  expect_step(res12, x = 1.2531, f = -0.35098, df = -0.033723, alpha = 1.2531, nfev = 5)
  res13 <- rls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e1, c1 = 0.001, c2 = 0.1)
  expect_step(res13, x = 10, f = -0.098039, df =  0.0094195, alpha = 10, nfev = 1)
  res14 <- rls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e3, c1 = 0.001, c2 = 0.1)
  expect_step(res14, x = 37.054, f = -0.026948, df = 7.2516e-004, alpha = 37.054, nfev = 4)
})

# Table 2
test_that("Table 2", {
  pv2 <- -f2$gr(0)/abs(f2$gr(0))
  res21 <- rls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(res21, x = 1.5960, f = -2.6214, df = 0, alpha = 1.5960, nfev = 55)
  res22 <- rls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(res22, x = 1.5960, f = -2.6214, df = -4.2819e-010, alpha = 1.5960, nfev = 38)
  res23 <- rls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(res23, x = 1.5960, f = -2.6214, df = -1.2233e-010, alpha = 1.5960, nfev = 10)
  res24 <- rls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(res24, x = 1.5960, f = -2.6214, df = -1.2233e-010, alpha = 1.5960, nfev = 12)
})

# Table 3
test_that("Table 3", {
  pv3 <- -f3$gr(0)/abs(f3$gr(0))
  res31 <- rls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(res31, x = 1.0, f = -0.011160, df = 2.3645e-006, alpha = 1.0, nfev = 18)
  res32 <- rls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(res32, x = 1.0, f = -0.011160, df = -3.8131e-006, alpha = 1.0, nfev = 15)
  res33 <- rls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(res33, x = 1.0, f = -0.011160, df = -4.1227e-015, alpha = 1.0, nfev = 2)
  res34 <- rls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(res34, x = 1.0, f = -0.011160, df = -4.1227e-005, alpha = 1.0, nfev = 4)
})

# Table 4
test_that("Table 4", {
  pv4 <- -f4$gr(0)/abs(f4$gr(0))
  res41 <- rls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(res41, x = 0.023344, f = 0.99902, df = -9.1485e-004, alpha = 0.023344, nfev = 13)
  res42 <- rls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(res42, x = 0.1, f = 0.99901, df = -4.9330e-005, alpha = 0.1, nfev = 1)
  res43 <- rls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(res43, x = 0.47507, f = 0.999002, df = -4.0043e-007, alpha = 0.47507, nfev = 3)
  res44 <- rls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(res44, x = 0.9235, f = 0.99901, df = 8.4666e-005, alpha = 0.9235, nfev = 5)
})

# Table 5
test_that("Table 5", {
  pv5 <- -f5$gr(0)/abs(f5$gr(0))
  res51 <- rls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(res51, x = 0.074201, f = 0.99138, df = 5.3016e-007, alpha = 0.074201, nfev = 8)
  res52 <- rls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(res52, x = 0.07175, f = 0.99139, df = -6.1309e-004, alpha = 0.07175, nfev = 3)
  res53 <- rls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(res53, x = 0.074060, f = 0.99138, df = -3.3023e-005, alpha = 0.074060, nfev = 7)
  res54 <- rls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(res54, x = 0.073447, f = 0.99138, df = -1.8157e-004, alpha = 0.073447, nfev = 9)
})

# Table 6
test_that("Table 6", {
  pv6 <- -f6$gr(0)/abs(f6$gr(0))
  res61 <- rls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(res61, x = 0.9296, f = 0.99139, df = 9.7500e-004, alpha = 0.9296, nfev = 56)
  res62 <- rls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(res62, x = 0.92662, f = 0.99138, df = 1.9698e-004, alpha = 0.92662, nfev = 30)
  res63 <- rls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(res63, x = 0.92966, f = 0.99139, df = 9.9413e-004, alpha = 0.92966, nfev = 7)
  res64 <- rls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(res64, x = 0.92436, f = 0.99139, df = -3.3353e-004, alpha = 0.92434, nfev = 8)
})

test_that("Function modification", {
  res4m <- rls(fn = f4$fn, gr = f4$gr, x = 1, pv = -f4$gr(1)/abs(f4$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(res4m, x = 0.99278, f = 0.99907, df = 0.009454, alpha = 0.0072168, nfev = 6)
  res5m <- rls(fn = f5$fn, gr = f5$gr, x = 1, pv = -f5$gr(1)/abs(f5$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(res5m, x = 0.99243, f = 0.99905, df = 0.017425, alpha = 0.0075707, nfev = 6)
  res6m <- rls(fn = f6$fn, gr = f6$gr, x = 1, pv = -f6$gr(1)/abs(f6$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(res6m, x = 0.936501, f = 0.99140, df = 0.0032111, alpha = 0.063499, nfev = 4)
})
