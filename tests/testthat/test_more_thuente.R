context("More'-Thuente Line Search")

# This test uses input parameters, directions and step sizes from using the MT
# line search with a few steps of the CG solver. The expected values come
# from plugging the input values into Dianne O'Leary's Matlab code (running
# under GNU Octave).

expect_step <- function(actual, x, f, df, alpha, nfev, tolerance = 1e-4) {
  expect_equal(actual$step$par, x, tolerance = tolerance)
  expect_equal(actual$step$f, f, tolerance = tolerance)
  expect_equal(actual$step$df, df, tolerance = tolerance)
  expect_equal(actual$step$alpha, alpha, tolerance = tolerance)
  expect_equal(actual$nfn, nfev)
}

step <- function(x, f, df, alpha, nfev) {
  list(x = x, step = list(f = f, df = df, alpha = alpha), nfev = nfev)
}

rfn <- rosenbrock_fg$fn
rgr <- rosenbrock_fg$gr

mtls <- function(fn, gr, x, pv, alpha, c1, c2) {
  cvsrch(phi = make_phi(fn, gr, x, pv),
          step0 = make_step0(fn, gr, x, pv), alpha = alpha, c1 = c1, c2 = c2)
}

## These tests are designed to reproduce the data in Tables 1-6
## of the More'-Thuente paper. They do so apart from a small number of minor
## differences (what value differs and when are indicated by comments before the
## test). Note that the values here are also reproduced by the Matlab code by
## O'Leary, i.e. where the R result differs from the published data in the
## original More'-Thuente paper, so does the Matlab code.

# Table 1
test_that("Table 1", {
  pv1 <- -f1$gr(0)/abs(f1$gr(0))
  res11 <- mtls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e-3, c1 = 0.001, c2 = 0.1)
  expect_step(res11, x = 1.3650, f = -0.35333, df = -0.0091645, alpha = 1.3650, nfev = 6)
  res12 <- mtls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e-1, c1 = 0.001, c2 = 0.1)
  expect_step(res12, x = 1.4414, f = -0.35349, df = 0.0046645, alpha = 1.4414, nfev = 3)
  res13 <- mtls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e1, c1 = 0.001, c2 = 0.1)
  expect_step(res13, x = 10, f = -0.098039, df =  0.0094195, alpha = 10, nfev = 1)
  res14 <- mtls(fn = f1$fn, gr = f1$gr, x = 0, pv = pv1, alpha = 1e3, c1 = 0.001, c2 = 0.1)
  expect_step(res14, x =  36.888, f = -0.027070, df = 7.3169e-004, alpha = 36.888, nfev = 4)
})

# Table 2
test_that("Table 2", {
  pv2 <- -f2$gr(0)/abs(f2$gr(0))
  # # gradient 7.1e-9
  res21 <- mtls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(res21, x = 1.5960, f = -2.6214, df = 3.8113e-009, alpha = 1.5960, nfev = 12)
  # gradient 10e-10 could be a typo in the paper and should be 1.0e-10?
  res22 <- mtls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(res22, x = 1.5960, f = -2.6214, df = 1.0106e-010, alpha = 1.5960, nfev = 8)
  res23 <- mtls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(res23, x = 1.5960, f = -2.6214, df = -4.9725e-009, alpha = 1.5960, nfev = 8)
  res24 <- mtls(fn = f2$fn, gr = f2$gr, x = 0, pv = pv2, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(res24, x = 1.5960, f = -2.6214, df = -2.3091e-008, alpha = 1.5960, nfev = 11)
})

# Table 3
test_that("Table 3", {
pv3 <- -f3$gr(0)/abs(f3$gr(0))
res31 <- mtls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
expect_step(res31, x = 1.0, f = -0.011160, df = -5.1440e-005, alpha = 1.0, nfev = 12)
res32 <- mtls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
expect_step(res32, x = 1.0, f = -0.011160, df = -1.9224e-004, alpha = 1.0, nfev = 12)
res33 <- mtls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e1, c1 = 0.1, c2 = 0.1)
expect_step(res33, x = 1.0, f = -0.011160, df = -1.9892e-006, alpha = 1.0, nfev = 10)
res34 <- mtls(fn = f3$fn, gr = f3$gr, x = 0, pv = pv3, alpha = 1e3, c1 = 0.1, c2 = 0.1)
expect_step(res34, x = 1.0, f = -0.011160, df = -1.5789e-005, alpha = 1.0, nfev = 13)
})

# Table 4
test_that("Table 4", {
pv4 <- -f4$gr(0)/abs(f4$gr(0))
# alpha = 0.08
res41 <- mtls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
expect_step(res41, x = 0.085, f = 0.99901, df = -6.8531e-005, alpha = 0.085, nfev = 4)
res42 <- mtls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
expect_step(res42, x = 0.1, f = 0.99901, df = -4.9330e-005, alpha = 0.1, nfev = 1)
res43 <- mtls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e1, c1 = 0.001, c2 = 0.001)
expect_step(res43, x = 0.34910, f = 0.999, df = -2.9195e-006, alpha = 0.34910, nfev = 3)
res44 <- mtls(fn = f4$fn, gr = f4$gr, x = 0, pv = pv4, alpha = 1e3, c1 = 0.001, c2 = 0.001)
expect_step(res44, x = 0.8294, f = 0.999, df = 1.6436e-005, alpha = 0.8294, nfev = 4)
})

# Table 5
test_that("Table 5", {
pv5 <- -f5$gr(0)/abs(f5$gr(0))
res51 <- mtls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
expect_step(res51, x = 0.075011, f = 0.99138, df = 1.9025e-004, alpha = 0.075011, nfev = 6)
res52 <- mtls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
expect_step(res52, x = 0.07751, f = 0.99139, df = 7.3935e-004, alpha = 0.07751, nfev = 3)
res53 <- mtls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e1, c1 = 0.001, c2 = 0.001)
expect_step(res53, x = 0.073142, f = 0.99138, df = -2.5691e-004, alpha = 0.073142, nfev = 7)
res54 <- mtls(fn = f5$fn, gr = f5$gr, x = 0, pv = pv5, alpha = 1e3, c1 = 0.001, c2 = 0.001)
expect_step(res54, x = 0.076159, f = 0.99139, df = 4.4913e-004, alpha = 0.076159, nfev = 8)
})

# Table 6
test_that("Table 6", {
pv6 <- -f6$gr(0)/abs(f6$gr(0))
res61 <- mtls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
expect_step(res61, x = 0.9279, f = 0.99139, df = 5.2203e-004, alpha = 0.9279, nfev = 13)
res62 <- mtls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
expect_step(res62, x = 0.92615, f = 0.99138, df = 8.3588e-005, alpha = 0.92615, nfev = 11)
res63 <- mtls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e1, c1 = 0.001, c2 = 0.001)
expect_step(res63, x = 0.92478, f = 0.99138, df = -2.3788e-004, alpha = 0.92478, nfev = 8)
res64 <- mtls(fn = f6$fn, gr = f6$gr, x = 0, pv = pv6, alpha = 1e3, c1 = 0.001, c2 = 0.001)
expect_step(res64, x = 0.92440, f = 0.99139, df = -3.2498e-004, alpha = 0.92440, nfev = 11)
})

# The above tests don't enter the code path where the function is modified much. The
# tests below do exercise that part.
test_that("Function modification", {
  res4m <- mtls(fn = f4$fn, gr = f4$gr, x = 1, pv = -f4$gr(1)/abs(f4$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9)
expect_step(res4m, x = 0.99615, f = 0.99913, df = 0.032049, alpha = 0.0038522, nfev = 6)
res5m <- mtls(fn = f5$fn, gr = f5$gr, x = 1, pv = -f5$gr(1)/abs(f5$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9)
expect_step(res5m, x = 0.99599, f = 0.99914, df = 0.038284, alpha = 0.0040126, nfev = 6)
res6m <- mtls(fn = f6$fn, gr = f6$gr, x = 1, pv = -f6$gr(1)/abs(f6$gr(1)), alpha = 1, c1 = 0.1, c2 = 0.9)
expect_step(res6m, x = 0.95655, f = 0.99157, df = 0.016504, alpha = 0.043447, nfev = 4)
})
