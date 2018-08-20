context("Truncated Newton")

# Numbers from Octave minfunc with settings:
# options.c1 = 1e-10, options.c2 = 1e-9, options.Method = 'newton0'
# options.maxFunEvals = 31
res <- mize(c(0 ,0), rosenbrock_fg, method = "TN", c1 = 1e-10, c2 = 1e-9,
            store_progress = TRUE, max_fn = 27)

expect_equal(res$nf, 27)
expect_equal(res$ng, 34)
expect_equal(res$f, 0.35, tol = 1e-3)
expect_equal(res$par, c(0.4446, 0.1772), tol = 1e-3)

expect_equal(res$progress$alpha,
             c(0, 1.61262e-01, 1.00000e+00, 1.79967e-01, 1.00453e+00,
               2.73325e-01),
             tol = 1e-5)
expect_equal(res$progress$f,
             c(1, 7.71110e-01, 7.03481e-01, 5.24969e-01, 4.81177e-01,
               3.50106e-01),
             tol = 1e-5)
