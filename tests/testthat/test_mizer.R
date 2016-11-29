test_that("steepest descent with constant step size", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(),
        step_size = constant_step_size(
          value = 0.0001)),
    verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 0, 0, 0)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 19.18, 15.52, 12.81)
  g2ns <- c(232.87, 198.56, 170.43, 147.11)
  steps <- c(0, 0.0233, 0.0199, 0.017) # step size is multiplied by gradient!
  par <- c(-1.144, 1.023)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("steepest descent with constant step size and normalized direction", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = constant_step_size(
          value = 0.0001)),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 0, 0, 0)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 24.18, 24.15, 24.13)
  g2ns <- c(232.87, 232.72, 232.57, 232.42)
  steps <- c(0, 1e-4, 1e-4, 1e-4)
  par <- c(-1.2, 1.0)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("steepest descent with bold driver", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver(
          init_step_size = 1)),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 7, 12)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 4.12, 4.11)
  g2ns <- c(232.87, 64.72, 2.90, 2.39)
  steps <- c(0, 0.25, 0.069, 0.0047)
  par <- c(-1.024, 1.060)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})


test_that("classical momentum with 0 step size should be like using no momentum", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = constant_step_size(
          value = 0.01)),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0
        )),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 0, 0, 0)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 21.95, 19.84, 17.88)
  g2ns <- c(232.87, 217.96, 203.31, 188.93)
  steps <- c(0, 1e-2, 1e-2, 1e-2)
  par <- c(-1.172, 1.011)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("classical momentum with constant step size", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = constant_step_size(
          value = 0.01)),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 0, 0, 0)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 21.95, 19.44, 17.06)
  g2ns <- c(232.87, 217.96, 200.42, 182.69)
  steps <- c(0, 0.01, 0.012, 0.0124)
  par <- c(-1.168, 1.013)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})


test_that("eager classical momentum with constant step size should give same results as non-eager", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = constant_step_size(
          value = 0.01)),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  opt$eager_update <- TRUE
  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 0, 0, 0)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 21.95, 19.44, 17.06)
  g2ns <- c(232.87, 217.96, 200.42, 182.69)
  steps <- c(0, 0.01, 0.012, 0.0124)
  par <- c(-1.168, 1.013)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})


test_that("classical momentum with bold driver", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 8, 10)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 5.25, 4.62)
  g2ns <- c(232.87, 64.72, 47.19, 33.88)
  steps <- c(0, 0.25, 0.020, 0.0795)
  par <- c(-1.051, 1.040)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("eager classical momentum with bold driver same as 'lazy' result", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  opt$eager_update <- TRUE
  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 8, 10)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 5.25, 4.62)
  g2ns <- c(232.87, 64.72, 47.19, 33.88)
  steps <- c(0, 0.25, 0.020, 0.0795)
  par <- c(-1.051, 1.040)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("linear weighted classical momentum with bold driver", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  opt$stages <- c(opt$stages, momentum_correction_stage())

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 10, 12)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 4.27, 5.04, 4.67)
  g2ns <- c(232.87, 17.16, 42.78, 33.27)
  steps <- c(0, 0.2, 0.027, 0.01)
  par <- c(-0.998, 1.078)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("linear weighted eager classical momentum with bold driver", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  opt$stages <- c(opt$stages, momentum_correction_stage())
  opt$eager_update <- TRUE

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 10, 12)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 4.27, 5.04, 4.67)
  g2ns <- c(232.87, 17.16, 42.78, 33.27)
  steps <- c(0, 0.2, 0.027, 0.01)
  par <- c(-0.998, 1.078)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("bold classical momentum with bold driver", {
  # Not a very good idea - momentum component easily shrinks to zero
  # so you waste a lot of time trying to find a non-existent acceptable step
  # size - but tests that you can use the same step size method in different
  # stages and they don't interfere with each other.
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(normalize = TRUE),
        step_size = bold_driver()),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 6, 19, 45)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 4.12, 4.10)
  g2ns <- c(232.87, 64.72, 2.81, 2.41)
  steps <- c(0, 0.25, 0.064, 0.0047)
  par <- c(-1.027, 1.059)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("bold classical momentum with bold driver without cache gives same results, but requires extra work", {
  # Checks that the caching of function calls works correctly
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(normalize = TRUE),
        step_size = bold_driver()),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE,
                 invalidate_cache = TRUE)

  nfs <- c(0, 6, 20, 47) # extra function evaluations
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 4.12, 4.10)
  g2ns <- c(232.87, 64.72, 2.81, 2.41)
  steps <- c(0, 0.25, 0.064, 0.0047)
  par <- c(-1.027, 1.059)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})


test_that("classical momentum with bold driver and fn adaptive restart, same results as without when everything is ok", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  opt <- adaptive_restart(opt, "fn")

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  # Have to carry out one extra fn evaluation when doing the first validation
  # After that, bold driver doesn't need to do a calculation for f0 on the
  # subsequent step
  nfs <- c(0, 5, 9, 11)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 5.25, 4.62)
  g2ns <- c(232.87, 64.72, 47.19, 33.88)
  steps <- c(0, 0.25, 0.020, 0.0795)
  par <- c(-1.051, 1.040)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("classical momentum with bold driver and gr adaptive restart, same results as without when everything is ok", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  opt <- adaptive_restart(opt, "gr")

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  # Get adaptive update check for free if using gradient!
  nfs <- c(0, 4, 8, 10)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 5.25, 4.62)
  g2ns <- c(232.87, 64.72, 47.19, 33.88)
  steps <- c(0, 0.25, 0.020, 0.0795)
  par <- c(-1.051, 1.040)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("classical momentum with bold driver aggressive momentum can cause cost increase", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.4
        )),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 8, 10)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 8.71, 4.69)
  g2ns <- c(232.87, 64.72, 91.13, 34.39)
  steps <- c(0, 0.25, 0.033, 0.064)
  par <- c(-0.989, 1.064)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("classical momentum with bold driver adaptive gr momentum prevents cost increase", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.4
        )),
      verbose = FALSE))

  opt <- adaptive_restart(opt, "gr")

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 8, 10)
  ngs <- c(0, 1, 2, 2) # no grad calc needed on repeated step
  fs <- c(24.2, 6.32, 6.32, 4.16)
  g2ns <- c(232.87, 64.72, 64.72, 9.70)
  steps <- c(0, 0.25, 0, 0.076)
  par <- c(-1.035, 1.058)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("classical momentum with bold driver adaptive fn momentum prevents cost increase", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.4
        )),
      verbose = FALSE))

  opt <- adaptive_restart(opt, "fn")

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 5, 9, 11)
  ngs <- c(0, 1, 2, 2) # no grad calc needed on repeated step
  fs <- c(24.2, 6.32, 6.32, 4.16)
  g2ns <- c(232.87, 64.72, 64.72, 9.70)
  steps <- c(0, 0.25, 0, 0.076)
  par <- c(-1.035, 1.058)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("nesterov momentum with bold driver", {

  opt <- make_opt(
    make_stages(
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      verbose = FALSE))

  opt$eager_update <- TRUE

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 4, 7, 10)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 4.33, 4.78)
  g2ns <- c(232.87, 64.72, 22.22, 36.89)
  steps <- c(0, 0.25, 0.088, 0.0584)
  par <- c(-0.987, 1.065)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("nesterov momentum with bold driver and adaptive fn", {

  opt <- make_opt(
    make_stages(
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = bold_driver()),
      verbose = FALSE))
  opt$eager_update <- TRUE

  opt <- adaptive_restart(opt, "fn")

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 5, 8, 11)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 6.32, 4.33, 4.33)
  g2ns <- c(232.87, 64.72, 22.22, 22.22)
  steps <- c(0, 0.25, 0.088, 0)
  par <- c(-1.042, 1.046)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("nesterov accelerated gradient with wolfe line search", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(),
        step_size = more_thuente_ls(c2 = 1.e-9)),
      momentum_stage(
        direction = nesterov_momentum_direction(),
        step_size = nesterov_convex_step()
      ),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 9, 17, 22)
  ngs <- c(0, 9, 17, 22)
  fs <- c(24.2, 4.128, 3.913, 3.558)
  g2ns <- c(232.87, 1.777, 23.908, 7.200)
  steps <- c(0, 0.184, 0.301, 0.048)
  par <- c(-0.869, 0.781)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

# Delta-Bar-Delta
test_that("classical momentum with approximate delta-bar-delta", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(normalize = TRUE),
        step_size = jacobs(epsilon = 0.1)),
      momentum_stage(
        direction = momentum_direction(),
        step_size = constant_step_size(
          value = 0.2
        )),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 0, 0, 0)
  ngs <- c(0, 1, 2, 3)
  fs <- c(24.2, 7.10, 6.53, 4.84)
  g2ns <- c(232.87, 83.18, 67.62, 37.87)
  steps <- c(0, 0.11, 0.143, 0.032)
  par <- c(-0.993, 1.080)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

# Wolfe line search
test_that("Polak Ribiere CG with Rasmussen LS", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = cg_direction(cg_update = pr_update),
        step_size = rasmussen_ls(initial_step_length = 1.8440535542657e-05)),
      verbose = FALSE))


  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 6, 10, 12)
  ngs <- c(0, 6, 10, 12)
  fs <- c(24.2, 4.13, 3.84, 3.52)
  g2ns <- c(232.87, 1.87, 19.06, 25.24)
  steps <- c(0, 0.184, 0.273, 0.311)
  par <- c(-0.777, 0.543)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})

test_that("BFGS with More-Thuente LS", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = bfgs_direction(),
        step_size = more_thuente_ls(
          c2 = 1e-9,
          initial_step_length = 0.0043372269,
          initializer = "q")),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 6, 11, 17)
  ngs <- c(0, 6, 11, 17)
  fs <- c(24.2, 4.13, 3.85, 3.53)
  g2ns <- c(232.87, 1.78, 18.62, 24.98)
  steps <- c(0, 0.184, 0.261, 0.307)
  par <- c(-0.785, 0.558)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})


test_that("L-BFGS with More-Thuente LS gives same results as BFGS with sufficiently high memory", {

  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = lbfgs_direction(),
        step_size = more_thuente_ls(
          c2 = 1e-9,
          initial_step_length = 0.0043372269,
          initializer = "q")),
      verbose = FALSE))

  res <- optloop(opt, out0, rosenbrock_fg$fn, rosenbrock_fg$gr, 3,
                 store_progress = TRUE, verbose = FALSE)

  nfs <- c(0, 6, 11, 17)
  ngs <- c(0, 6, 11, 17)
  fs <- c(24.2, 4.13, 3.85, 3.53)
  g2ns <- c(232.87, 1.78, 18.62, 24.98)
  steps <- c(0, 0.184, 0.261, 0.307)
  par <- c(-0.785, 0.558)

  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tol = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tol = 1e-3)
  expect_equal(res$progress$step, steps, tol = 1e-3)
  expect_equal(res$par, par, tol = 1e-3)
})


