context("Newton Method")

# From https://math.stackexchange.com/questions/1130002/newton-optimization-algorithm-with-non-positive-definite-hessian
# At starting location, Hessian is not positive definite and cholesky fails
# Details of minimization probably aren't important, but it should reach
# the minimum at (3, 0.5) without exploding
test_that("Newton method can survive non-positive definite Hessian", {
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = newton_direction(),
        step_size = more_thuente_ls(initial_step_length = 1)),
      verbose = FALSE))

  res <- optloop(opt, c(4, 1), tricky_fg(), 9,
                 store_progress = TRUE, verbose = FALSE)

  par <- c(3, 0.5)

  expect_equal(res$par, par, tol = 1e-3)
})


