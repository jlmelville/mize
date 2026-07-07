test_that("Newton with basic backtracking works ok", {
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = newton_direction(),
        step_size = backtracking(rho = 0.5, c1 = 1e-4)
      ),
      verbose = FALSE
    )
  )
  opt$count_res_fg <- FALSE

  res <- opt_loop(
    opt,
    rb0,
    rosenbrock_fg,
    20,
    store_progress = TRUE,
    verbose = FALSE,
    grad_tol = 1e-5
  )

  nfs <- c(
    0,
    2,
    6,
    7,
    8,
    9,
    12,
    13,
    14,
    15,
    17,
    18,
    19,
    20,
    22,
    23,
    24,
    25,
    26,
    27,
    28
  )
  ngs <- c(
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20
  )
  fs <- c(
    24.2,
    4.73,
    4.09,
    3.23,
    3.21,
    1.94,
    1.60,
    1.18,
    0.922,
    0.597,
    0.453,
    0.281,
    0.211,
    0.089,
    0.0515,
    0.0200,
    0.00717,
    0.00107,
    7.78e-5,
    2.82e-7,
    8.42e-12
  )
  g2ns <- c(
    232.87,
    4.64,
    28.55,
    11.57,
    30.33,
    3.60,
    9.25,
    4.92,
    8.66,
    1.78,
    5.88,
    2.18,
    9.40,
    0.492,
    3.92,
    1.24,
    2.53,
    0.238,
    0.348,
    0.00387,
    1.187e-4
  )
  steps <- c(
    0,
    0.381,
    0.619,
    0.268,
    0.580,
    0.0793,
    0.233,
    0.146,
    0.208,
    0.112,
    0.165,
    0.143,
    0.248,
    0.107,
    0.174,
    0.125,
    0.160,
    0.0608,
    0.0621,
    0.00801,
    0.00117
  )
  alphas <- c(
    0,
    1,
    0.125,
    1,
    1,
    1,
    0.25,
    1,
    1,
    1,
    0.5,
    1,
    1,
    1,
    0.5,
    1,
    1,
    1,
    1,
    1,
    1
  )
  par <- c(1, 1)
  expect_equal(res$progress$nf, nfs)
  expect_equal(res$progress$ng, ngs)
  expect_equal(res$progress$f, fs, tolerance = 1e-3)
  expect_equal(res$progress$g2n, g2ns, tolerance = 1e-3)
  expect_equal_abs(res$progress$step, steps, tolerance = 1e-3)
  expect_equal(res$progress$alpha, alphas, tolerance = 1e-3)
  expect_equal(res$par, par, tolerance = 1e-3)
})

test_that("internal backtracking uses a per-search function budget", {
  # The public backtracking line search uses Schmidt Armijo. This test covers
  # the legacy internal helper because it is still used by direct optimizer
  # construction tests and should not confuse cumulative counts with local
  # line-search budget.
  quadratic_fg <- list(
    fn = function(x) {
      x[1]^2
    },
    gr = function(x) {
      2 * x
    }
  )
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(),
        step_size = backtracking(rho = 0.5, c1 = 1e-4, max_fn = 1)
      )
    )
  )
  opt$convergence <- list(max_fn = 103, max_fg = Inf)
  opt$counts$fn <- 100
  opt <- set_gr_curr(opt, quadratic_fg$gr(c(1)), 1)

  stage <- opt$stages$gradient_descent
  stage$direction$value <- -opt$cache$gr_curr
  sub_stage <- stage$step_size

  res <- sub_stage$calculate(opt, stage, sub_stage, c(1), quadratic_fg, 1)

  expect_equal(res$sub_stage$value, 0.5)
  expect_equal(res$opt$counts$fn, 103)
})

# From https://math.stackexchange.com/questions/1130002/newton-optimization-algorithm-with-non-positive-definite-hessian
# At starting location, Hessian is not positive definite and cholesky fails
# Details of minimization probably aren't important, but it should reach
# the minimum at (3, 0.5) without exploding
test_that("Newton method can survive non-positive definite Hessian", {
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = newton_direction(),
        step_size = more_thuente_ls(initial_step_length = 1)
      ),
      verbose = FALSE
    )
  )

  res <- opt_loop(
    opt,
    c(4, 1),
    tricky_fg(),
    9,
    store_progress = TRUE,
    verbose = FALSE
  )

  par <- c(3, 0.5)

  expect_equal(res$par, par, tolerance = 1e-3)
})

test_that("Newton safe Cholesky path repairs indefinite Hessians", {
  gradient <- c(1, 1)
  indefinite_fg <- list(
    hs = function(x) {
      matrix(c(1, 0, 0, -1), nrow = 2)
    }
  )
  opt <- list(cache = list(gr_curr = gradient))

  safe_direction <- newton_direction(try_safe_chol = TRUE)
  safe_res <- safe_direction$calculate(
    opt,
    stage = list(),
    sub_stage = safe_direction,
    par = c(0, 0),
    fg = indefinite_fg,
    iter = 1
  )

  plain_direction <- newton_direction(try_safe_chol = FALSE)
  plain_res <- plain_direction$calculate(
    opt,
    stage = list(),
    sub_stage = plain_direction,
    par = c(0, 0),
    fg = indefinite_fg,
    iter = 1
  )

  expect_equal(plain_res$sub_stage$value, -gradient)
  expect_equal(safe_res$sub_stage$value, c(-1, -1e10), tolerance = 1e-12)
  expect_lt(dot(gradient, safe_res$sub_stage$value), 0)
})
