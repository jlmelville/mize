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

test_that("safe Cholesky repairs non-positive eigenvalues", {
  rotation <- matrix(
    c(
      cos(pi / 6),
      sin(pi / 6),
      -sin(pi / 6),
      cos(pi / 6)
    ),
    nrow = 2
  )
  eps <- 1e-6
  spectra <- list(
    zero = c(2, 0),
    near_zero = c(2, -1e-12),
    indefinite = c(2, -3)
  )
  hessians <- list(
    zero = diag(spectra$zero),
    near_zero = rotation %*% diag(spectra$near_zero) %*% t(rotation),
    indefinite = rotation %*% diag(spectra$indefinite) %*% t(rotation)
  )

  for (case_name in names(hessians)) {
    repaired_values <- pmax(spectra[[case_name]], eps)
    expected <- if (case_name == "zero") {
      diag(repaired_values)
    } else {
      rotation %*% diag(repaired_values) %*% t(rotation)
    }
    factor <- safe_chol(hessians[[case_name]], eps = eps)

    expect_false(is.null(factor), info = case_name)
    if (!is.null(factor)) {
      expect_equal(
        crossprod(factor),
        expected,
        tolerance = 1e-12,
        info = case_name
      )
    }
  }
})

test_that("safe Cholesky honors a non-default eigenvalue floor", {
  eps <- 0.25
  expected <- diag(c(4, eps))
  factor <- safe_chol(diag(c(4, -2)), eps = eps)

  expect_false(is.null(factor))
  expect_equal(crossprod(factor), expected, tolerance = 1e-12)
  expect_equal(crossprod(factor)[2, 2] / eps, 1, tolerance = 1e-12)
})

test_that("safe Cholesky preserves the ordinary positive-definite path", {
  hessian <- matrix(c(4, 1, 1, 3), nrow = 2)
  expected <- chol(hessian)
  factor <- safe_chol(hessian, eps = 10)

  expect_identical(factor, expected)
  expect_equal(crossprod(factor), hessian, tolerance = 1e-12)
})

test_that("Newton safe Cholesky path repairs indefinite Hessians", {
  gradient <- c(1, 1)
  indefinite_fg <- list(
    hs = function(x) {
      matrix(c(1, 0, 0, -1), nrow = 2)
    }
  )
  opt <- list(cache = list(gr_curr = gradient, gr_curr_iter = 1))

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
