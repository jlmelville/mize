run_schmidt_wolfe_oracle <- function(
  objective,
  initial_parameter,
  initial_alpha,
  armijo_constant,
  curvature_constant,
  direction = -objective$gr(initial_parameter) /
    abs(objective$gr(initial_parameter)),
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE
) {
  initial_step <- make_step0(objective, initial_parameter, direction)
  result <- new_schmidt_wolfe_search(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    max_evaluations = 10000,
    approximation_tolerance = approximation_tolerance,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature
  )(
    phi = make_phi_alpha(
      initial_parameter,
      objective,
      direction,
      calc_gradient_default = TRUE
    ),
    step0 = initial_step,
    alpha = initial_alpha,
    pm = direction
  )

  result$step$par <- initial_parameter + result$step$alpha * direction
  result
}

expect_schmidt_wolfe_oracle <- function(result, expected, info) {
  expect_equal(
    result$step$par,
    expected$parameter,
    tolerance = 1e-4,
    info = info
  )
  expect_equal(
    result$step$f,
    expected$value,
    tolerance = 1e-4,
    info = info
  )
  expect_equal_abs(
    result$step$df,
    expected$gradient,
    tolerance = 1e-4,
    info = info
  )
  expect_equal(
    result$step$alpha,
    expected$alpha,
    tolerance = 1e-4,
    info = info
  )
  expect_equal(result$nfn, expected$evaluations, info = info)
}

schmidt_wolfe_table_oracles <- list(
  table_1 = list(
    objective = f1,
    armijo_constant = 0.001,
    curvature_constant = 0.1,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(1.4218, 1.4218, 10, 333.34),
      value = c(-0.35355, -0.35355, -0.098039, -0.0029999),
      gradient = c(0.0013375, 0.0013220, 0.0094195, 8.9994e-6),
      evaluations = c(6L, 4L, 1L, 37L)
    )
  ),
  table_2 = list(
    objective = f2,
    armijo_constant = 0.1,
    curvature_constant = 0.1,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = rep(1.5960, 4),
      value = rep(-2.6214, 4),
      gradient = c(2.1828e-14, -9.4587e-13, 7.2760e-14, 6.8758e-12),
      evaluations = c(13L, 10L, 9L, 14L)
    )
  ),
  table_3 = list(
    objective = f3,
    armijo_constant = 0.1,
    curvature_constant = 0.1,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = rep(1, 4),
      value = rep(-0.011160, 4),
      gradient = c(8.6906e-6, -1.9129e-5, -2.3751e-6, -4.1259e-6),
      evaluations = c(14L, 11L, 9L, 12L)
    )
  ),
  table_4 = list(
    objective = f4,
    armijo_constant = 0.001,
    curvature_constant = 0.001,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(0.030526, 0.1, 0.34963, 0.90475),
      value = c(0.99902, 0.99901, 0.99900, 0.99901),
      gradient = c(-5.3509e-4, -4.9330e-5, -2.9053e-6, 5.4438e-5),
      evaluations = c(4L, 1L, 3L, 4L)
    )
  ),
  table_5 = list(
    objective = f5,
    armijo_constant = 0.001,
    curvature_constant = 0.001,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(0.074271, 0.071748, 0.075308, 0.073117),
      value = c(0.99138, 0.99139, 0.99138, 0.99138),
      gradient = c(1.7174e-5, -6.1309e-4, 2.5835e-4, -2.6323e-4),
      evaluations = c(7L, 3L, 6L, 8L)
    )
  ),
  table_6 = list(
    objective = f6,
    armijo_constant = 0.001,
    curvature_constant = 0.001,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(0.92579, 0.92706, 0.92593, 0.92918),
      value = c(0.99138, 0.99139, 0.99138, 0.99139),
      gradient = c(-3.5017e-6, 3.0799e-4, 3.1671e-5, 8.6261e-4),
      evaluations = c(8L, 11L, 7L, 10L)
    )
  )
)

# These rounded outputs characterize the supported Schmidt Wolfe search. They
# are numerical oracles, not contracts for any private policy or helper shape.
test_that("Schmidt Wolfe retains its benchmark output oracles", {
  for (table_name in names(schmidt_wolfe_table_oracles)) {
    oracle <- schmidt_wolfe_table_oracles[[table_name]]
    for (row in seq_len(nrow(oracle$expected))) {
      expected <- oracle$expected[row, ]
      expected$alpha <- expected$parameter
      info <- paste(table_name, "initial alpha", expected$initial_alpha)
      result <- run_schmidt_wolfe_oracle(
        objective = oracle$objective,
        initial_parameter = 0,
        initial_alpha = expected$initial_alpha,
        armijo_constant = oracle$armijo_constant,
        curvature_constant = oracle$curvature_constant
      )

      expect_schmidt_wolfe_oracle(result, expected, info)
    }
  }
})

test_that("Schmidt Wolfe retains its modified-function output oracles", {
  cases <- list(
    function_4 = list(
      objective = f4,
      parameter = 0.5,
      value = 0.999,
      gradient = 1.3132e-11,
      alpha = 0.5,
      evaluations = 12L
    ),
    function_5 = list(
      objective = f5,
      parameter = 0.50112,
      value = 0.99464,
      gradient = 0.0087536,
      alpha = 0.49888,
      evaluations = 14L
    ),
    function_6 = list(
      objective = f6,
      parameter = 0.82848,
      value = 0.99188,
      gradient = -0.0072578,
      alpha = 0.17152,
      evaluations = 14L
    )
  )

  for (name in names(cases)) {
    expected <- cases[[name]]
    result <- run_schmidt_wolfe_oracle(
      objective = expected$objective,
      initial_parameter = 1,
      initial_alpha = 1,
      armijo_constant = 0.1,
      curvature_constant = 0.9
    )

    expect_schmidt_wolfe_oracle(result, expected, name)
  }
})

test_that("Schmidt Wolfe recovers from an undefined cubic expansion", {
  objective <- function(alpha) {
    1 - alpha + 2.5 * alpha^2 - 2 * alpha^3
  }
  gradient <- function(alpha) {
    -1 + 5 * alpha - 6 * alpha^2
  }
  evaluated_alphas <- numeric()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated_alphas <<- c(evaluated_alphas, alpha)
    list(
      alpha = alpha,
      f = objective(alpha),
      df = gradient(alpha),
      d = gradient(alpha),
      par = alpha
    )
  }
  initial_step <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)

  result <- new_schmidt_wolfe_search(
    armijo_constant = 1e-4,
    curvature_constant = 0.9,
    max_evaluations = 2
  )(
    phi = phi,
    step0 = initial_step,
    alpha = 1,
    pm = 1
  )

  expect_equal(evaluated_alphas, c(1, 5.505))
  expect_equal(result$step$alpha, 5.505)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})
