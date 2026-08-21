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
  initial_point <- make_initial_line_point(
    objective,
    initial_parameter,
    direction
  )
  result <- make_schmidt_wolfe_search(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    max_evaluations = 10000,
    approximation_tolerance = approximation_tolerance,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature
  )(
    evaluate_line = make_line_function(
      initial_parameter,
      objective,
      direction,
      calc_gradient_default = TRUE
    ),
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    search_direction = direction
  )

  result$line_point$parameters <- initial_parameter +
    result$line_point$alpha * direction
  result
}

expect_schmidt_wolfe_oracle <- function(result, expected, info) {
  expect_equal(
    result$line_point$parameters,
    expected$parameter,
    tolerance = 1e-4,
    info = info
  )
  expect_equal(
    result$line_point$value,
    expected$value,
    tolerance = 1e-4,
    info = info
  )
  expect_equal_abs(
    result$line_point$gradient,
    expected$gradient,
    tolerance = 1e-4,
    info = info
  )
  expect_equal(
    result$line_point$alpha,
    expected$alpha,
    tolerance = 1e-4,
    info = info
  )
  expect_equal(result$function_evaluations, expected$evaluations, info = info)
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

# Retain the full rounded reference tables above, but exercise only the smallest
# and largest initial scales here. Public condition tests and focused transition
# tests own the behavior between those characterization endpoints.
test_that("Schmidt Wolfe retains representative benchmark output oracles", {
  for (table_name in names(schmidt_wolfe_table_oracles)) {
    oracle <- schmidt_wolfe_table_oracles[[table_name]]
    representative_rows <- unique(c(1L, nrow(oracle$expected)))
    for (row in representative_rows) {
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
      value = objective(alpha),
      gradient = gradient(alpha),
      slope = gradient(alpha),
      parameters = alpha
    )
  }
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )

  result <- make_schmidt_wolfe_search(
    armijo_constant = 1e-4,
    curvature_constant = 0.9,
    max_evaluations = 2
  )(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = 1
  )

  expect_equal(evaluated_alphas, c(1, 5.505))
  expect_equal(result$line_point$alpha, 5.505)
  expect_equal(result$function_evaluations, 2)
  expect_equal(result$gradient_evaluations, 2)
})

test_that("Schmidt Wolfe zooms when expansion stops improving", {
  evaluated_steps <- list()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluation <- length(evaluated_steps) + 1L
    step <- if (evaluation <= 2L) {
      list(
        alpha = alpha,
        value = 0,
        gradient = -1,
        slope = -1,
        parameters = alpha
      )
    } else {
      list(
        alpha = alpha,
        value = -0.1,
        gradient = 0,
        slope = 0,
        parameters = alpha
      )
    }
    evaluated_steps[[evaluation]] <<- step
    step
  }
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )

  result <- make_schmidt_wolfe_search(
    armijo_constant = 0.05,
    curvature_constant = 0.1,
    max_evaluations = 4
  )(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = 1
  )

  evaluated_alphas <- vapply(evaluated_steps, `[[`, numeric(1L), "alpha")
  expect_equal(evaluated_alphas[1:2], c(1, 5.505))
  expect_gt(evaluated_alphas[[3L]], evaluated_alphas[[1L]])
  expect_lt(evaluated_alphas[[3L]], evaluated_alphas[[2L]])
  expect_equal(result$line_point, evaluated_steps[[3L]])
  expect_identical(result$function_evaluations, 3L)
  expect_identical(result$gradient_evaluations, 3L)
})

test_that("Schmidt Wolfe disables parameter tolerance without a direction", {
  evaluated_steps <- list()
  evaluate_line <- function(alpha, calc_gradient = TRUE) {
    evaluation <- length(evaluated_steps) + 1L
    point <- if (evaluation <= 3L) {
      list(
        alpha = alpha,
        value = if (evaluation <= 2L) 0 else -0.1,
        gradient = -1,
        slope = -1,
        parameters = alpha
      )
    } else {
      list(
        alpha = alpha,
        value = -0.2,
        gradient = 0,
        slope = 0,
        parameters = alpha
      )
    }
    evaluated_steps[[evaluation]] <<- point
    point
  }
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )

  result <- make_schmidt_wolfe_search(
    armijo_constant = 0.05,
    curvature_constant = 0.1,
    max_evaluations = 5
  )(
    evaluate_line = evaluate_line,
    initial_point = initial_point,
    initial_alpha = 1
  )

  expect_identical(result$termination_reason, "wolfe")
  expect_identical(result$function_evaluations, 4L)
  expect_identical(result$gradient_evaluations, 4L)
  expect_equal(result$line_point, evaluated_steps[[4L]])
})
