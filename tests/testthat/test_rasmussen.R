run_rasmussen_oracle <- function(
  objective,
  initial_parameter,
  direction = -objective$gr(initial_parameter) /
    abs(objective$gr(initial_parameter)),
  initial_alpha,
  armijo_constant,
  curvature_constant,
  relative_interval_tolerance = 1e-6,
  approximation_tolerance = 1e-6,
  approximate_armijo = FALSE,
  strong_curvature = TRUE
) {
  initial_step <- make_step0(objective, initial_parameter, direction)
  search <- new_rasmussen_wolfe_search(
    armijo_constant = armijo_constant,
    curvature_constant = curvature_constant,
    max_evaluations = 10000,
    approximation_tolerance = approximation_tolerance,
    approximate_armijo = approximate_armijo,
    strong_curvature = strong_curvature,
    relative_interval_tolerance = relative_interval_tolerance
  )
  result <- search(
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
  result$step0 <- initial_step
  result
}

expect_rasmussen_oracle <- function(result, expected, info) {
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

rasmussen_table_oracles <- list(
  table_1 = list(
    objective = f1,
    armijo_constant = 0.001,
    curvature_constant = 0.1,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(1.2132, 1.2531, 10, 37.054),
      value = c(-0.34944, -0.35098, -0.098039, -0.026948),
      gradient = c(-0.043803, -0.033723, 0.0094195, 7.2516e-4),
      evaluations = c(9L, 5L, 1L, 4L)
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
      gradient = c(0, -4.2819e-10, -1.2233e-10, -1.2233e-10),
      evaluations = c(55L, 38L, 10L, 12L)
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
      gradient = c(2.3645e-6, -3.8131e-6, -4.1227e-15, -4.1227e-5),
      evaluations = c(18L, 15L, 2L, 4L)
    )
  ),
  table_4 = list(
    objective = f4,
    armijo_constant = 0.001,
    curvature_constant = 0.001,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(0.023344, 0.1, 0.47507, 0.9235),
      value = c(0.99902, 0.99901, 0.999002, 0.99901),
      gradient = c(-9.1485e-4, -4.9330e-5, -4.0043e-7, 8.4666e-5),
      evaluations = c(13L, 1L, 3L, 5L)
    )
  ),
  table_5 = list(
    objective = f5,
    armijo_constant = 0.001,
    curvature_constant = 0.001,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(0.074201, 0.07175, 0.074060, 0.073447),
      value = c(0.99138, 0.99139, 0.99138, 0.99138),
      gradient = c(5.3016e-7, -6.1309e-4, -3.3023e-5, -1.8157e-4),
      evaluations = c(8L, 3L, 7L, 9L)
    )
  ),
  table_6 = list(
    objective = f6,
    armijo_constant = 0.001,
    curvature_constant = 0.001,
    expected = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      parameter = c(0.9296, 0.92662, 0.92966, 0.92436),
      value = c(0.99139, 0.99138, 0.99139, 0.99139),
      gradient = c(9.7500e-4, 1.9698e-4, 9.9413e-4, -3.3353e-4),
      evaluations = c(56L, 30L, 7L, 8L)
    )
  )
)

test_that("Rasmussen validates its supported safeguards", {
  expect_error(
    new_rasmussen_wolfe_search(0.05, 0.1, expansion_factor = 1),
    "expansion_factor"
  )
  expect_error(
    new_rasmussen_wolfe_search(0.05, 0.1, interior_fraction = 0.5),
    "interior_fraction"
  )
  expect_error(
    new_rasmussen_wolfe_search(
      0.05,
      0.1,
      relative_interval_tolerance = -1
    ),
    "relative_interval_tolerance"
  )
  expect_error(
    new_rasmussen_wolfe_search(0.05, 0.1, max_evaluations = 1.5),
    "whole number"
  )
})

test_that("Rasmussen safeguards undefined cubic proposals", {
  first_step <- list(alpha = 0.9493646, f = 1.04033, d = 0.6949782)
  second_step <- list(alpha = 0.1898729, f = 0.8138197, d = 0.03981831)

  expect_warning(
    expansion_alpha <- propose_rasmussen_expansion_cubic_alpha(
      first_step,
      second_step
    ),
    NA
  )
  expect_warning(
    zoom_alpha <- propose_rasmussen_zoom_cubic_alpha(first_step, second_step),
    NA
  )
  expect_true(is.nan(expansion_alpha))
  expect_true(is.nan(zoom_alpha))
  expect_equal(
    safeguard_rasmussen_zoom_alpha(
      zoom_alpha,
      lower_alpha = second_step$alpha,
      upper_alpha = first_step$alpha,
      interior_fraction = 0.1
    ),
    mean(c(first_step$alpha, second_step$alpha))
  )
})

test_that("Rasmussen checks zoom progress without skipping the next proposal", {
  coefficients <- c(
    -1.2529485822,
    -0.2158857423,
    -3.1815714967,
    -0.6465067016,
    0.0001063759
  )
  objective <- list(
    fn = function(alpha) sum(coefficients * alpha^(0:4)),
    gr = function(alpha) sum((1:4) * coefficients[2:5] * alpha^(0:3))
  )
  result <- run_rasmussen_oracle(
    objective = objective,
    initial_parameter = 0,
    initial_alpha = 8.111436,
    armijo_constant = 0.0013907,
    curvature_constant = 0.5281278
  )
  conditions <- new_line_condition_policy(0.0013907, 0.5281278)

  expect_equal(result$nfn, 16)
  expect_true(conditions$wolfe(result$step0, result$step))
})

# These rounded outputs characterize the supported Rasmussen search. They are
# numerical oracles, not contracts for any private policy or helper shape.
test_that("Rasmussen retains its benchmark output oracles", {
  for (table_name in names(rasmussen_table_oracles)) {
    oracle <- rasmussen_table_oracles[[table_name]]
    for (row in seq_len(nrow(oracle$expected))) {
      expected <- oracle$expected[row, ]
      expected$alpha <- expected$parameter
      info <- paste(table_name, "initial alpha", expected$initial_alpha)
      result <- run_rasmussen_oracle(
        objective = oracle$objective,
        initial_parameter = 0,
        initial_alpha = expected$initial_alpha,
        armijo_constant = oracle$armijo_constant,
        curvature_constant = oracle$curvature_constant
      )

      expect_rasmussen_oracle(result, expected, info)
    }
  }
})

test_that("Rasmussen retains its modified-function output oracles", {
  cases <- list(
    function_4 = list(
      objective = f4,
      parameter = 0.99278,
      value = 0.99907,
      gradient = 0.009454,
      alpha = 0.0072168,
      evaluations = 6L
    ),
    function_5 = list(
      objective = f5,
      parameter = 0.99243,
      value = 0.99905,
      gradient = 0.017425,
      alpha = 0.0075707,
      evaluations = 6L
    ),
    function_6 = list(
      objective = f6,
      parameter = 0.936501,
      value = 0.99140,
      gradient = 0.0032111,
      alpha = 0.063499,
      evaluations = 4L
    )
  )

  for (name in names(cases)) {
    expected <- cases[[name]]
    result <- run_rasmussen_oracle(
      objective = expected$objective,
      initial_parameter = 1,
      initial_alpha = 1,
      armijo_constant = 0.1,
      curvature_constant = 0.9
    )

    expect_rasmussen_oracle(result, expected, name)
  }
})
