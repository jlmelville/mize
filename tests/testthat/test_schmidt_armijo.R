# These values were produced by the original minFunc implementation for the two
# Armijo policies that mize exposes. They are output oracles, not an internal
# API contract: the R implementation is free to use different names and
# structure.

schmidt_armijo_oracle_cases <- list(
  `function 1` = list(
    fg = f1,
    armijo_constant = 0.001,
    cubic = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 1e1, 37.054),
      value = c(-5e-4, -0.049751, -0.098039, -0.026948),
      gradient = c(-0.5, -0.49256, 0.0094195, 7.2516e-4),
      evaluations = c(1L, 1L, 1L, 4L)
    ),
    fixed = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 1e1, 31.25),
      value = c(-5e-4, -0.049751, -0.098039, -0.031935),
      evaluations = c(1L, 1L, 1L, 6L)
    )
  ),
  `function 2` = list(
    fg = f2,
    armijo_constant = 0.1,
    cubic = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 1.3546, 1.8850),
      value = c(-1.2469e-9, -2.2181e-4, -2.1852, -1.4135),
      gradient = c(-9.9688e-7, -0.008414, -3.0268, 9.7396),
      evaluations = c(1L, 1L, 4L, 9L)
    ),
    fixed = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 1.25, 1.9531),
      value = c(-1.2469e-9, -2.2181e-4, -1.8447, -0.62904),
      evaluations = c(1L, 1L, 4L, 10L)
    )
  ),
  `function 3` = list(
    fg = f3,
    armijo_constant = 0.1,
    cubic = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.020790, 0.015942),
      value = c(0.99999, 0.89747, 0.99466, 0.99745),
      gradient = c(-0.011857, -0.022189, -0.71010, -0.44579),
      evaluations = c(1L, 1L, 2L, 3L)
    ),
    fixed = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 1.25, 1.9531),
      value = c(0.99999, 0.89747, 0.26493, 0.95744),
      evaluations = c(1L, 1L, 4L, 10L)
    )
  ),
  `function 4` = list(
    fg = f4,
    armijo_constant = 0.001,
    cubic = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.34963, 0.8294),
      value = c(0.99941, 0.99901, 0.99900, 0.99900),
      gradient = c(-0.29260, -4.9330e-5, -2.9053e-6, 1.6436e-5),
      evaluations = c(1L, 1L, 3L, 4L)
    ),
    fixed = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.625, 0.97656),
      value = c(0.99941, 0.99901, 0.99900, 0.99902),
      evaluations = c(1L, 1L, 5L, 11L)
    )
  ),
  `function 5` = list(
    fg = f5,
    armijo_constant = 0.001,
    cubic = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.33678, 0.82237),
      value = c(0.99910, 0.99144, 0.99321, 0.99747),
      gradient = c(-0.89065, 0.0039933, 0.0085115, 0.0088923),
      evaluations = c(1L, 1L, 3L, 4L)
    ),
    fixed = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.625, 0.97656),
      value = c(0.99910, 0.99144, 0.99573, 0.99886),
      evaluations = c(1L, 1L, 5L, 11L)
    )
  ),
  `function 6` = list(
    fg = f6,
    armijo_constant = 0.001,
    cubic = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.51544, 0.83698),
      value = c(0.99945, 0.99817, 0.99449, 0.99182),
      gradient = c(-0.29888, -0.0089383, -0.0087397, -0.0070770),
      evaluations = c(1L, 1L, 3L, 4L)
    ),
    fixed = data.frame(
      initial_alpha = c(1e-3, 1e-1, 1e1, 1e3),
      selected_alpha = c(1e-3, 1e-1, 0.625, 0.97656),
      value = c(0.99945, 0.99817, 0.99354, 0.99230),
      evaluations = c(1L, 1L, 5L, 11L)
    )
  )
)

run_schmidt_armijo_oracle_case <- function(
  fg,
  initial_alpha,
  armijo_constant,
  step_down
) {
  initial_parameters <- 0
  direction <- -fg$gr(initial_parameters) / abs(fg$gr(initial_parameters))
  initial_point <- make_initial_line_point(fg, initial_parameters, direction)
  search <- make_schmidt_armijo_search(
    armijo_constant = armijo_constant,
    step_down = step_down
  )

  search(
    evaluate_line = make_line_function(
      initial_parameters,
      fg,
      direction,
      calc_gradient_default = TRUE
    ),
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    search_direction = direction
  )
}

test_that("supported Schmidt cubic Armijo outputs match their oracle", {
  for (case_name in names(schmidt_armijo_oracle_cases)) {
    case <- schmidt_armijo_oracle_cases[[case_name]]
    for (row in seq_len(nrow(case$cubic))) {
      expected <- case$cubic[row, ]
      result <- run_schmidt_armijo_oracle_case(
        fg = case$fg,
        initial_alpha = expected$initial_alpha,
        armijo_constant = case$armijo_constant,
        step_down = NULL
      )
      info <- paste(case_name, "row", row)

      expect_equal(
        result$line_point$alpha,
        expected$selected_alpha,
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
        result$function_evaluations,
        expected$evaluations,
        info = info
      )
      expect_equal(
        result$gradient_evaluations,
        expected$evaluations,
        info = info
      )
      expect_true(result$gradient_is_current, info = info)
    }
  }
})

test_that("supported Schmidt fixed Armijo outputs match their oracle", {
  for (case_name in names(schmidt_armijo_oracle_cases)) {
    case <- schmidt_armijo_oracle_cases[[case_name]]
    for (row in seq_len(nrow(case$fixed))) {
      expected <- case$fixed[row, ]
      result <- run_schmidt_armijo_oracle_case(
        fg = case$fg,
        initial_alpha = expected$initial_alpha,
        armijo_constant = case$armijo_constant,
        step_down = 0.5
      )
      info <- paste(case_name, "row", row)

      expect_equal(
        result$line_point$alpha,
        expected$selected_alpha,
        tolerance = 1e-4,
        info = info
      )
      expect_equal(
        result$line_point$value,
        expected$value,
        tolerance = 1e-4,
        info = info
      )
      expect_null(result$line_point$gradient, info = info)
      expect_equal(
        result$function_evaluations,
        expected$evaluations,
        info = info
      )
      expect_equal(result$gradient_evaluations, 0, info = info)
      expect_false(result$gradient_is_current, info = info)
    }
  }
})
