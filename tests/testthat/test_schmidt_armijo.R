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
    }
  }
})

test_that("public backtracking reaches an exactly representable minimizer", {
  initial_parameter <- 2^-32
  fg <- list(
    fn = function(x) 16 * x^2,
    gr = function(x) 32 * x
  )

  result <- mize(
    par = initial_parameter,
    fg = fg,
    method = "SD",
    line_search = "backtracking",
    step0 = 1 / 8,
    step_next_init = "slope ratio",
    step_down = 0.5,
    ls_max_fn = 20,
    max_iter = 3,
    max_fg = 100,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = 1e-12,
    step_tol = NULL,
    check_conv_every = 1,
    store_progress = TRUE
  )

  final_progress <- result$progress[nrow(result$progress), , drop = FALSE]
  expect_identical(result$par, 0)
  expect_identical(result$terminate$what, "ginf_tol")
  expect_identical(final_progress$alpha_init, 1 / 8)
  expect_identical(final_progress$alpha, 1 / 32)
  expect_identical(final_progress$ls_nf, 3)
  expect_identical(final_progress$ls_outcome, "armijo")
})

test_that("Schmidt Armijo stops before evaluating the starting parameters", {
  # The trial rounds to the starting parameters. A callback that falsely reports
  # improvement exposes any failure to detect this before evaluation.
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 1
  )
  callback_count <- 0L
  evaluate_line <- function(alpha, calc_gradient = TRUE) {
    callback_count <<- callback_count + 1L
    list(
      alpha = alpha,
      value = 0,
      gradient = -1,
      slope = -1,
      parameters = 1
    )
  }

  result <- make_schmidt_armijo_search(max_evaluations = Inf)(
    evaluate_line = evaluate_line,
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = .Machine$double.eps / 4
  )

  expect_identical(callback_count, 0L)
  expect_identical(result$line_point, initial_point)
  expect_identical(result$function_evaluations, 0L)
  expect_identical(result$gradient_evaluations, 0L)
  expect_identical(result$termination_reason, "rounding_stagnation")
  expect_identical(result$outcome, "no_step")
})

test_that("Schmidt Armijo rechecks repeated parameters away from the start", {
  # Both step lengths round to the same parameters, but the shorter step relaxes
  # the Armijo bound enough to accept their objective value.
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 1
  )
  evaluated <- numeric()
  search_direction <- .Machine$double.eps
  evaluate_line <- function(alpha, calc_gradient = FALSE) {
    evaluated <<- c(evaluated, alpha)
    list(
      alpha = alpha,
      value = 0.93,
      gradient = NULL,
      slope = NULL,
      parameters = project_line_parameters(1, alpha, search_direction)
    )
  }

  result <- make_schmidt_armijo_search(
    armijo_constant = 0.1,
    step_down = 0.6,
    max_evaluations = Inf
  )(
    evaluate_line = evaluate_line,
    initial_point = initial_point,
    initial_alpha = 1,
    search_direction = search_direction
  )

  expect_identical(evaluated, c(1, 0.6))
  expect_identical(result$line_point$alpha, 0.6)
  expect_identical(result$termination_reason, "armijo")
  expect_identical(result$outcome, "armijo")
})

test_that("Schmidt Armijo stops when safeguarded alpha cannot contract", {
  initial_alpha <- .Machine$double.xmin * .Machine$double.eps
  search_direction <- 1e300
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -1,
    slope = -1,
    parameters = 0
  )
  evaluate_line <- function(alpha, calc_gradient = TRUE) {
    list(
      alpha = alpha,
      value = 2,
      gradient = -1,
      slope = -1,
      parameters = alpha * search_direction
    )
  }

  for (step_down in list(0.5, NULL)) {
    result <- make_schmidt_armijo_search(
      step_down = step_down,
      max_evaluations = Inf
    )(
      evaluate_line = evaluate_line,
      initial_point = initial_point,
      initial_alpha = initial_alpha,
      search_direction = search_direction
    )

    expect_identical(result$line_point, initial_point)
    expect_identical(result$function_evaluations, 1L)
    expect_identical(result$termination_reason, "rounding_stagnation")
    expect_identical(result$outcome, "no_step")
  }
})
