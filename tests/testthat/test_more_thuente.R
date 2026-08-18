# This test uses input parameters, directions and step sizes from using the MT
# line search with a few steps of the CG solver. The expected values come
# from plugging the input values into Dianne O'Leary's Matlab code (running
# under GNU Octave).

mtls <- function(
  fg,
  x,
  pv = -fg$gr(x) / abs(fg$gr(x)),
  alpha,
  c1,
  c2,
  eps = 1e-6,
  approx_armijo = FALSE,
  strong_curvature = TRUE,
  safeguard_cubic = FALSE
) {
  search <- make_wolfe_line_search(
    more_thuente_core,
    armijo_constant = c1,
    curvature_constant = c2,
    approximation_tolerance = eps,
    approximate_armijo = approx_armijo,
    strong_curvature = strong_curvature,
    method_policy = make_more_thuente_policy(
      safeguard_cubic = safeguard_cubic
    )
  )
  res <- search(
    evaluate_line = make_line_function(x, fg, pv, calc_gradient_default = TRUE),
    alpha,
    initial_point = make_initial_line_point(fg, x, pv),
    search_direction = pv
  )
  res$line_point$parameters <- x + res$line_point$alpha * pv
  res
}

## These tests are designed to reproduce the data in Tables 1-6
## of the More'-Thuente paper. They do so apart from a small number of minor
## differences (what value differs and when are indicated by comments before the
## test). Note that the values here are also reproduced by the Matlab code by
## O'Leary, i.e. where the R result differs from the published data in the
## original More'-Thuente paper, so does the Matlab code.

# Table 1
test_that("Table 1", {
  res11 <- mtls(fg = f1, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.1)
  expect_step(
    res11,
    x = 1.3650,
    value = -0.35333,
    gradient = -0.0091645,
    nfev = 6
  )
  res12 <- mtls(fg = f1, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.1)
  expect_step(
    res12,
    x = 1.4414,
    value = -0.35349,
    gradient = 0.0046645,
    nfev = 3
  )
  res13 <- mtls(fg = f1, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.1)
  expect_step(res13, x = 10, value = -0.098039, gradient = 0.0094195, nfev = 1)
  res14 <- mtls(fg = f1, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.1)
  expect_step(
    res14,
    x = 36.888,
    value = -0.027070,
    gradient = 7.3169e-004,
    nfev = 4
  )
})

# Table 2
test_that("Table 2", {
  # # gradient 7.1e-9
  res21 <- mtls(fg = f2, x = 0, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(
    res21,
    x = 1.5960,
    value = -2.6214,
    gradient = 3.8113e-009,
    nfev = 12
  )
  # gradient 10e-10 could be a typo in the paper and should be 1.0e-10?
  res22 <- mtls(fg = f2, x = 0, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(
    res22,
    x = 1.5960,
    value = -2.6214,
    gradient = 1.0106e-010,
    nfev = 8
  )
  res23 <- mtls(fg = f2, x = 0, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(
    res23,
    x = 1.5960,
    value = -2.6214,
    gradient = -4.9725e-009,
    nfev = 8
  )
  res24 <- mtls(fg = f2, x = 0, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(
    res24,
    x = 1.5960,
    value = -2.6214,
    gradient = -2.3091e-008,
    nfev = 11
  )
})

# Table 3
test_that("Table 3", {
  res31 <- mtls(fg = f3, x = 0, alpha = 1e-3, c1 = 0.1, c2 = 0.1)
  expect_step(
    res31,
    x = 1.0,
    value = -0.011160,
    gradient = -5.1440e-005,
    nfev = 12
  )
  res32 <- mtls(fg = f3, x = 0, alpha = 1e-1, c1 = 0.1, c2 = 0.1)
  expect_step(
    res32,
    x = 1.0,
    value = -0.011160,
    gradient = -1.9224e-004,
    nfev = 12
  )
  res33 <- mtls(fg = f3, x = 0, alpha = 1e1, c1 = 0.1, c2 = 0.1)
  expect_step(
    res33,
    x = 1.0,
    value = -0.011160,
    gradient = -1.9892e-006,
    nfev = 10
  )
  res34 <- mtls(fg = f3, x = 0, alpha = 1e3, c1 = 0.1, c2 = 0.1)
  expect_step(
    res34,
    x = 1.0,
    value = -0.011160,
    gradient = -1.5789e-005,
    nfev = 13
  )
})

# Table 4
test_that("Table 4", {
  # alpha = 0.08
  res41 <- mtls(fg = f4, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(
    res41,
    x = 0.085,
    value = 0.99901,
    gradient = -6.8531e-005,
    nfev = 4
  )
  res42 <- mtls(fg = f4, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(
    res42,
    x = 0.1,
    value = 0.99901,
    gradient = -4.9330e-005,
    nfev = 1
  )
  res43 <- mtls(fg = f4, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(
    res43,
    x = 0.34910,
    value = 0.999,
    gradient = -2.9195e-006,
    nfev = 3
  )
  res44 <- mtls(fg = f4, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(
    res44,
    x = 0.8294,
    value = 0.999,
    gradient = 1.6436e-005,
    nfev = 4
  )
})

# Table 5
test_that("Table 5", {
  res51 <- mtls(fg = f5, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(
    res51,
    x = 0.075011,
    value = 0.99138,
    gradient = 1.9025e-004,
    nfev = 6
  )
  res52 <- mtls(fg = f5, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(
    res52,
    x = 0.07751,
    value = 0.99139,
    gradient = 7.3935e-004,
    nfev = 3
  )
  res53 <- mtls(fg = f5, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(
    res53,
    x = 0.073142,
    value = 0.99138,
    gradient = -2.5691e-004,
    nfev = 7
  )
  res54 <- mtls(fg = f5, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(
    res54,
    x = 0.076159,
    value = 0.99139,
    gradient = 4.4913e-004,
    nfev = 8
  )
})

# Table 6
test_that("Table 6", {
  res61 <- mtls(fg = f6, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.001)
  expect_step(
    res61,
    x = 0.9279,
    value = 0.99139,
    gradient = 5.2203e-004,
    nfev = 13
  )
  res62 <- mtls(fg = f6, x = 0, alpha = 1e-1, c1 = 0.001, c2 = 0.001)
  expect_step(
    res62,
    x = 0.92615,
    value = 0.99138,
    gradient = 8.3588e-005,
    nfev = 11
  )
  res63 <- mtls(fg = f6, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
  expect_step(
    res63,
    x = 0.92478,
    value = 0.99138,
    gradient = -2.3788e-004,
    nfev = 8
  )
  res64 <- mtls(fg = f6, x = 0, alpha = 1e3, c1 = 0.001, c2 = 0.001)
  expect_step(
    res64,
    x = 0.92440,
    value = 0.99139,
    gradient = -3.2498e-004,
    nfev = 11
  )
})

# Test line search modification in
# Xie, D., & Schlick, T. (2002).
# A more lenient stopping rule for line search algorithms.
# Optimization Methods and Software, 17(4), 683-700.
test_that("Safeguard Cubic", {
  # Only test examples that give different results
  res32c <- mtls(
    fg = f3,
    x = 0,
    alpha = 1e-1,
    c1 = 0.1,
    c2 = 0.1,
    safeguard_cubic = TRUE
  )
  expect_step(
    res32c,
    x = 1.0,
    value = -0.011160,
    gradient = -1.5842e-10,
    nfev = 13
  )
  res64c <- mtls(
    fg = f6,
    x = 0,
    alpha = 1e3,
    c1 = 0.001,
    c2 = 0.001,
    safeguard_cubic = TRUE
  )
  expect_step(
    res64c,
    x = 0.92525,
    value = 0.99138,
    gradient = -1.2989e-4,
    nfev = 10
  )
})


# The above tests don't enter the code path where the function is modified much. The
# tests below do exercise that part.
test_that("Function modification", {
  res4m <- mtls(fg = f4, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(
    res4m,
    x = 0.99615,
    value = 0.99913,
    gradient = 0.032049,
    alpha = 0.0038522,
    nfev = 6
  )
  res5m <- mtls(fg = f5, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(
    res5m,
    x = 0.99599,
    value = 0.99914,
    gradient = 0.038284,
    alpha = 0.0040126,
    nfev = 6
  )
  res6m <- mtls(fg = f6, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  expect_step(
    res6m,
    x = 0.95655,
    value = 0.99157,
    gradient = 0.016504,
    alpha = 0.043447,
    nfev = 4
  )
})

test_that("More-Thuente interval updates reject invalid state without mutation", {
  make_interval_point <- function(alpha, f, d) {
    list(alpha = alpha, value = f, slope = d, gradient = d)
  }
  policy <- make_more_thuente_policy(alpha_max = 4)

  invalid_cases <- list(
    trial_outside_bracket = list(
      best_endpoint = make_interval_point(1, 0, -1),
      other_endpoint = make_interval_point(3, 1, 1),
      trial_point = make_interval_point(4, 2, 0.5),
      is_bracketed = TRUE
    ),
    zero_best_slope = list(
      best_endpoint = make_interval_point(0, 0, 0),
      other_endpoint = make_interval_point(0, 0, 0),
      trial_point = make_interval_point(1, -1, -1),
      is_bracketed = FALSE
    )
  )

  # This private transition owns a safety guard that normal search invariants
  # cannot exercise without first constructing an invalid state.
  for (case_name in names(invalid_cases)) {
    case <- invalid_cases[[case_name]]
    state <- initialize_more_thuente_search_state(
      case$best_endpoint,
      case$trial_point$alpha,
      policy
    )
    state$best_endpoint <- case$best_endpoint
    state$other_endpoint <- case$other_endpoint
    state$trial_point <- case$trial_point
    state$is_bracketed <- case$is_bracketed
    state$trial_lower_bound <- 0
    state$trial_upper_bound <- 4
    result <- tryCatch(
      update_more_thuente_interval(state, policy),
      error = identity
    )

    expect_false(inherits(result, "error"), info = case_name)
    if (!inherits(result, "error")) {
      expect_identical(result$state, state, info = case_name)
      expect_identical(result$classification, "invalid", info = case_name)
    }
  }
})

test_that("More-Thuente interval updates cover all four mathematical cases", {
  make_interval_point <- function(alpha, f, d) {
    list(alpha = alpha, value = f, slope = d, gradient = d)
  }
  policy <- make_more_thuente_policy(alpha_max = 4)
  initial <- make_interval_point(0, 0, -1)
  cases <- list(
    higher_value = list(
      trial = make_interval_point(1, 1, -0.5),
      classification = "higher_trial_value",
      is_bracketed = TRUE,
      best_endpoint = initial,
      other_endpoint = make_interval_point(1, 1, -0.5),
      next_alpha = 0.10056217060402
    ),
    opposite_slope = list(
      trial = make_interval_point(1, -0.5, 0.5),
      classification = "lower_value_opposite_slope",
      is_bracketed = TRUE,
      best_endpoint = make_interval_point(1, -0.5, 0.5),
      other_endpoint = initial,
      next_alpha = 2 / 3
    ),
    reduced_slope_magnitude = list(
      trial = make_interval_point(1, -0.5, -0.25),
      classification = "lower_value_reduced_slope_magnitude",
      is_bracketed = FALSE,
      best_endpoint = make_interval_point(1, -0.5, -0.25),
      other_endpoint = initial,
      next_alpha = 4
    ),
    unreduced_slope_magnitude = list(
      trial = make_interval_point(1, -0.5, -2),
      classification = "lower_value_unreduced_slope_magnitude",
      is_bracketed = FALSE,
      best_endpoint = make_interval_point(1, -0.5, -2),
      other_endpoint = initial,
      next_alpha = 4
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    state <- initialize_more_thuente_search_state(
      initial,
      case$trial$alpha,
      policy
    )
    state$trial_point <- case$trial
    state$trial_lower_bound <- 0
    state$trial_upper_bound <- 4
    result <- update_more_thuente_interval(state, policy)

    expect_identical(
      result$classification,
      case$classification,
      info = case_name
    )
    expect_identical(
      result$state$is_bracketed,
      case$is_bracketed,
      info = case_name
    )
    expect_equal(
      result$state$best_endpoint,
      case$best_endpoint,
      info = case_name
    )
    expect_equal(
      result$state$other_endpoint,
      case$other_endpoint,
      info = case_name
    )
    expect_equal(
      result$state$trial_point$alpha,
      case$next_alpha,
      info = case_name
    )
    expect_equal(
      result$state$trial_point[c("value", "slope", "gradient")],
      case$trial[c("value", "slope", "gradient")],
      info = case_name
    )
  }
})

test_that("More-Thuente termination guard reports named reasons", {
  initial_point <- list(alpha = 0, value = 1, slope = -1, gradient = -1)
  cases <- list(
    alpha_min = list(
      expected = "alpha_min",
      point = list(alpha = 0, value = 2, slope = -1, gradient = -1),
      is_bracketed = FALSE,
      interval_update_case = "initial",
      trial_lower_bound = 0,
      trial_upper_bound = 4,
      alpha_min = 0,
      alpha_max = 10,
      evaluation_count = 1,
      max_evaluations = 10,
      relative_interval_tolerance = 1e-8
    ),
    alpha_max = list(
      expected = "alpha_max",
      point = list(alpha = 10, value = 0.5, slope = -2, gradient = -2),
      is_bracketed = FALSE,
      interval_update_case = "initial",
      trial_lower_bound = 0,
      trial_upper_bound = 10,
      alpha_min = 0,
      alpha_max = 10,
      evaluation_count = 1,
      max_evaluations = 10,
      relative_interval_tolerance = 1e-8
    ),
    rounding = list(
      expected = "rounding_stagnation",
      point = list(alpha = 2, value = 2, slope = -1, gradient = -1),
      is_bracketed = FALSE,
      interval_update_case = "invalid",
      trial_lower_bound = 0,
      trial_upper_bound = 4,
      alpha_min = 0,
      alpha_max = 10,
      evaluation_count = 1,
      max_evaluations = 10,
      relative_interval_tolerance = 1e-8
    ),
    narrow_bracket = list(
      expected = "relative_interval_tolerance",
      point = list(alpha = 1 + 5e-13, value = 2, slope = -1, gradient = -1),
      is_bracketed = TRUE,
      interval_update_case = "initial",
      trial_lower_bound = 1,
      trial_upper_bound = 1 + 1e-12,
      alpha_min = 0,
      alpha_max = 10,
      evaluation_count = 1,
      max_evaluations = 10,
      relative_interval_tolerance = 1e-6
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    policy <- make_more_thuente_policy(
      relative_interval_tolerance = case$relative_interval_tolerance,
      alpha_min = case$alpha_min,
      alpha_max = case$alpha_max
    )
    state <- initialize_more_thuente_search_state(
      initial_point,
      case$point$alpha,
      policy
    )
    state$trial_point <- case$point
    state$is_bracketed <- case$is_bracketed
    state$interval_update_case <- case$interval_update_case
    state$trial_lower_bound <- case$trial_lower_bound
    state$trial_upper_bound <- case$trial_upper_bound
    reason <- classify_more_thuente_termination(
      state = state,
      initial_point = initial_point,
      condition_policy = make_line_condition_policy(1e-4, 0.9),
      policy = policy,
      evaluation_count = case$evaluation_count,
      max_evaluations = case$max_evaluations
    )

    expect_identical(reason, case$expected, info = case_name)
  }
})

test_that("More-Thuente policy owns named algorithm defaults", {
  policy <- make_more_thuente_policy()

  expect_equal(policy$relative_interval_tolerance, .Machine$double.eps)
  expect_equal(policy$contraction_factor, 0.66)
  expect_equal(policy$expansion_factor, 4)
  expect_equal(policy$alpha_min, 0)
  expect_equal(policy$alpha_max, Inf)
  expect_false(policy$safeguard_cubic)
  expect_equal(policy$cubic_interior_fraction, 0.001)

  state <- initialize_more_thuente_search_state(
    list(alpha = 0, value = 1, slope = -1, gradient = -1),
    initial_alpha = 1,
    policy = policy
  )
  expect_named(
    state,
    c(
      "best_endpoint",
      "other_endpoint",
      "trial_point",
      "is_bracketed",
      "modified_function_stage",
      "current_interval_width",
      "previous_interval_width",
      "trial_lower_bound",
      "trial_upper_bound",
      "interval_update_case",
      "termination_reason"
    )
  )
})

test_that("More-Thuente termination reasons retain their precedence", {
  initial <- list(alpha = 0, value = 1, slope = -1, gradient = -1)
  policy <- make_more_thuente_policy(
    relative_interval_tolerance = 1,
    alpha_max = 10
  )
  conditions <- make_line_condition_policy(1e-4, 0.9)
  state <- initialize_more_thuente_search_state(initial, 1, policy)
  state$is_bracketed <- TRUE
  state$trial_lower_bound <- 0.5
  state$trial_upper_bound <- 1
  state$interval_update_case <- "invalid"

  state$trial_point <- list(alpha = 1, value = 0, slope = 0, gradient = 0)
  expect_identical(
    classify_more_thuente_termination(
      state,
      initial,
      conditions,
      policy,
      evaluation_count = 1,
      max_evaluations = 1
    ),
    "wolfe"
  )

  state$trial_point <- list(alpha = 1, value = 2, slope = -1, gradient = -1)
  expect_identical(
    classify_more_thuente_termination(
      state,
      initial,
      conditions,
      policy,
      evaluation_count = 1,
      max_evaluations = 1
    ),
    "relative_interval_tolerance"
  )

  state$is_bracketed <- FALSE
  state$trial_point <- list(alpha = 0, value = 2, slope = -1, gradient = -1)
  expect_identical(
    classify_more_thuente_termination(
      state,
      initial,
      conditions,
      policy,
      evaluation_count = 1,
      max_evaluations = 1
    ),
    "budget_exhausted"
  )
})

test_that("More-Thuente termination reasons select the same endpoint roles", {
  policy <- make_more_thuente_policy()
  best <- list(alpha = 0.5, value = 0.25, slope = -0.5, gradient = -0.5)
  trial <- list(alpha = 1, value = 0, slope = 0, gradient = 0)
  state <- initialize_more_thuente_search_state(best, trial$alpha, policy)
  state$best_endpoint <- best
  state$trial_point <- trial

  for (reason in c(
    "budget_exhausted",
    "relative_interval_tolerance",
    "rounding_stagnation"
  )) {
    state$termination_reason <- reason
    expect_identical(
      select_more_thuente_candidate(state),
      best,
      info = reason
    )
  }

  for (reason in c("wolfe", "alpha_min", "alpha_max")) {
    state$termination_reason <- reason
    expect_identical(
      select_more_thuente_candidate(state),
      trial,
      info = reason
    )
  }
})
