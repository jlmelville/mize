test_that("bracket-and-zoom policies expose descriptive supported controls", {
  rasmussen <- new_rasmussen_bracket_zoom_policy()
  schmidt <- new_schmidt_bracket_zoom_policy()

  expect_equal(rasmussen$expansion_factor, 3)
  expect_equal(rasmussen$interior_fraction, 0.1)
  expect_equal(rasmussen$relative_interval_tolerance, 1e-6)

  expect_equal(schmidt$expansion_factor, 10)
  expect_equal(schmidt$minimum_expansion_fraction, 0.01)
  expect_equal(schmidt$interior_fraction, 0.1)
  expect_equal(schmidt$progress_tolerance, 1e-9)

  previous_point <- list(alpha = 2, f = 1, df = -1, d = -1)
  expect_equal(
    rasmussen$expansion_recovery_lower_bound(previous_point, 1L),
    0
  )
  expect_equal(
    schmidt$expansion_recovery_lower_bound(previous_point, 1L),
    2
  )

  expect_error(
    new_rasmussen_bracket_zoom_policy(expansion_factor = 1),
    "expansion_factor"
  )
  expect_error(
    new_schmidt_bracket_zoom_policy(interior_fraction = 0.5),
    "interior_fraction"
  )
})

test_that("Schmidt zoom proposals stay safeguarded inside their bracket", {
  policy <- new_schmidt_bracket_zoom_policy()
  initial_point <- list(alpha = 0, f = 1, df = -1, d = -1)
  other_point <- list(alpha = 1, f = 2, df = 1, d = 1)
  state <- policy$initialize_zoom(
    initial_point,
    list(initial_point, other_point),
    other_point
  )
  state$insufficient_progress <- TRUE
  state <- policy$prepare_zoom(
    state,
    other_point,
    initial_point,
    new_line_condition_policy(0.1, 0.5)
  )
  proposal <- policy$propose_zoom(state, initial_point)

  expect_gte(proposal$alpha, 0.1)
  expect_lte(proposal$alpha, 0.9)
})

test_that("Schmidt zoom updates preserve the lower-value role and contract", {
  policy <- new_schmidt_bracket_zoom_policy()
  conditions <- new_line_condition_policy(0.1, 0.5)
  initial_point <- list(alpha = 0, f = 1, df = -1, d = -1)
  upper_point <- list(alpha = 2, f = 2, df = 1, d = 1)
  trial_point <- list(alpha = 1, f = 0.5, df = -1, d = -1)
  state <- policy$initialize_zoom(
    initial_point,
    list(initial_point, upper_point),
    upper_point
  )
  state <- policy$prepare_zoom(
    state,
    upper_point,
    initial_point,
    conditions
  )
  old_width <- abs(state$second_point$alpha - state$first_point$alpha)
  state <- policy$update_zoom(
    state,
    trial_point,
    initial_point,
    conditions
  )
  new_width <- abs(state$second_point$alpha - state$first_point$alpha)
  values <- c(state$first_point$f, state$second_point$f)

  expect_lt(new_width, old_width)
  expect_equal(min(values), trial_point$f)
  expect_true(all(
    c(
      state$first_point$alpha,
      state$second_point$alpha
    ) >=
      0
  ))
})
