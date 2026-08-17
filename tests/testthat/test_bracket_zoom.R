test_that("Schmidt Wolfe validates its numerical safeguards", {
  expect_error(
    new_schmidt_wolfe_policy(expansion_factor = 1),
    "expansion_factor"
  )
  expect_error(
    new_schmidt_wolfe_policy(interior_fraction = 0.5),
    "interior_fraction"
  )
})

test_that("bracket-and-zoom validates entry conditions before callbacks", {
  searches <- list(
    rasmussen = new_rasmussen_wolfe_search(0.05, 0.1),
    schmidt = new_schmidt_wolfe_search(0.05, 0.1)
  )
  phi <- function(alpha, calc_gradient = TRUE) {
    stop("phi should not be called")
  }

  for (name in names(searches)) {
    search <- searches[[name]]
    descent_step <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)
    expect_error(
      search(phi, descent_step, alpha = 0, pm = 1),
      "initial_alpha",
      info = name
    )

    non_descent_step <- list(alpha = 0, f = 1, df = 1, d = 1, par = 0)
    result <- search(phi, non_descent_step, alpha = 1, pm = 1)
    expect_equal(result$step, non_descent_step, info = name)
    expect_identical(result$nfn, 0L, info = name)
    expect_identical(result$ngr, 0L, info = name)
  }
})

test_that("Schmidt zoom proposals stay safeguarded inside their bracket", {
  zoom_state <- list(
    first_step = list(alpha = 0, f = 1, df = -1, d = -1),
    second_step = list(alpha = 1, f = 2, df = 1, d = 1),
    insufficient_progress = TRUE
  )

  proposal <- propose_schmidt_zoom_alpha(
    zoom_state,
    interior_fraction = 0.1
  )

  expect_gte(proposal$alpha, 0.1)
  expect_lte(proposal$alpha, 0.9)
})

test_that("Schmidt zoom updates preserve the lower-value role", {
  conditions <- new_line_condition_policy(0.1, 0.5)
  initial_step <- list(alpha = 0, f = 1, df = -1, d = -1)
  trial_step <- list(alpha = 1, f = 0.5, df = -1, d = -1)
  zoom_state <- list(
    first_step = initial_step,
    second_step = list(alpha = 2, f = 2, df = 1, d = 1),
    insufficient_progress = FALSE
  )
  old_width <- abs(
    zoom_state$second_step$alpha - zoom_state$first_step$alpha
  )

  zoom_state <- update_schmidt_zoom(
    zoom_state,
    trial_step,
    initial_step,
    conditions
  )
  new_width <- abs(
    zoom_state$second_step$alpha - zoom_state$first_step$alpha
  )
  values <- c(zoom_state$first_step$f, zoom_state$second_step$f)

  expect_lt(new_width, old_width)
  expect_equal(min(values), trial_step$f)
})
