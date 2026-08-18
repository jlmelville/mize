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

test_that("Schmidt brackets when the expanded objective stops decreasing", {
  initial_step <- list(alpha = 0, f = 1, df = -1, d = -1)
  previous_step <- list(alpha = 1, f = 0.5, df = -1, d = -1)
  trial_step <- list(alpha = 2, f = 0.5, df = 0, d = 0)
  expansion_state <- list(
    initial_step = initial_step,
    previous_step = previous_step,
    iteration = 1L
  )
  conditions <- new_line_condition_policy(0.05, 0.1)

  result <- classify_schmidt_expansion(
    expansion_state,
    trial_step,
    conditions
  )

  expect_true(conditions$wolfe(initial_step, trial_step))
  expect_false(result$accepted)
  expect_equal(result$bracket, list(previous_step, trial_step))
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
    for (alpha in c(0, -1, Inf, -Inf, NA_real_, NaN)) {
      expect_error(
        search(phi, descent_step, alpha = alpha, pm = 1),
        "initial_alpha",
        info = paste(name, alpha)
      )
    }

    non_descent_step <- list(alpha = 0, f = 1, df = 1, d = 1, par = 0)
    result <- search(phi, non_descent_step, alpha = 1, pm = 1)
    expect_equal(result$step, non_descent_step, info = name)
    expect_identical(result$nfn, 0L, info = name)
    expect_identical(result$ngr, 0L, info = name)
  }
})

test_that("bracket-and-zoom expansion never evaluates an infinite proposal", {
  searches <- list(
    rasmussen = new_rasmussen_wolfe_search(0.05, 0.1),
    schmidt = new_schmidt_wolfe_search(0.05, 0.1)
  )
  initial_step <- list(alpha = 0, f = 0, df = -1, d = -1, par = 0)

  for (name in names(searches)) {
    evaluated_alphas <- numeric()
    phi <- function(alpha, calc_gradient = TRUE) {
      if (!is.finite(alpha)) {
        stop("non-finite alpha reached the callback")
      }
      evaluated_alphas <<- c(evaluated_alphas, alpha)
      list(alpha = alpha, f = -alpha, df = -1, d = -1, par = alpha)
    }

    result <- searches[[name]](
      phi,
      step0 = initial_step,
      alpha = 1e308,
      pm = 1
    )

    expect_equal(
      evaluated_alphas,
      c(1e308, .Machine$double.xmax),
      info = name
    )
    expect_equal(result$step$alpha, .Machine$double.xmax, info = name)
    expect_identical(result$nfn, 2L, info = name)
    expect_identical(result$ngr, 2L, info = name)
  }
})

test_that("Rasmussen zoom stops when no interior alpha is representable", {
  smallest_positive <- .Machine$double.xmin * .Machine$double.eps
  evaluated_alphas <- numeric()
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated_alphas <<- c(evaluated_alphas, alpha)
    if (length(evaluated_alphas) > 1L) {
      stop("zoom repeated without representable progress")
    }
    list(alpha = alpha, f = 2, df = 1, d = 1, par = alpha)
  }
  initial_step <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)

  result <- new_rasmussen_wolfe_search(0.05, 0.1)(
    phi,
    step0 = initial_step,
    alpha = smallest_positive,
    pm = 1
  )

  expect_equal(evaluated_alphas, smallest_positive)
  expect_equal(result$step, initial_step)
  expect_identical(result$nfn, 1L)
  expect_identical(result$ngr, 1L)
})

test_that("Schmidt zoom proposals stay safeguarded inside their bracket", {
  zoom_state <- initialize_schmidt_zoom_state(list(
    list(alpha = 0, f = 1, df = -1, d = -1),
    list(alpha = 1, f = 2, df = 1, d = 1)
  ))
  zoom_state$previous_proposal_not_safely_interior <- TRUE

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
  other_step <- list(alpha = 2, f = 2, df = 1, d = 1)
  zoom_state <- initialize_schmidt_zoom_state(list(other_step, initial_step))
  expect_equal(zoom_state$endpoints, list(other_step, initial_step))

  old_width <- abs(
    zoom_state$endpoints[[2L]]$alpha - zoom_state$endpoints[[1L]]$alpha
  )

  zoom_state <- update_schmidt_zoom(
    zoom_state,
    trial_step,
    initial_step,
    conditions
  )
  new_width <- abs(
    zoom_state$endpoints[[2L]]$alpha - zoom_state$endpoints[[1L]]$alpha
  )
  endpoint_values <- vapply(zoom_state$endpoints, `[[`, numeric(1L), "f")

  expect_lt(new_width, old_width)
  expect_equal(min(endpoint_values), trial_step$f)
})

test_that("Schmidt zoom preserves endpoint identity across value ties", {
  conditions <- new_line_condition_policy(0.1, 0.5)
  initial_step <- list(alpha = 0, f = 1, df = -1, d = -1)
  best_step <- list(alpha = 2, f = 0, df = -1, d = -1)
  tied_trial <- list(alpha = 1, f = 0, df = -1, d = -1)
  zoom_state <- initialize_schmidt_zoom_state(list(initial_step, best_step))

  zoom_state <- update_schmidt_zoom(
    zoom_state,
    tied_trial,
    initial_step,
    conditions
  )

  expect_equal(zoom_state$endpoints, list(tied_trial, best_step))
  expect_identical(
    which.min(vapply(zoom_state$endpoints, `[[`, numeric(1L), "f")),
    1L
  )
})
