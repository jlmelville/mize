record_default_line_search_trace <- function(factory, initial_alpha = 12.5) {
  trace <- list()
  phi <- function(alpha, calc_gradient = TRUE) {
    point <- list(
      alpha = alpha,
      value = (alpha - 1)^2,
      gradient = 2 * (alpha - 1),
      slope = 2 * (alpha - 1),
      parameters = alpha
    )
    trace[[length(trace) + 1L]] <<- unlist(point[c("alpha", "value", "slope")])
    point
  }
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -2,
    slope = -2,
    parameters = 0
  )

  result <- factory()(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = initial_alpha,
    search_direction = 1
  )

  list(
    trials = do.call(rbind, trace),
    selected = unlist(result$line_point[c("alpha", "value", "slope")]),
    callback_counts = c(
      fn = result$function_evaluations,
      gr = result$gradient_evaluations
    )
  )
}

test_that("default Wolfe search profiles retain their trial traces", {
  cases <- list(
    `more-thuente` = list(
      factory = function() {
        make_wolfe_line_search(
          more_thuente_core,
          armijo_constant = 1e-4,
          curvature_constant = 0.1,
          method_policy = make_more_thuente_policy(
            alpha_max = Inf,
            safeguard_cubic = FALSE
          )
        )
      },
      trials = matrix(
        c(12.5, 132.25, 23, 1, 0, 0),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(NULL, c("alpha", "value", "slope"))
      ),
      callbacks = c(fn = 2, gr = 2)
    ),
    rasmussen = list(
      factory = function() {
        make_rasmussen_wolfe_search(
          armijo_constant = 0.05,
          curvature_constant = 0.1,
          interior_fraction = 0.1,
          expansion_factor = 3,
          relative_interval_tolerance = 1e-6
        )
      },
      trials = matrix(
        c(12.5, 132.25, 23, 1.25, 0.0625, 0.5, 1, 0, 0),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(NULL, c("alpha", "value", "slope"))
      ),
      callbacks = c(fn = 3, gr = 3)
    ),
    schmidt = list(
      factory = function() {
        make_schmidt_wolfe_search(
          armijo_constant = 0.05,
          curvature_constant = 0.1
        )
      },
      trials = matrix(
        c(12.5, 132.25, 23, 1, 0, 0),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(NULL, c("alpha", "value", "slope"))
      ),
      callbacks = c(fn = 2, gr = 2)
    )
  )
  expected_selection <- c(alpha = 1, value = 0, slope = 0)

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    result <- record_default_line_search_trace(case$factory)

    expect_equal(result$trials, case$trials, info = case_name)
    expect_equal(result$selected, expected_selection, info = case_name)
    expect_equal(result$callback_counts, case$callbacks, info = case_name)
  }
})

record_schmidt_armijo_trace <- function(step_down) {
  trace <- list()
  phi <- function(alpha, calc_gradient = TRUE) {
    point <- list(alpha = alpha, value = (alpha - 1)^2, parameters = alpha)
    if (calc_gradient) {
      point$gradient <- 2 * (alpha - 1)
      point$slope <- point$gradient
    }
    trace[[length(trace) + 1L]] <<- c(
      alpha = point$alpha,
      value = point$value,
      slope = if (is.null(point$slope)) NA_real_ else point$slope
    )
    point
  }
  initial_point <- list(
    alpha = 0,
    value = 1,
    gradient = -2,
    slope = -2,
    parameters = 0
  )

  result <- make_schmidt_armijo_search(
    armijo_constant = 0.05,
    step_down = step_down
  )(
    evaluate_line = phi,
    initial_point = initial_point,
    initial_alpha = 12.5,
    search_direction = 1
  )

  list(
    trials = do.call(rbind, trace),
    selected = unlist(result$line_point[c("alpha", "value", "slope")]),
    callback_counts = c(
      fn = result$function_evaluations,
      gr = result$gradient_evaluations
    ),
    gradient_is_current = result$gradient_is_current
  )
}

test_that("supported Schmidt Armijo policies retain their trial traces", {
  cubic <- record_schmidt_armijo_trace(step_down = NULL)
  expect_equal(
    cubic$trials,
    matrix(
      c(12.5, 132.25, 23, 1, 0, 0),
      ncol = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("alpha", "value", "slope"))
    )
  )
  expect_equal(cubic$selected, c(alpha = 1, value = 0, slope = 0))
  expect_equal(cubic$callback_counts, c(fn = 2, gr = 2))
  expect_true(cubic$gradient_is_current)

  fixed <- record_schmidt_armijo_trace(step_down = 0.5)
  expect_equal(
    fixed$trials,
    matrix(
      c(
        12.5,
        132.25,
        NA,
        6.25,
        27.5625,
        NA,
        3.125,
        4.515625,
        NA,
        1.5625,
        0.31640625,
        NA
      ),
      ncol = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("alpha", "value", "slope"))
    )
  )
  expect_equal(
    fixed$selected,
    c(alpha = 1.5625, value = 0.31640625)
  )
  expect_equal(fixed$callback_counts, c(fn = 4, gr = 0))
  expect_false(fixed$gradient_is_current)
})
