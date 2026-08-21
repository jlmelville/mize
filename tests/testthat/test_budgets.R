make_budget_witness <- function(
  fn = function(x) sum(x^2),
  gr = function(x) 2 * x
) {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0
  calls$gr <- 0

  list(
    fg = list(
      fn = function(x) {
        calls$fn <- calls$fn + 1
        fn(x)
      },
      gr = function(x) {
        calls$gr <- calls$gr + 1
        gr(x)
      }
    ),
    counts = function() c(fn = calls$fn, gr = calls$gr)
  )
}

make_bold_driver_witness <- function(
  fn = function(x) x^4,
  gr = function(x) 4 * x^3
) {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0
  calls$gr <- 0
  calls$fn_par <- numeric()

  list(
    fg = list(
      fn = function(x) {
        calls$fn <- calls$fn + 1
        calls$fn_par <- c(calls$fn_par, x)
        fn(x)
      },
      gr = function(x) {
        calls$gr <- calls$gr + 1
        gr(x)
      }
    ),
    counts = function() c(fn = calls$fn, gr = calls$gr),
    fn_par = function() calls$fn_par
  )
}

make_schmidt_armijo_witness <- function(
  fn = function(x) x^4,
  gr = function(x) 4 * x^3
) {
  calls <- new.env(parent = emptyenv())
  calls$fn_par <- numeric()
  calls$gr_par <- numeric()

  list(
    fg = list(
      fn = function(x) {
        calls$fn_par <- c(calls$fn_par, x)
        fn(x)
      },
      gr = function(x) {
        calls$gr_par <- c(calls$gr_par, x)
        gr(x)
      }
    ),
    counts = function() c(fn = length(calls$fn_par), gr = length(calls$gr_par)),
    fn_par = function() calls$fn_par,
    gr_par = function() calls$gr_par
  )
}

make_wolfe_overflow_witness <- function() {
  calls <- new.env(parent = emptyenv())
  calls$fn_par <- numeric()
  calls$gr_par <- numeric()

  list(
    fg = list(
      fn = function(x) {
        calls$fn_par <- c(calls$fn_par, x)
        if (all(is.finite(x))) 1 else 0
      },
      gr = function(x) {
        calls$gr_par <- c(calls$gr_par, x)
        if (all(is.finite(x))) -1 else 0
      }
    ),
    counts = function() c(fn = length(calls$fn_par), gr = length(calls$gr_par)),
    fn_par = function() calls$fn_par,
    gr_par = function() calls$gr_par
  )
}

make_wolfe_witness <- function(fn, gr) {
  calls <- new.env(parent = emptyenv())
  calls$fn_par <- numeric()
  calls$gr_par <- numeric()

  list(
    fg = list(
      fn = function(x) {
        calls$fn_par <- c(calls$fn_par, x)
        fn(x)
      },
      gr = function(x) {
        calls$gr_par <- c(calls$gr_par, x)
        gr(x)
      }
    ),
    counts = function() c(fn = length(calls$fn_par), gr = length(calls$gr_par)),
    fn_par = function() calls$fn_par,
    gr_par = function() calls$gr_par
  )
}

make_combined_wolfe_witness <- function(fn, gr) {
  calls <- new.env(parent = emptyenv())
  calls$fn_parameters <- numeric()
  calls$gr_parameters <- numeric()
  calls$combined_parameters <- numeric()

  list(
    fg = list(
      fn = function(x) {
        calls$fn_parameters <- c(calls$fn_parameters, x)
        fn(x)
      },
      gr = function(x) {
        calls$gr_parameters <- c(calls$gr_parameters, x)
        gr(x)
      },
      fg = function(x) {
        calls$combined_parameters <- c(calls$combined_parameters, x)
        list(fn = fn(x), gr = gr(x))
      }
    ),
    counts = function() {
      c(
        fn = length(calls$fn_parameters) + length(calls$combined_parameters),
        gr = length(calls$gr_parameters) + length(calls$combined_parameters)
      )
    },
    fn_par = function() c(calls$fn_parameters, calls$combined_parameters),
    gr_par = function() c(calls$gr_parameters, calls$combined_parameters),
    combined_par = function() calls$combined_parameters
  )
}

run_wolfe_witness <- function(
  witness,
  line_search,
  par = 0,
  step0 = 1,
  c1 = 0.1,
  c2 = 0.5,
  default_ls_max_fn = 1,
  ls_max_fn = default_ls_max_fn,
  ls_max_gr = ls_max_fn,
  ls_max_fg = 2 * ls_max_fn,
  max_fn = Inf,
  max_gr = Inf,
  max_fg = Inf
) {
  opt <- make_mize(
    method = "SD",
    line_search = line_search,
    step0 = step0,
    c1 = c1,
    c2 = c2,
    ls_max_fn = ls_max_fn,
    ls_max_gr = ls_max_gr,
    ls_max_fg = ls_max_fg,
    max_fn = max_fn,
    max_gr = max_gr,
    max_fg = max_fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, par, witness$fg)
  mize_step(opt, par, witness$fg)
}

run_rasmussen_witness <- function(...) {
  run_wolfe_witness(..., line_search = "rasmussen")
}

run_hager_zhang_witness <- function(...) {
  run_wolfe_witness(..., line_search = "hager-zhang")
}

run_schmidt_wolfe_witness <- function(...) {
  run_wolfe_witness(..., line_search = "schmidt", default_ls_max_fn = 2)
}

Ops.ascent_gradient <- function(e1, e2) {
  if (.Generic == "-" && missing(e2)) {
    return(unclass(e1))
  }
  if (.Generic == "*") {
    return(unclass(e1) * e2)
  }
  NextMethod()
}

expect_budget_counts <- function(result, witness) {
  expect_equal(c(fn = result$nf, gr = result$ng), witness$counts())
}

budget_test_args <- function(fg, ...) {
  c(
    list(
      par = 1,
      fg = fg,
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    ),
    list(...)
  )
}

test_that("public zero budgets invoke no callbacks and return no objective", {
  budgets <- c("max_fn", "max_gr", "max_fg")

  for (budget in budgets) {
    witness <- make_budget_witness()
    args <- budget_test_args(witness$fg)
    args[[budget]] <- 0

    result <- do.call(mize, args)

    expect_identical(witness$counts(), c(fn = 0, gr = 0), info = budget)
    expect_budget_counts(result, witness)
    expect_equal(result$terminate$what, budget, info = budget)
    expect_equal(result$par, 1, info = budget)
    expect_false("f" %in% names(result), info = budget)
    expect_false("best_f" %in% names(result), info = budget)
  }
})

test_that("TN zero gradient and combined budgets invoke no callbacks", {
  for (budget in c("max_gr", "max_fg")) {
    witness <- make_budget_witness()
    args <- list(
      par = c(1, -1),
      fg = witness$fg,
      method = "TN",
      line_search = "constant",
      step0 = 1,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )
    args[[budget]] <- 0

    result <- do.call(mize, args)

    expect_identical(witness$counts(), c(fn = 0, gr = 0), info = budget)
    expect_budget_counts(result, witness)
    expect_equal(result$terminate$what, budget, info = budget)
    expect_equal(result$terminate$val, 0, info = budget)
    expect_equal(result$par, c(1, -1), info = budget)
    expect_false("f" %in% names(result), info = budget)
  }
})

test_that("TN does not exceed a one-callback combined budget", {
  witness <- make_budget_witness()

  result <- mize(
    c(1, -1),
    witness$fg,
    method = "TN",
    line_search = "constant",
    step0 = 1,
    max_iter = 1,
    max_fg = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_identical(witness$counts(), c(fn = 0, gr = 1))
  expect_budget_counts(result, witness)
  expect_equal(result$terminate, list(what = "max_fg", val = 1))
  expect_equal(result$par, c(1, -1))
  expect_false("f" %in% names(result))
})

test_that("TN blocks a subsequent finite difference at exact boundaries", {
  rosenbrock_gr <- function(x) {
    c(
      -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
      200 * (x[2] - x[1]^2)
    )
  }

  for (budget in c("max_gr", "max_fg")) {
    witness <- make_budget_witness(gr = rosenbrock_gr)
    args <- list(
      par = c(-1.2, 1),
      fg = witness$fg,
      method = "TN",
      line_search = "constant",
      step0 = 1,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )
    args[[budget]] <- 2

    result <- do.call(mize, args)

    expect_identical(witness$counts(), c(fn = 0, gr = 2), info = budget)
    expect_budget_counts(result, witness)
    expect_equal(
      result$terminate,
      list(what = budget, val = 2),
      info = budget
    )
    # The allowed finite difference completes a usable TN direction, so the
    # constant step is applied before the ordinary lifecycle records the cap.
    expect_equal(
      result$par,
      c(-1.056697, 1.058491),
      tolerance = 1e-6,
      info = budget
    )
    expect_false("f" %in% names(result), info = budget)
  }
})

test_that("TN previous initialization cannot bypass the combined cap", {
  witness <- make_budget_witness()
  par <- c(1, -1)
  opt <- make_mize(
    method = "TN",
    tn_init = "previous",
    line_search = "constant",
    step0 = 0.1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, par, witness$fg)

  first <- mize_step(opt, par, witness$fg)
  expect_false(all(first$opt$stages$gradient_descent$direction$value == 0))

  calls_before <- witness$counts()
  first$opt$convergence$max_fg <- sum(calls_before) + 1
  second <- mize_step(first$opt, first$par, witness$fg)

  expect_identical(
    witness$counts(),
    calls_before + c(fn = 0, gr = 1)
  )
  expect_equal(c(fn = second$nf, gr = second$ng), witness$counts())
  expect_equal(
    second$opt$terminate,
    list(what = "max_fg", val = sum(calls_before) + 1)
  )
  expect_equal(second$par, first$par)
})

test_that("public budgets stop exactly after one callback", {
  cases <- list(
    list(budget = "max_fn", max_iter = 0, expected = c(fn = 1, gr = 0)),
    list(budget = "max_gr", max_iter = 1, expected = c(fn = 0, gr = 1)),
    list(budget = "max_fg", max_iter = 1, expected = c(fn = 0, gr = 1))
  )

  for (case in cases) {
    witness <- make_budget_witness()
    args <- budget_test_args(witness$fg)
    args$max_iter <- case$max_iter
    args[[case$budget]] <- 1

    result <- do.call(mize, args)

    expect_identical(witness$counts(), case$expected, info = case$budget)
    expect_budget_counts(result, witness)
    expect_equal(result$terminate$what, case$budget, info = case$budget)
    if (case$budget == "max_fg") {
      expect_false("f" %in% names(result))
    }
  }
})

test_that("line-search callbacks stop at exact function and combined caps", {
  cases <- list(
    list(budget = "max_fn", expected = c(fn = 2, gr = 2)),
    list(budget = "max_fg", expected = c(fn = 1, gr = 1))
  )

  for (case in cases) {
    witness <- make_budget_witness()
    args <- budget_test_args(witness$fg)
    args$line_search <- "backtracking"
    args$max_iter <- 3
    args[[case$budget]] <- 2

    result <- do.call(mize, args)

    expect_identical(witness$counts(), case$expected, info = case$budget)
    expect_budget_counts(result, witness)
    expect_equal(result$terminate$what, case$budget, info = case$budget)
  }
})

test_that("budget exhaustion before a Bold Driver trial leaves par unchanged", {
  witness <- make_budget_witness()

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "bold",
    max_iter = 3,
    max_fg = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_identical(witness$counts(), c(fn = 0, gr = 1))
  expect_budget_counts(result, witness)
  expect_equal(result$terminate$what, "max_fg")
  expect_equal(result$par, 1)
  expect_false("f" %in% names(result))
})

test_that("Bold Driver rejects trials outside zero and one callback budgets", {
  quartic <- function(x) x^4
  expected_calls <- list(`0` = 1, `1` = c(1, 1))

  for (limit in 0:1) {
    witness <- make_bold_driver_witness()

    result <- mize(
      1,
      witness$fg,
      method = "SD",
      line_search = "bold",
      max_iter = 1,
      ls_max_fn = limit,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )

    expect_equal(witness$fn_par(), expected_calls[[as.character(limit)]])
    expect_budget_counts(result, witness)
    expect_equal(result$par, 1)
    expect_equal(quartic(result$par), quartic(1))
    expect_equal(result$f, quartic(result$par))
  }
})

test_that("Bold Driver charges an uncached starting cost to its local budget", {
  quartic <- function(x) x^4
  witness <- make_bold_driver_witness()

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "bold",
    max_iter = 1,
    ls_max_fn = 3,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1, -3, -1, 1))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 1)
  expect_equal(quartic(result$par), quartic(1))
  expect_equal(result$f, quartic(result$par))
})

test_that("Bold Driver failure does not cache a rejected candidate", {
  quartic <- function(x) x^4
  witness <- make_bold_driver_witness()
  opt <- make_mize(
    method = "SD",
    line_search = "bold",
    ls_max_fn = 3,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, 1, witness$fg)

  result <- mize_step(opt, 1, witness$fg)

  expect_equal(witness$fn_par(), c(1, -3, -1))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 1)
  expect_equal(quartic(result$par), quartic(1))
  expect_false(identical(result$opt$cache$fn_new_iter, 1))
  expect_false("f" %in% names(result))
})

test_that("Bold Driver accepts a strict decrease on the final callback", {
  quartic <- function(x) x^4
  witness <- make_bold_driver_witness()

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "bold",
    max_iter = 1,
    ls_max_fn = 4,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1, -3, -1, 0))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 0)
  expect_lt(quartic(result$par), quartic(1))
  expect_equal(result$f, quartic(result$par))
})

test_that("Bold Driver does not charge a cached starting cost", {
  quartic <- function(x) x^4
  witness <- make_bold_driver_witness()

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "bold",
    max_iter = 1,
    ls_max_fn = 1,
    check_conv_every = 1,
    store_progress = TRUE,
    abs_tol = 0,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1, -3, 1))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 1)
  expect_equal(quartic(result$par), quartic(1))
  expect_equal(result$f, quartic(result$par))
})

test_that("Bold Driver rejects nonfinite candidates", {
  objective <- function(x) {
    if (x == -3) Inf else x^4
  }
  witness <- make_bold_driver_witness(fn = objective)

  result <- suppressMessages(mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "bold",
    max_iter = 1,
    ls_max_fn = 2,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  ))

  expect_equal(witness$fn_par(), c(1, -3, 1))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 1)
  expect_equal(objective(result$par), objective(1))
  expect_equal(result$f, objective(result$par))
})

test_that("Bold Driver non-descent shortcut returns an exact zero step", {
  objective <- function(x) x^2
  witness <- make_bold_driver_witness(
    fn = objective,
    gr = function(x) structure(2 * x, class = "ascent_gradient")
  )

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "bold",
    max_iter = 1,
    ls_max_fn = 0,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), 1)
  expect_budget_counts(result, witness)
  expect_identical(as.numeric(result$par), 1)
  expect_equal(objective(result$par), objective(1))
  expect_equal(result$f, objective(result$par))
})

test_that("fixed Schmidt Armijo budgets only objective evaluations", {
  cases <- list(
    local_gradient = list(ls_max_gr = 0),
    local_combined = list(ls_max_fg = 1),
    global_gradient = list(max_gr = 1),
    global_combined = list(max_fg = 3)
  )

  for (case_name in names(cases)) {
    counts <- c(fn = 0L, gr = 0L)
    fg <- list(
      fn = function(x) {
        counts[["fn"]] <<- counts[["fn"]] + 1L
        x^2
      },
      gr = function(x) {
        counts[["gr"]] <<- counts[["gr"]] + 1L
        2 * x
      }
    )
    arguments <- c(
      list(
        par = 1,
        fg = fg,
        method = "SD",
        line_search = "backtracking",
        step0 = 0.25,
        step_down = 0.5,
        max_iter = 1,
        check_conv_every = NULL,
        abs_tol = NULL,
        rel_tol = NULL,
        grad_tol = NULL,
        ginf_tol = NULL,
        step_tol = NULL
      ),
      cases[[case_name]]
    )

    result <- do.call(mize, arguments)

    expect_identical(counts, c(fn = 2L, gr = 1L), info = case_name)
    expect_equal(result$par, 0.5, info = case_name)
    expect_equal(result$f, 0.25, info = case_name)
  }
})

test_that("Schmidt Armijo rejects one-trial quartic failures", {
  quartic <- function(x) x^4
  cases <- list(
    interpolation = list(step_down = NULL, gr_par = c(1, -3)),
    fixed = list(step_down = 0.5, gr_par = 1)
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    witness <- make_schmidt_armijo_witness()

    result <- mize(
      1,
      witness$fg,
      method = "SD",
      line_search = "backtracking",
      step0 = 1,
      step_down = case$step_down,
      ls_max_fn = 1,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )

    expect_equal(witness$fn_par(), c(1, -3), info = name)
    expect_equal(witness$gr_par(), case$gr_par, info = name)
    expect_budget_counts(result, witness)
    expect_equal(result$par, 1, info = name)
    expect_equal(quartic(result$par), quartic(1), info = name)
    expect_equal(result$f, quartic(result$par), info = name)
  }
})

test_that("Schmidt Armijo rejects an Armijo trial at overflowed parameters", {
  objective <- function(x) if (is.finite(x)) 1 else 0
  gradient <- function(x) if (is.finite(x)) -1 else 0
  witness <- make_schmidt_armijo_witness(objective, gradient)

  result <- mize(
    1e308,
    witness$fg,
    method = "SD",
    line_search = "backtracking",
    step0 = 1e308,
    c1 = 1e-320,
    ls_max_fn = 1,
    max_iter = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1e308, Inf))
  expect_equal(witness$gr_par(), c(1e308, Inf))
  expect_budget_counts(result, witness)
  expect_true(is.finite(result$par))
  expect_equal(result$par, 1e308)
  expect_equal(result$f, objective(result$par))
})

test_that("Wolfe searches reject condition-satisfying overflowed parameters", {
  searches <- c(
    `more-thuente` = "more-thuente",
    rasmussen = "rasmussen",
    schmidt = "schmidt",
    `hager-zhang` = "hager-zhang"
  )

  for (name in names(searches)) {
    witness <- make_wolfe_overflow_witness()
    opt <- make_mize(
      method = "SD",
      line_search = searches[[name]],
      step0 = 1e308,
      c1 = 1e-308,
      c2 = 0.5,
      ls_max_fn = 1,
      ls_max_gr = 1,
      ls_max_fg = 2,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )
    opt <- mize_init(opt, 1e308, witness$fg)

    result <- mize_step(opt, 1e308, witness$fg)

    expect_equal(result$par, 1e308, info = name)
    expect_equal(result$f, 1, info = name)
    expect_equal(result$g, -1, info = name)
    expect_identical(witness$counts(), c(fn = 2L, gr = 2L), info = name)
    expect_budget_counts(result, witness)
    expect_equal(witness$fn_par(), c(1e308, Inf), info = name)
    expect_equal(witness$gr_par(), c(1e308, Inf), info = name)
  }
})

test_that("Wolfe searches retain the current point after initializer overflow", {
  searches <- c(rasmussen = "rasmussen", schmidt = "schmidt")

  for (name in names(searches)) {
    witness <- make_wolfe_witness(
      fn = function(x) {
        if (!all(is.finite(x))) {
          stop("non-finite parameter reached the callback")
        }
        if (x == 0) .Machine$double.xmax else -.Machine$double.xmax
      },
      gr = function(x) {
        if (!all(is.finite(x))) {
          stop("non-finite parameter reached the callback")
        }
        if (x == 0) -2 else -0.5
      }
    )

    result <- mize(
      0,
      witness$fg,
      method = "SD",
      line_search = searches[[name]],
      step0 = 1,
      step_next_init = "quadratic",
      c1 = 0.1,
      c2 = 0.5,
      max_iter = 2,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )

    expect_equal(result$par, 2, info = name)
    expect_equal(result$f, -.Machine$double.xmax, info = name)
    expect_equal(witness$fn_par(), c(0, 2), info = name)
    expect_equal(witness$gr_par(), c(0, 2), info = name)
    expect_identical(witness$counts(), c(fn = 2L, gr = 2L), info = name)
    expect_budget_counts(result, witness)
  }
})

test_that("Rasmussen rejects an unsafe final extrapolation trial", {
  witness <- make_wolfe_witness(
    fn = function(x) x^4,
    gr = function(x) 4 * x^3
  )
  result <- run_rasmussen_witness(witness, par = 1)

  expect_equal(result$par, 1)
  expect_equal(result$f, 1)
  expect_equal(result$g, 4)
  expect_equal(witness$fn_par(), c(1, -3))
  expect_equal(witness$gr_par(), c(1, -3))
  expect_identical(witness$counts(), c(fn = 2L, gr = 2L))
  expect_budget_counts(result, witness)
})

test_that("Rasmussen preserves safe trials at exact exhaustion", {
  cases <- list(
    strict_decrease = list(
      witness = make_wolfe_witness(
        fn = function(x) 1 - x,
        gr = function(x) -1
      ),
      step0 = 1,
      ls_max_fn = 1,
      par = 1,
      f = 0,
      g = -1,
      trace = c(0, 1)
    ),
    wolfe = list(
      witness = make_wolfe_witness(
        fn = function(x) (x - 1)^2,
        gr = function(x) 2 * (x - 1)
      ),
      step0 = 0.5,
      ls_max_fn = 1,
      par = 1,
      f = 0,
      g = 0,
      trace = c(0, 1)
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    result <- run_rasmussen_witness(
      case$witness,
      par = 0,
      step0 = case$step0,
      ls_max_fn = case$ls_max_fn
    )

    expect_equal(result$par, case$par, info = name)
    expect_equal(result$f, case$f, info = name)
    expect_equal(result$g, case$g, info = name)
    expect_equal(case$witness$fn_par(), case$trace, info = name)
    expect_equal(case$witness$gr_par(), case$trace, info = name)
    expect_budget_counts(result, case$witness)
  }
})

test_that("Rasmussen rejects equal and nonfinite exhausted trials", {
  cases <- list(
    equal = list(
      witness = make_wolfe_witness(
        fn = function(x) (x - 1)^2,
        gr = function(x) 2 * (x - 1)
      ),
      g = -2,
      trial = 2
    ),
    nonfinite_derivative = list(
      witness = make_wolfe_witness(
        fn = function(x) 1 - x,
        gr = function(x) if (x == 0) -1 else Inf
      ),
      g = -1,
      trial = 1
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    witness <- case$witness
    result <- run_rasmussen_witness(witness, par = 0)

    expect_equal(result$par, 0, info = name)
    expect_equal(result$f, 1, info = name)
    expect_equal(result$g, case$g, info = name)
    expect_equal(witness$fn_par(), c(0, case$trial), info = name)
    expect_equal(witness$gr_par(), c(0, case$trial), info = name)
    expect_budget_counts(result, witness)
  }
})

test_that("Rasmussen honors each local and global evaluation limit", {
  cases <- list(
    ls_max_fn = list(ls_max_fn = 2),
    ls_max_gr = list(ls_max_gr = 2),
    ls_max_fg = list(ls_max_fg = 4),
    max_fn = list(max_fn = 3),
    max_gr = list(max_gr = 3),
    max_fg = list(max_fg = 6)
  )

  witness_factories <- list(
    separate = make_wolfe_witness,
    combined = make_combined_wolfe_witness
  )

  for (interface in names(witness_factories)) {
    for (name in names(cases)) {
      info <- paste(interface, name)
      witness <- witness_factories[[interface]](
        fn = function(x) if (x == 0) 1 else Inf,
        gr = function(x) -1
      )
      args <- modifyList(
        list(
          witness = witness,
          par = 0,
          ls_max_fn = Inf,
          ls_max_gr = Inf,
          ls_max_fg = Inf
        ),
        cases[[name]]
      )

      result <- do.call(run_rasmussen_witness, args)

      expect_equal(witness$fn_par(), c(0, 1, 0.5), info = info)
      expect_equal(witness$gr_par(), c(0, 1, 0.5), info = info)
      if (interface == "combined") {
        expect_equal(witness$combined_par(), c(1, 0.5), info = info)
      }
      expect_identical(witness$counts(), c(fn = 3L, gr = 3L), info = info)
      expect_budget_counts(result, witness)
      expect_equal(result$par, 0, info = info)
      expect_equal(result$f, 1, info = info)
      expect_equal(result$g, -1, info = info)
    }
  }
})

test_that("Schmidt Wolfe stops before an exhausted Armijo fallback", {
  cases <- list(
    nonfinite_objective = list(
      fn = function(x) if (x == 0) 1 else Inf,
      gr = function(x) -1
    ),
    nonfinite_derivative = list(
      fn = function(x) 1 - x,
      gr = function(x) if (x == 0) -1 else Inf
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    witness <- make_wolfe_witness(case$fn, case$gr)
    result <- run_schmidt_wolfe_witness(
      witness,
      max_fn = 3,
      max_gr = 3,
      max_fg = 6
    )

    expect_equal(witness$fn_par(), c(0, 1, 0.5), info = name)
    expect_equal(witness$gr_par(), c(0, 1, 0.5), info = name)
    expect_identical(witness$counts(), c(fn = 3L, gr = 3L), info = name)
    expect_budget_counts(result, witness)
    expect_equal(result$par, 0, info = name)
    expect_equal(result$f, 1, info = name)
    expect_equal(result$g, -1, info = name)
  }
})

test_that("Schmidt Wolfe honors each local and global evaluation limit", {
  cases <- list(
    ls_max_fn = list(ls_max_fn = 2),
    ls_max_gr = list(ls_max_gr = 2),
    ls_max_fg = list(ls_max_fg = 4),
    max_fn = list(max_fn = 3),
    max_gr = list(max_gr = 3),
    max_fg = list(max_fg = 6)
  )

  for (name in names(cases)) {
    witness <- make_wolfe_witness(
      fn = function(x) if (x == 0) 1 else Inf,
      gr = function(x) -1
    )
    args <- modifyList(
      list(
        witness = witness,
        ls_max_fn = Inf,
        ls_max_gr = Inf,
        ls_max_fg = Inf
      ),
      cases[[name]]
    )
    result <- do.call(run_schmidt_wolfe_witness, args)

    expect_equal(witness$fn_par(), c(0, 1, 0.5), info = name)
    expect_equal(witness$gr_par(), c(0, 1, 0.5), info = name)
    expect_identical(witness$counts(), c(fn = 3L, gr = 3L), info = name)
    expect_budget_counts(result, witness)
    expect_equal(result$par, 0, info = name)
    expect_equal(result$f, 1, info = name)
    expect_equal(result$g, -1, info = name)
  }
})

test_that("Schmidt Wolfe preserves safe final evaluations", {
  cases <- list(
    wolfe = list(
      fn = function(x) (x - 1)^2,
      gr = function(x) 2 * (x - 1),
      ls_max_fn = 2,
      par = 1,
      f = 0,
      g = 0,
      trace = c(0, 2, 1)
    ),
    strict_decrease = list(
      fn = function(x) 1 - x,
      gr = function(x) -1,
      ls_max_fn = 1,
      par = 1,
      f = 0,
      g = -1,
      trace = c(0, 1)
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    witness <- make_wolfe_witness(case$fn, case$gr)
    result <- run_schmidt_wolfe_witness(
      witness,
      ls_max_fn = case$ls_max_fn,
      ls_max_gr = case$ls_max_fn,
      ls_max_fg = 2 * case$ls_max_fn,
      max_fn = case$ls_max_fn + 1,
      max_gr = case$ls_max_fn + 1,
      max_fg = 2 * (case$ls_max_fn + 1)
    )

    expect_equal(witness$fn_par(), case$trace, info = name)
    expect_equal(witness$gr_par(), case$trace, info = name)
    expect_budget_counts(result, witness)
    expect_equal(result$par, case$par, info = name)
    expect_equal(result$f, case$f, info = name)
    expect_equal(result$g, case$g, info = name)
  }
})

test_that("Hager-Zhang rejects an unsafe initial exhausted trial", {
  witness <- make_wolfe_witness(
    fn = function(x) x^4,
    gr = function(x) 4 * x^3
  )
  step_result <- run_hager_zhang_witness(
    witness,
    par = 1,
    max_fn = 2,
    max_gr = 2,
    max_fg = 4
  )
  result <- mize_step(
    step_result$opt,
    step_result$par,
    witness$fg
  )

  expect_equal(witness$fn_par(), c(1, -3))
  expect_equal(witness$gr_par(), c(1, -3))
  expect_identical(witness$counts(), c(fn = 2L, gr = 2L))
  expect_budget_counts(result, witness)
  expect_equal(result$opt$terminate$what, "max_fn")
  expect_equal(result$par, 1)
  expect_equal(step_result$f, 1)
  expect_equal(step_result$g, 4)
})

test_that("Hager-Zhang makes no line-search callback at zero allowance", {
  witness <- make_wolfe_witness(
    fn = function(x) 1 - x,
    gr = function(x) -1
  )
  result <- run_hager_zhang_witness(witness, ls_max_fn = 0)

  expect_equal(witness$fn_par(), 0)
  expect_equal(witness$gr_par(), 0)
  expect_identical(witness$counts(), c(fn = 1L, gr = 1L))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 0)
  expect_equal(result$f, 1)
  expect_equal(result$g, -1)
})

test_that("Hager-Zhang preserves safe initial exhausted trials", {
  cases <- list(
    condition = list(
      witness = make_wolfe_witness(
        fn = function(x) (x - 1)^2,
        gr = function(x) 2 * (x - 1)
      ),
      step0 = 0.5,
      g = 0
    ),
    strict_decrease = list(
      witness = make_wolfe_witness(
        fn = function(x) 1 - x,
        gr = function(x) -1
      ),
      step0 = 1,
      g = -1
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    result <- run_hager_zhang_witness(
      case$witness,
      step0 = case$step0
    )

    expect_equal(case$witness$fn_par(), c(0, 1), info = name)
    expect_equal(case$witness$gr_par(), c(0, 1), info = name)
    expect_budget_counts(result, case$witness)
    expect_equal(result$par, 1, info = name)
    expect_equal(result$f, 0, info = name)
    expect_equal(result$g, case$g, info = name)
  }
})

test_that("Hager-Zhang rejects unusable initial exhausted trials", {
  cases <- list(
    equal = list(
      fn = function(x) (x - 1)^2,
      gr = function(x) 2 * (x - 1),
      trial = 2,
      g = -2
    ),
    nonfinite_objective = list(
      fn = function(x) if (x == 0) 1 else Inf,
      gr = function(x) -1,
      trial = 1,
      g = -1
    ),
    nonfinite_derivative = list(
      fn = function(x) 1 - x,
      gr = function(x) if (x == 0) -1 else Inf,
      trial = 1,
      g = -1
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    witness <- make_wolfe_witness(case$fn, case$gr)
    result <- run_hager_zhang_witness(witness)

    expect_equal(witness$fn_par(), c(0, case$trial), info = name)
    expect_equal(witness$gr_par(), c(0, case$trial), info = name)
    expect_identical(
      witness$counts(),
      c(fn = 2L, gr = 2L),
      info = name
    )
    expect_budget_counts(result, witness)
    expect_equal(result$par, 0, info = name)
    expect_equal(result$f, 1, info = name)
    expect_equal(result$g, case$g, info = name)
  }
})

test_that("Hager-Zhang preserves a later exhausted decrease", {
  witness <- make_wolfe_witness(
    fn = function(x) 1 - x,
    gr = function(x) -1
  )
  result <- run_hager_zhang_witness(
    witness,
    ls_max_fn = 2,
    ls_max_gr = 2,
    ls_max_fg = 4
  )

  expect_equal(witness$fn_par(), c(0, 1, 5))
  expect_equal(witness$gr_par(), c(0, 1, 5))
  expect_identical(witness$counts(), c(fn = 3L, gr = 3L))
  expect_budget_counts(result, witness)
  expect_equal(result$par, 5)
  expect_equal(result$f, -4)
  expect_equal(result$g, -1)
})

test_that("Schmidt Armijo rejects a failure at the exact global function cap", {
  quartic <- function(x) x^4
  witness <- make_schmidt_armijo_witness()

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "backtracking",
    step0 = 1,
    ls_max_fn = 1,
    max_fn = 2,
    max_iter = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1, -3))
  expect_equal(witness$gr_par(), c(1, -3))
  expect_budget_counts(result, witness)
  expect_equal(result$terminate$what, "max_fn")
  expect_equal(result$par, 1)
  expect_equal(quartic(result$par), quartic(1))
  expect_equal(result$f, quartic(result$par))
})

test_that("Schmidt Armijo rejects equal and nonfinite terminal trials", {
  quartic <- function(x) x^4
  cases <- list(
    equal = list(
      fn = quartic,
      ls_max_fn = 2,
      fn_par = c(1, -3, -1)
    ),
    nonfinite = list(
      fn = function(x) if (x == -3) Inf else quartic(x),
      ls_max_fn = 1,
      fn_par = c(1, -3)
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    witness <- make_schmidt_armijo_witness(fn = case$fn)

    result <- suppressMessages(mize(
      1,
      witness$fg,
      method = "SD",
      line_search = "backtracking",
      step0 = 1,
      step_down = 0.5,
      ls_max_fn = case$ls_max_fn,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    ))

    expect_equal(witness$fn_par(), case$fn_par, info = name)
    expect_equal(witness$gr_par(), 1, info = name)
    expect_budget_counts(result, witness)
    expect_equal(result$par, 1, info = name)
    expect_equal(quartic(result$par), quartic(1), info = name)
    expect_equal(result$f, quartic(result$par), info = name)
  }
})

test_that("Schmidt Armijo returns an earlier decreasing fallback", {
  objective <- function(x) {
    0.628125 * x^3 + 1.621875 * x^2 - 1.128125 * x - 0.121875
  }
  gradient <- function(x) {
    3 * 0.628125 * x^2 + 2 * 1.621875 * x - 1.128125
  }
  witness <- make_schmidt_armijo_witness(objective, gradient)

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "backtracking",
    step0 = 1,
    step_down = 0.5,
    c1 = 0.1,
    ls_max_fn = 2,
    max_iter = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1, -3, -1))
  expect_equal(witness$gr_par(), 1)
  expect_budget_counts(result, witness)
  expect_equal(result$par, -3)
  expect_lt(objective(result$par), objective(1))
  expect_equal(result$f, objective(result$par))
})

test_that("Schmidt Armijo accepts an Armijo step on the final callback", {
  quartic <- function(x) x^4
  witness <- make_schmidt_armijo_witness()

  result <- mize(
    1,
    witness$fg,
    method = "SD",
    line_search = "backtracking",
    step0 = 1,
    step_down = 0.5,
    c1 = 0.1,
    ls_max_fn = 3,
    max_fn = 4,
    max_iter = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(witness$fn_par(), c(1, -3, -1, 0))
  expect_equal(witness$gr_par(), 1)
  expect_budget_counts(result, witness)
  expect_equal(result$terminate$what, "max_fn")
  expect_equal(result$par, 0)
  expect_lt(quartic(result$par), quartic(1))
  expect_equal(result$f, quartic(result$par))
})

test_that("Schmidt Armijo caches only its returned parameter value", {
  quartic <- function(x) x^4
  rejected_witness <- make_schmidt_armijo_witness()
  rejected_opt <- make_mize(
    method = "SD",
    line_search = "backtracking",
    step0 = 1,
    ls_max_fn = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  rejected_opt <- mize_init(rejected_opt, 1, rejected_witness$fg)

  rejected <- mize_step(rejected_opt, 1, rejected_witness$fg)

  expect_equal(rejected$par, 1)
  expect_equal(rejected$opt$cache$fn_new, quartic(rejected$par))
  expect_false(identical(rejected$opt$cache$fn_new, quartic(-3)))
  expect_equal(rejected$opt$cache$gr_curr, 4)
  expect_equal(rejected$opt$cache$gr_curr_iter, 2)
  expect_budget_counts(rejected, rejected_witness)

  objective <- function(x) {
    0.628125 * x^3 + 1.621875 * x^2 - 1.128125 * x - 0.121875
  }
  gradient <- function(x) {
    3 * 0.628125 * x^2 + 2 * 1.621875 * x - 1.128125
  }
  fallback_witness <- make_schmidt_armijo_witness(objective, gradient)
  fallback_opt <- make_mize(
    method = "SD",
    line_search = "backtracking",
    step0 = 1,
    step_down = 0.5,
    c1 = 0.1,
    ls_max_fn = 2,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  fallback_opt <- mize_init(fallback_opt, 1, fallback_witness$fg)

  fallback <- mize_step(fallback_opt, 1, fallback_witness$fg)

  expect_equal(fallback$par, -3)
  expect_equal(fallback$opt$cache$fn_new, objective(fallback$par))
  expect_equal(fallback$opt$cache$fn_curr, objective(fallback$par))
  expect_false(identical(fallback$opt$cache$gr_curr_iter, 2))
  expect_budget_counts(fallback, fallback_witness)
})

test_that("budget rollback preserves a stationary gradient across searches", {
  for (search in c("bold", "More-Thuente")) {
    witness <- make_budget_witness()

    result <- mize(
      0,
      witness$fg,
      method = "SD",
      line_search = search,
      max_iter = 3,
      max_fg = 1,
      check_conv_every = 1,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = 0,
      ginf_tol = NULL,
      step_tol = NULL
    )

    expect_identical(
      witness$counts(),
      c(fn = 0, gr = 1),
      info = search
    )
    expect_budget_counts(result, witness)
    expect_equal(result$par, 0, info = search)
    expect_equal(result$terminate$what, "grad_tol", info = search)
    expect_true(result$converged, info = search)
  }
})

test_that("the last permitted callback keeps convergence and failure precedence", {
  gradient_witness <- make_budget_witness()
  gradient_result <- mize(
    1,
    gradient_witness$fg,
    method = "SD",
    line_search = "constant",
    step0 = 0.5,
    max_iter = 3,
    max_gr = 2,
    check_conv_every = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = 0,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_identical(gradient_witness$counts()[["gr"]], 2)
  expect_budget_counts(gradient_result, gradient_witness)
  expect_equal(gradient_result$terminate$what, "grad_tol")
  expect_true(gradient_result$converged)

  function_witness <- make_budget_witness(fn = function(x) Inf)
  function_result <- mize(
    1,
    function_witness$fg,
    method = "SD",
    max_iter = 0,
    max_fn = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    step_tol = NULL
  )

  expect_identical(function_witness$counts(), c(fn = 1, gr = 0))
  expect_budget_counts(function_result, function_witness)
  expect_equal(function_result$terminate$what, "fn_inf")
  expect_equal(function_result$status, "failed")
})

test_that("an initial zero gradient keeps convergence precedence at its cap", {
  # The zero gradient makes the accepted constant-step update exactly zero, so
  # the initial gradient remains valid at the unchanged parameter.
  for (budget in c("max_gr", "max_fg")) {
    witness <- make_budget_witness()
    expected <- if (budget == "max_gr") {
      c(fn = 1, gr = 1)
    } else {
      c(fn = 0, gr = 1)
    }
    args <- list(
      par = 0,
      fg = witness$fg,
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      max_iter = 3,
      check_conv_every = 1,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = 0,
      ginf_tol = NULL,
      step_tol = NULL
    )
    args[[budget]] <- 1

    result <- do.call(mize, args)

    expect_identical(
      witness$counts(),
      expected,
      info = budget
    )
    expect_budget_counts(result, witness)
    expect_equal(result$terminate$what, "grad_tol", info = budget)
    expect_true(result$converged, info = budget)
    expect_equal(result$par, 0, info = budget)
  }
})

test_that("a changed parameter does not inherit the starting gradient", {
  witness <- make_budget_witness()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, 1, witness$fg)

  result <- mize_step(opt, 1, witness$fg)

  expect_identical(witness$counts(), c(fn = 0, gr = 1))
  expect_equal(result$par, 0.8)
  expect_true(has_gr_curr(result$opt, 1))
  expect_false(has_gr_curr(result$opt, 2))
})

test_that("initial and periodic summaries share the global function budget", {
  witness <- make_budget_witness()

  expect_message(
    result <- mize(
      1,
      witness$fg,
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      max_iter = 3,
      max_fn = 3,
      check_conv_every = 1,
      log_every = 1,
      abs_tol = 0.3,
      rel_tol = NULL,
      grad_tol = NULL,
      step_tol = NULL,
      verbose = TRUE,
      store_progress = TRUE
    ),
    "iter"
  )

  expect_identical(witness$counts(), c(fn = 3, gr = 2))
  expect_budget_counts(result, witness)
  expect_equal(result$terminate$what, "abs_tol")
  expect_true(result$converged)
  expect_equal(result$progress$nf, c(1, 2, 3))
})

test_that("progress and final reporting do not spend an exhausted gradient budget", {
  witness <- make_budget_witness()

  expect_message(
    result <- mize(
      1,
      witness$fg,
      method = "SD",
      line_search = "constant",
      step0 = 0.1,
      max_iter = 2,
      max_gr = 1,
      check_conv_every = 1,
      log_every = 2,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL,
      verbose = TRUE,
      store_progress = TRUE
    ),
    "iter 1"
  )

  expect_identical(witness$counts(), c(fn = 0, gr = 1))
  expect_budget_counts(result, witness)
  expect_equal(result$terminate$what, "max_gr")
  expect_equal(max(result$progress$ng), 1)
})

test_that("stateful summaries honor zero and exact callback budgets", {
  cases <- list(
    list(budget = "max_fn", calc_fn = TRUE, calc_gr = FALSE),
    list(budget = "max_gr", calc_fn = FALSE, calc_gr = TRUE),
    list(budget = "max_fg", calc_fn = FALSE, calc_gr = TRUE)
  )

  for (case in cases) {
    for (limit in 0:1) {
      witness <- make_budget_witness()
      opt <- make_mize(
        method = "SD",
        line_search = "constant",
        step0 = 0.1,
        abs_tol = NULL,
        rel_tol = NULL,
        grad_tol = NULL,
        ginf_tol = NULL,
        step_tol = NULL
      )
      init_args <- list(opt = opt, par = 1, fg = witness$fg)
      init_args[[case$budget]] <- limit
      opt <- do.call(mize_init, init_args)

      if (limit == 0) {
        step <- mize_step(opt, 1, witness$fg)
        expect_equal(step$par, 1)
        expect_identical(witness$counts(), c(fn = 0, gr = 0))
        opt <- step$opt
      }

      summary <- mize_step_summary(
        opt,
        1,
        witness$fg,
        calc_fn = case$calc_fn,
        calc_gr = case$calc_gr
      )
      checked <- check_mize_convergence(summary)

      expected <- if (limit == 0) 0 else 1
      actual <- if (case$calc_fn) "fn" else "gr"
      expect_identical(witness$counts()[[actual]], expected)
      expect_equal(summary$nf, witness$counts()[["fn"]])
      expect_equal(summary$ng, witness$counts()[["gr"]])
      expect_equal(checked$terminate$what, case$budget)
    }
  }
})

test_that("cached summary values consume no remaining budget", {
  witness <- make_budget_witness()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    max_fn = 0,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, 1, witness$fg)
  opt <- set_fn_curr(opt, 1, 1)

  summary <- mize_step_summary(opt, 1, witness$fg, calc_fn = TRUE)

  expect_identical(witness$counts(), c(fn = 0, gr = 0))
  expect_equal(summary$f, 1)
  expect_equal(summary$nf, 0)
})

test_that("Hager-Zhang initializer respects local and global search budgets", {
  run_two_steps <- function(ls_max_fn, ls_max_gr, ls_max_fg, ...) {
    witness <- make_wolfe_witness(
      fn = function(x) 0.5 * x^2,
      gr = function(x) x
    )
    opt <- make_mize(
      method = "SD",
      line_search = "hager-zhang",
      step0 = 0.1,
      step_next_init = "hz",
      ls_max_fn = ls_max_fn,
      ls_max_gr = ls_max_gr,
      ls_max_fg = ls_max_fg,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL,
      ...
    )
    opt <- mize_init(opt, 1, witness$fg)
    first <- mize_step(opt, 1, witness$fg)
    calls_before <- witness$counts()
    second <- mize_step(first$opt, first$par, witness$fg)

    list(
      witness = witness,
      first = first,
      second = second,
      second_counts = witness$counts() - calls_before
    )
  }

  local <- run_two_steps(
    ls_max_fn = 1,
    ls_max_gr = 1,
    ls_max_fg = 2
  )
  expect_equal(local$first$par, 0.9)
  expect_equal(local$second$par, 0.72)
  expect_identical(local$second_counts, c(fn = 1L, gr = 1L))
  expect_equal(local$witness$fn_par(), c(1, 0.9, 0.72))
  expect_equal(local$witness$gr_par(), c(1, 0.9, 0.72))
  expect_budget_counts(local$second, local$witness)

  global <- run_two_steps(
    ls_max_fn = 2,
    ls_max_gr = 1,
    ls_max_fg = 3,
    max_fn = 3,
    max_gr = 3,
    max_fg = 6
  )
  expect_equal(global$second$par, 0.72)
  expect_identical(global$second_counts, c(fn = 1L, gr = 1L))
  expect_equal(global$witness$fn_par(), c(1, 0.9, 0.72))
  expect_equal(global$witness$gr_par(), c(1, 0.9, 0.72))
  expect_budget_counts(global$second, global$witness)
})

test_that("Hager-Zhang initializer recovers from a nonfinite objective probe", {
  witness <- make_wolfe_witness(
    fn = function(x) {
      if (isTRUE(all.equal(x, 0.891))) NA_real_ else 0.5 * x^2
    },
    gr = function(x) x
  )
  opt <- make_mize(
    method = "SD",
    line_search = "hager-zhang",
    step0 = 0.1,
    step_next_init = "hz",
    ls_max_fn = 2,
    ls_max_gr = 1,
    ls_max_fg = 3,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
  opt <- mize_init(opt, 1, witness$fg)
  first <- mize_step(opt, 1, witness$fg)
  calls_before <- witness$counts()

  second <- mize_step(first$opt, first$par, witness$fg)

  expect_equal(first$par, 0.9)
  expect_equal(second$par, 0.72)
  expect_identical(
    witness$counts() - calls_before,
    c(fn = 2L, gr = 1L)
  )
  expect_equal(witness$fn_par(), c(1, 0.9, 0.891, 0.72))
  expect_equal(witness$gr_par(), c(1, 0.9, 0.72))
  expect_budget_counts(second, witness)
})
