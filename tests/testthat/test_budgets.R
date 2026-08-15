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

test_that("budget exhaustion before an internal backtracking trial is safe", {
  witness <- make_budget_witness()
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = sd_direction(),
        step_size = backtracking()
      )
    )
  )
  opt <- mize_init(
    opt,
    1,
    witness$fg,
    max_fg = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  result <- mize_step(opt, 1, witness$fg)

  expect_identical(witness$counts(), c(fn = 0, gr = 1))
  expect_equal(c(fn = result$nf, gr = result$ng), witness$counts())
  expect_equal(result$opt$terminate$what, "max_fg")
  expect_equal(result$par, 1)
  expect_false("f" %in% names(result))
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
