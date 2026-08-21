optimizer_diagnostic_quadratic <- function() {
  hessian <- diag(c(1, 4))
  list(
    fn = function(x) {
      drop(0.5 * crossprod(x, hessian %*% x))
    },
    gr = function(x) {
      as.vector(hessian %*% x)
    },
    fg = function(x) {
      list(
        fn = drop(0.5 * crossprod(x, hessian %*% x)),
        gr = as.vector(hessian %*% x)
      )
    }
  )
}

test_that("Wolfe diagnostics have common public semantics", {
  fg <- optimizer_diagnostic_quadratic()
  line_searches <- c(
    "More-Thuente",
    "Rasmussen",
    "Schmidt",
    "Hager-Zhang"
  )
  numeric_diagnostics <- c("alpha_init", "slope_init", "ls_nf", "ls_ng")
  categorical_diagnostics <- c("ls_reason", "ls_outcome")

  for (line_search in line_searches) {
    result <- mize(
      par = c(2, -1),
      fg = fg,
      method = "CG",
      cg_update = "PR+",
      line_search = line_search,
      max_iter = 10,
      store_progress = TRUE,
      grad_tol = 1e-8
    )
    progress <- result$progress
    search_rows <- as.integer(rownames(progress)) > 0

    expect_true(
      all(c(numeric_diagnostics, categorical_diagnostics) %in% names(progress)),
      info = line_search
    )
    expect_true(
      all(vapply(progress[numeric_diagnostics], is.numeric, logical(1))),
      info = line_search
    )
    expect_true(
      all(vapply(progress[categorical_diagnostics], is.character, logical(1))),
      info = line_search
    )
    expect_true(
      all(is.na(progress[1, c(numeric_diagnostics, categorical_diagnostics)])),
      info = line_search
    )
    expect_true(all(progress$alpha_init[search_rows] > 0), info = line_search)
    expect_true(all(progress$slope_init[search_rows] < 0), info = line_search)
    expect_true(all(progress$alpha[search_rows] > 0), info = line_search)
    expect_identical(
      progress$ls_reason[search_rows],
      rep("wolfe", sum(search_rows)),
      info = line_search
    )
    expect_identical(
      progress$ls_outcome[search_rows],
      rep("wolfe", sum(search_rows)),
      info = line_search
    )
    expect_equal(result$nf, 1 + sum(progress$ls_nf, na.rm = TRUE))
    expect_equal(result$ng, 1 + sum(progress$ls_ng, na.rm = TRUE))
    expect_identical(result$terminate$what, "grad_tol")
    expect_true(result$converged)
  }
})

test_that("a nonstationary no-step search is a public failure", {
  fg <- optimizer_diagnostic_quadratic()

  result <- expect_no_warning(mize(
    par = c(2, -1),
    fg = fg,
    method = "CG",
    line_search = "More-Thuente",
    step0 = 1,
    ls_max_fn = 0,
    max_iter = 1,
    store_progress = TRUE,
    grad_tol = 1e-8
  ))
  final_progress <- result$progress[nrow(result$progress), , drop = FALSE]

  expect_identical(result$terminate$what, "line_search_failed")
  expect_identical(result$terminate$val, "budget_exhausted")
  expect_identical(result$status, "failed")
  expect_false(result$converged)
  expect_match(result$message, "budget_exhausted", fixed = TRUE)
  expect_identical(final_progress$ls_reason, "budget_exhausted")
  expect_identical(final_progress$ls_outcome, "no_step")
  expect_equal(final_progress$ls_nf, 0)
  expect_equal(final_progress$ls_ng, 0)
  expect_lt(final_progress$slope_init, 0)
  expect_equal(result$nf, 1)
  expect_equal(result$ng, 1)

  result_without_progress <- expect_no_warning(mize(
    par = c(2, -1),
    fg = fg,
    method = "CG",
    line_search = "More-Thuente",
    step0 = 1,
    ls_max_fn = 0,
    max_iter = 1,
    store_progress = FALSE,
    grad_tol = 1e-8
  ))
  expect_identical(result_without_progress$terminate, result$terminate)
  expect_identical(result_without_progress$status, result$status)
  expect_identical(result_without_progress$message, result$message)
  expect_equal(result_without_progress$nf, result$nf)
  expect_equal(result_without_progress$ng, result$ng)
  expect_null(result_without_progress$progress)
})

test_that("an improving non-Wolfe fallback remains usable", {
  fg <- optimizer_diagnostic_quadratic()
  initial_f <- fg$fn(c(2, -1))

  result <- expect_no_warning(mize(
    par = c(2, -1),
    fg = fg,
    method = "CG",
    line_search = "More-Thuente",
    step0 = 0.01,
    ls_max_fn = 1,
    max_iter = 1,
    store_progress = TRUE,
    grad_tol = 1e-8
  ))
  final_progress <- result$progress[nrow(result$progress), , drop = FALSE]

  expect_lt(result$f, initial_f)
  expect_identical(result$terminate$what, "max_iter")
  expect_identical(result$status, "budget_exhausted")
  expect_identical(final_progress$ls_reason, "budget_exhausted")
  expect_identical(final_progress$ls_outcome, "improving_fallback")
  expect_equal(final_progress$ls_nf, 1)
  expect_equal(final_progress$ls_ng, 1)
  expect_equal(result$nf, 2)
  expect_equal(result$ng, 2)
})

test_that("a stationary start does not fabricate line-search diagnostics", {
  fg <- optimizer_diagnostic_quadratic()

  result <- expect_no_warning(mize(
    par = c(0, 0),
    fg = fg,
    method = "CG",
    line_search = "Hager-Zhang",
    max_iter = 5,
    store_progress = TRUE,
    grad_tol = 1e-8
  ))

  expect_identical(result$terminate$what, "grad_tol")
  expect_true(result$converged)
  expect_equal(result$iter, 0)
  expect_false(any(
    c(
      "alpha_init",
      "slope_init",
      "ls_reason",
      "ls_outcome",
      "ls_nf",
      "ls_ng"
    ) %in%
      names(result$progress)
  ))
})

test_that("the exported step workflow classifies a realized no-step failure", {
  fg <- optimizer_diagnostic_quadratic()
  par0 <- c(2, -1)
  opt <- make_mize(
    method = "CG",
    line_search = "More-Thuente",
    step0 = 1,
    ls_max_fn = 0,
    par = par0,
    fg = fg,
    max_iter = Inf,
    grad_tol = 1e-8,
    step_tol = sqrt(.Machine$double.eps)
  )

  step_result <- mize_step(opt, par0, fg)
  step_info <- mize_step_summary(
    step_result$opt,
    step_result$par,
    fg,
    par0
  )
  checked <- check_mize_convergence(step_info)

  expect_equal(step_info$step, 0)
  expect_identical(step_info$ls_outcome, "no_step")
  expect_identical(checked$terminate$what, "line_search_failed")
  expect_identical(checked$terminate$val, "budget_exhausted")
  expect_identical(checked$status, "failed")
  expect_false(checked$converged)
})

test_that("summary callbacks retain global budget precedence over no-step", {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0L
  calls$gr <- 0L
  fg <- list(
    fn = function(x) {
      calls$fn <- calls$fn + 1L
      1
    },
    gr = function(x) {
      calls$gr <- calls$gr + 1L
      rep(1, length(x))
    }
  )

  result <- mize(
    par = 1,
    fg = fg,
    method = "MOM",
    line_search = "Rasmussen",
    step0 = 1,
    mom_schedule = 0.9,
    ls_max_fn = 7,
    max_iter = 2,
    max_fn = 9
  )

  expect_identical(result$terminate, list(what = "max_fn", val = 9))
  expect_identical(result$status, "budget_exhausted")
  expect_false(result$converged)
  expect_equal(result$nf, calls$fn)
  expect_equal(result$ng, calls$gr)
  expect_equal(result$nf, 9)
})

test_that("no-step budget wins over function-difference convergence", {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0L
  calls$gr <- 0L
  fg <- list(
    fn = function(x) {
      calls$fn <- calls$fn + 1L
      0.5 * x^2
    },
    gr = function(x) {
      calls$gr <- calls$gr + 1L
      x
    }
  )

  result <- mize(
    par = 1,
    fg = fg,
    method = "MOM",
    line_search = "More-Thuente",
    step0 = 0.02,
    mom_schedule = 0,
    ls_max_fn = 1,
    max_iter = 20,
    max_fn = 13,
    store_progress = TRUE
  )
  final_progress <- result$progress[nrow(result$progress), , drop = FALSE]
  previous_progress <- result$progress[
    nrow(result$progress) - 1,
    ,
    drop = FALSE
  ]

  expect_gte(result$iter, 2)
  expect_identical(final_progress$ls_outcome, "no_step")
  expect_equal(final_progress$step, 0)
  expect_equal(final_progress$f, previous_progress$f)
  expect_identical(result$terminate, list(what = "max_fn", val = 13))
  expect_identical(result$status, "budget_exhausted")
  expect_false(result$converged)
  expect_equal(result$nf, calls$fn)
  expect_equal(result$ng, calls$gr)
  expect_equal(result$nf, 13)
})

test_that("manual no-step adjudication materializes an exact budget", {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0L
  calls$gr <- 0L
  fg <- list(
    fn = function(x) {
      calls$fn <- calls$fn + 1L
      1
    },
    gr = function(x) {
      calls$gr <- calls$gr + 1L
      rep(1, length(x))
    }
  )
  par0 <- 1
  opt <- make_mize(
    method = "CG",
    line_search = "More-Thuente",
    step0 = 1,
    ls_max_fn = 1,
    par = par0,
    fg = fg,
    max_iter = Inf,
    max_fn = 2,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = sqrt(.Machine$double.eps)
  )

  step_result <- mize_step(opt, par0, fg)
  step_info <- mize_step_summary(
    step_result$opt,
    step_result$par,
    fg,
    par0
  )
  checked <- check_mize_convergence(step_info)

  expect_identical(step_info$ls_outcome, "no_step")
  expect_equal(step_info$step, 0)
  expect_identical(checked$terminate, list(what = "max_fn", val = 2))
  expect_identical(checked$status, "budget_exhausted")
  expect_false(checked$converged)
  expect_equal(step_info$nf, calls$fn)
  expect_equal(step_info$ng, calls$gr)
  expect_equal(step_info$nf, 2)
})

test_that("a no-step substage does not cancel a realized momentum transition", {
  hessian <- diag(c(1, 2))
  fg <- list(
    fn = function(x) {
      drop(0.5 * crossprod(x, hessian %*% x))
    },
    gr = function(x) {
      as.vector(hessian %*% x)
    }
  )

  result <- mize(
    par = rep(sqrt(8), 2),
    fg = fg,
    method = "MOM",
    line_search = "Rasmussen",
    step0 = 1 / sqrt(10),
    mom_schedule = 0.9,
    ls_max_fn = 1,
    max_iter = 2,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = sqrt(.Machine$double.eps),
    store_progress = TRUE
  )
  final_progress <- result$progress[nrow(result$progress), , drop = FALSE]

  expect_identical(final_progress$ls_outcome, "no_step")
  expect_equal(final_progress$alpha, 0)
  expect_equal(final_progress$step, 1.8)
  expect_identical(result$terminate$what, "max_iter")
  expect_identical(result$status, "budget_exhausted")
  expect_false(result$converged)
})

test_that("Schmidt Armijo reports selected-point outcomes", {
  fg <- optimizer_diagnostic_quadratic()
  cases <- list(
    list(
      step0 = 0.1,
      ls_max_fn = 1,
      outcome = "armijo",
      reason = "armijo",
      termination = "max_iter"
    ),
    list(
      step0 = 0.58823,
      ls_max_fn = 1,
      outcome = "improving_fallback",
      reason = "budget_exhausted",
      termination = "max_iter"
    ),
    list(
      step0 = 0.1,
      ls_max_fn = 0,
      outcome = "no_step",
      reason = "budget_exhausted",
      termination = "line_search_failed"
    )
  )

  for (case in cases) {
    result <- mize(
      par = c(2, -1),
      fg = fg,
      method = "CG",
      line_search = "backtracking",
      step0 = case$step0,
      ls_max_fn = case$ls_max_fn,
      max_iter = 1,
      store_progress = TRUE,
      grad_tol = 1e-8
    )
    final_progress <- result$progress[nrow(result$progress), , drop = FALSE]

    expect_identical(final_progress$ls_outcome, case$outcome)
    expect_identical(final_progress$ls_reason, case$reason)
    expect_identical(result$terminate$what, case$termination)
  }
})

test_that("exact Newton reports public direction provenance", {
  objective_hessian <- diag(c(2, 4))
  inverse_hessian <- diag(1 / diag(objective_hessian))
  cases <- list(
    hessian_matrix = list(
      hs = objective_hessian,
      reason = "hessian_solve"
    ),
    hessian_diagonal = list(
      hs = diag(objective_hessian),
      reason = "hessian_diagonal"
    ),
    inverse_hessian_matrix = list(
      hi = inverse_hessian,
      reason = "inverse_hessian_multiply"
    ),
    inverse_hessian_diagonal = list(
      hi = diag(inverse_hessian),
      reason = "inverse_hessian_diagonal"
    ),
    cholesky_failure = list(
      hs = diag(c(1, -1)),
      reason = "cholesky_fallback"
    ),
    direction_check_failure = list(
      hi = -diag(2),
      reason = "direction_check_fallback"
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    fg <- list(
      fn = function(x) {
        drop(0.5 * crossprod(x, objective_hessian %*% x))
      },
      gr = function(x) {
        as.vector(objective_hessian %*% x)
      }
    )
    if (!is.null(case$hs)) {
      fg$hs <- function(x) case$hs
    }
    if (!is.null(case$hi)) {
      fg$hi <- function(x) case$hi
    }

    result <- mize(
      par = c(1, -1),
      fg = fg,
      method = "NEWTON",
      line_search = "More-Thuente",
      step0 = 1,
      max_iter = 1,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL,
      store_progress = TRUE
    )
    progress <- result$progress

    expect_true(is.character(progress$direction_reason), info = case_name)
    expect_true(is.na(progress$direction_reason[1]), info = case_name)
    expect_identical(
      progress$direction_reason[nrow(progress)],
      case$reason,
      info = case_name
    )
  }
})

test_that("direction provenance is absent for non-Newton methods", {
  result <- mize(
    par = c(2, -1),
    fg = optimizer_diagnostic_quadratic(),
    method = "BFGS",
    max_iter = 1,
    store_progress = TRUE
  )

  expect_false("direction_reason" %in% names(result$progress))
})
