make_curvature_count_witness <- function(callback = NULL, value = NULL) {
  calls <- new.env(parent = emptyenv())
  calls$hs <- 0L
  calls$hi <- 0L

  fg <- list(
    fn = function(x) sum(x^2) / 2,
    gr = function(x) x
  )
  if (!is.null(callback)) {
    fg[[callback]] <- function(x) {
      calls[[callback]] <- calls[[callback]] + 1L
      if (is.function(value)) value(x) else value
    }
  }

  list(
    fg = fg,
    counts = function() c(hs = calls$hs, hi = calls$hi)
  )
}

curvature_count_args <- function(fg, method, max_iter = 3, ...) {
  list(
    par = c(2, -1),
    fg = fg,
    method = method,
    line_search = "constant",
    step0 = 0.25,
    max_iter = max_iter,
    check_conv_every = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL,
    ...
  )
}

test_that("curvature counts match callback cadence and progress", {
  cases <- list(
    newton_hs = list(
      method = "NEWTON",
      callback = "hs",
      value = diag(2),
      initial = c(nh = 0, nhi = 0),
      final = c(nh = 3, nhi = 0)
    ),
    newton_hi = list(
      method = "NEWTON",
      callback = "hi",
      value = diag(2),
      initial = c(nh = 0, nhi = 0),
      final = c(nh = 0, nhi = 3)
    ),
    phess = list(
      method = "PHESS",
      callback = "hs",
      value = diag(2),
      initial = c(nh = 1, nhi = 0),
      final = c(nh = 1, nhi = 0)
    ),
    bfgs = list(
      method = "BFGS",
      callback = "hi",
      value = diag(2),
      initial = c(nh = 0, nhi = 1),
      final = c(nh = 0, nhi = 1)
    ),
    sr1 = list(
      method = "SR1",
      callback = "hi",
      value = diag(2),
      initial = c(nh = 0, nhi = 1),
      final = c(nh = 0, nhi = 1)
    ),
    lbfgs_vector = list(
      method = "L-BFGS",
      callback = "hi",
      value = c(1, 1),
      initial = c(nh = 0, nhi = 0),
      final = c(nh = 0, nhi = 3)
    ),
    lbfgs_matrix = list(
      method = "L-BFGS",
      callback = "hi",
      value = diag(2),
      initial = c(nh = 0, nhi = 0),
      final = c(nh = 0, nhi = 3)
    ),
    steepest_descent = list(
      method = "SD",
      callback = NULL,
      value = NULL,
      initial = c(nh = 0, nhi = 0),
      final = c(nh = 0, nhi = 0)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    witness <- make_curvature_count_witness(case$callback, case$value)
    result <- do.call(
      mize,
      curvature_count_args(
        witness$fg,
        case$method,
        store_progress = TRUE
      )
    )

    expect_equal(
      c(nh = result$nh, nhi = result$nhi),
      case$final,
      info = case_name
    )
    expect_equal(
      witness$counts(),
      c(hs = case$final[["nh"]], hi = case$final[["nhi"]]),
      info = case_name
    )
    expect_equal(
      unlist(result$progress[1, c("nh", "nhi")]),
      case$initial,
      info = case_name
    )
    expect_equal(
      unlist(result$progress[nrow(result$progress), c("nh", "nhi")]),
      case$final,
      info = case_name
    )
  }
})

test_that("one-shot and stateful curvature counts agree", {
  cases <- list(
    newton = list(method = "NEWTON", callback = "hi", value = diag(2)),
    phess = list(method = "PHESS", callback = "hs", value = diag(2))
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    one_shot_witness <- make_curvature_count_witness(case$callback, case$value)
    one_shot <- do.call(
      mize,
      curvature_count_args(one_shot_witness$fg, case$method)
    )

    stateful_witness <- make_curvature_count_witness(case$callback, case$value)
    opt <- make_mize(
      method = case$method,
      line_search = "constant",
      step0 = 0.25
    )
    opt <- mize_init(opt, c(2, -1), stateful_witness$fg)
    par <- c(2, -1)
    for (iter in 1:3) {
      step <- mize_step(opt, par, stateful_witness$fg)
      opt <- step$opt
      par <- step$par
    }
    summary <- mize_step_summary(opt, par, stateful_witness$fg)

    expect_equal(
      c(nh = step$nh, nhi = step$nhi),
      c(nh = one_shot$nh, nhi = one_shot$nhi),
      info = case_name
    )
    expect_equal(
      c(nh = summary$nh, nhi = summary$nhi),
      c(nh = one_shot$nh, nhi = one_shot$nhi),
      info = case_name
    )
    expect_equal(par, one_shot$par, info = case_name)
  }
})

test_that("reinitialization preserves and extends curvature counts", {
  witness <- make_curvature_count_witness("hi", diag(2))
  opt <- make_mize(method = "BFGS", line_search = "constant", step0 = 0.25)

  opt <- mize_init(opt, c(2, -1), witness$fg)
  expect_equal(opt$counts[c("hs", "hi")], list(hs = 0, hi = 1))

  step <- mize_step(opt, c(2, -1), witness$fg)
  opt <- mize_init(step$opt, c(2, -1), witness$fg)
  summary <- mize_step_summary(opt, c(2, -1), witness$fg)

  expect_equal(opt$counts[c("hs", "hi")], list(hs = 0, hi = 2))
  expect_equal(c(nh = summary$nh, nhi = summary$nhi), c(nh = 0, nhi = 2))
  expect_equal(witness$counts(), c(hs = 0, hi = 2))
})

test_that("an already-terminal step does not call curvature callbacks", {
  witness <- make_curvature_count_witness("hi", diag(2))
  opt <- make_mize(
    method = "NEWTON",
    line_search = "constant",
    step0 = 0.25,
    max_gr = 0,
    par = c(2, -1),
    fg = witness$fg
  )

  result <- mize_step(opt, c(2, -1), witness$fg)

  expect_identical(result$opt$terminate$what, "max_gr")
  expect_equal(c(nh = result$nh, nhi = result$nhi), c(nh = 0, nhi = 0))
  expect_equal(witness$counts(), c(hs = 0, hi = 0))
})

test_that("curvature callback failures preserve errors and attempted counts", {
  cases <- list(
    throwing = list(
      value = function(x) stop("curvature exploded"),
      pattern = "curvature exploded"
    ),
    malformed = list(
      value = "not an inverse Hessian",
      pattern = "fg\\$hi\\(par\\) must return"
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    witness <- make_curvature_count_witness("hi", case$value)

    expect_error(
      do.call(mize, curvature_count_args(witness$fg, "NEWTON", max_iter = 1)),
      case$pattern,
      info = case_name
    )
    expect_equal(witness$counts(), c(hs = 0, hi = 1), info = case_name)
  }
})

test_that("private periodic PHESS refreshes use native accounting", {
  # The refresh cadence is private. This protects the public invariant that
  # every accepted curvature callback in production source is counted once.
  witness <- make_curvature_count_witness("hs", diag(2))
  direction <- partial_hessian_direction(hessian_every = 1)
  opt <- list(
    cache = list(gr_curr = c(1, -1), gr_curr_iter = 1),
    counts = make_counts()
  )

  result <- direction$calculate(
    opt,
    stage = list(),
    sub_stage = direction,
    par = c(1, -1),
    fg = witness$fg,
    iter = 1
  )

  expect_equal(result$opt$counts[c("hs", "hi")], list(hs = 1, hi = 0))
  expect_equal(witness$counts(), c(hs = 1, hi = 0))
})
