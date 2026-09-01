validation_fg <- function(calls = NULL) {
  list(
    fn = function(x) {
      if (!is.null(calls)) {
        calls$fn <- calls$fn + 1
      }
      0
    },
    gr = function(x) {
      if (!is.null(calls)) {
        calls$gr <- calls$gr + 1
      }
      rep(0, length(x))
    }
  )
}

validation_result_witness <- function(
  fn = function(x) sum(x^2),
  gr = function(x) 2 * x,
  fg = NULL
) {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0L
  calls$gr <- 0L
  calls$fg <- 0L

  callbacks <- list(
    fn = function(x) {
      calls$fn <- calls$fn + 1L
      fn(x)
    },
    gr = function(x) {
      calls$gr <- calls$gr + 1L
      gr(x)
    }
  )
  if (!is.null(fg)) {
    callbacks$fg <- function(x) {
      calls$fg <- calls$fg + 1L
      fg(x)
    }
  }

  list(
    fg = callbacks,
    counts = function() {
      c(fn = calls$fn, gr = calls$gr, fg = calls$fg)
    }
  )
}

validation_wolfe_args <- function(fg) {
  list(
    par = c(1, -1),
    fg = fg,
    method = "SD",
    line_search = "more-thuente",
    step0 = 0.1,
    max_iter = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
}

curvature_validation_witness <- function(
  callback,
  value,
  fn = function(x) sum(x^2) / 2,
  gr = function(x) x
) {
  calls <- new.env(parent = emptyenv())
  calls$hs <- 0L
  calls$hi <- 0L

  fg <- list(fn = fn, gr = gr)
  fg[[callback]] <- function(x) {
    calls[[callback]] <- calls[[callback]] + 1L
    value
  }

  list(
    fg = fg,
    counts = function() {
      c(hs = calls$hs, hi = calls$hi)
    }
  )
}

curvature_validation_args <- function(
  fg,
  method,
  par = c(1, -1, 2, -2),
  step0 = 1
) {
  list(
    par = par,
    fg = fg,
    method = method,
    line_search = "constant",
    step0 = step0,
    max_iter = 1,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )
}

expect_callback_validation_error <- function(code, pattern, info = NULL) {
  result <- tryCatch(code(), error = identity)
  expect_true(inherits(result, "error"), info = info)
  if (inherits(result, "error")) {
    expect_match(conditionMessage(result), pattern, info = info)
  }
  invisible(result)
}

test_that("objective callback results are numeric dimensionless scalars", {
  malformed <- list(
    vector = c(1, 2),
    dimensional = matrix(1, nrow = 1)
  )

  for (case_name in names(malformed)) {
    witness <- validation_result_witness(
      fn = function(x) malformed[[case_name]]
    )

    expect_callback_validation_error(
      function() {
        mize(
          c(1, -1),
          witness$fg,
          method = "SD",
          line_search = "constant",
          step0 = 0.1,
          max_iter = 0,
          check_conv_every = NULL
        )
      },
      "fg\\$fn\\(par\\).*numeric scalar.*no dimensions",
      info = case_name
    )
    expect_identical(
      witness$counts(),
      c(fn = 1L, gr = 0L, fg = 0L),
      info = case_name
    )
  }
})

test_that("gradient callback results match the parameter vector shape", {
  par <- c(1, -1)

  short <- validation_result_witness(gr = function(x) 1)
  expect_callback_validation_error(
    function() {
      mize(
        par,
        short$fg,
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
      )
    },
    "fg\\$gr\\(par\\).*numeric vector.*length.*par"
  )
  expect_identical(short$counts(), c(fn = 0L, gr = 1L, fg = 0L))

  long <- validation_result_witness(gr = function(x) c(1, 2, 3))
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = par,
    fg = long$fg
  )
  expect_callback_validation_error(
    function() mize_step(opt, par, long$fg),
    "fg\\$gr\\(par\\).*numeric vector.*length.*par"
  )
  expect_identical(long$counts(), c(fn = 0L, gr = 1L, fg = 0L))

  dimensional <- validation_result_witness(
    gr = function(x) matrix(c(1, 2), nrow = 1)
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = par,
    fg = dimensional$fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL
  )
  expect_callback_validation_error(
    function() {
      mize_step_summary(opt, par, dimensional$fg, calc_gr = TRUE)
    },
    "fg\\$gr\\(par\\).*numeric vector.*no dimensions"
  )
  expect_identical(dimensional$counts(), c(fn = 0L, gr = 1L, fg = 0L))
})

test_that("combined callback results contain valid exact components", {
  malformed <- list(
    non_list = function(x) 1,
    missing_fn = function(x) list(gr = 2 * x),
    missing_gr = function(x) list(fn = sum(x^2)),
    malformed_fn = function(x) list(fn = c(1, 2), gr = 2 * x),
    malformed_gr = function(x) list(fn = sum(x^2), gr = 1)
  )
  patterns <- c(
    non_list = "fg\\$fg\\(par\\).*return a list",
    missing_fn = "fg\\$fg\\(par\\).*exact.*fn",
    missing_gr = "fg\\$fg\\(par\\).*exact.*gr",
    malformed_fn = "fg\\$fg\\(par\\)\\$fn.*numeric scalar",
    malformed_gr = "fg\\$fg\\(par\\)\\$gr.*numeric vector.*length.*par"
  )

  for (case_name in names(malformed)) {
    witness <- validation_result_witness(fg = malformed[[case_name]])

    expect_callback_validation_error(
      function() do.call(mize, validation_wolfe_args(witness$fg)),
      patterns[[case_name]],
      info = case_name
    )
    expect_identical(
      witness$counts(),
      c(fn = 1L, gr = 1L, fg = 1L),
      info = case_name
    )
  }
})

test_that("TN validates each finite-difference gradient result once", {
  calls <- 0L
  witness <- validation_result_witness(gr = function(x) {
    calls <<- calls + 1L
    if (calls == 1L) 2 * x else 1
  })

  expect_callback_validation_error(
    function() {
      mize(
        c(1, -1),
        witness$fg,
        method = "TN",
        line_search = "constant",
        step0 = 0.1,
        max_iter = 1,
        check_conv_every = NULL,
        abs_tol = NULL,
        rel_tol = NULL,
        grad_tol = NULL,
        ginf_tol = NULL,
        step_tol = NULL
      )
    },
    "fg\\$gr\\(par\\).*numeric vector.*length.*par"
  )
  expect_identical(witness$counts(), c(fn = 0L, gr = 2L, fg = 0L))
})

test_that("uncounted summaries validate their direct callback results", {
  # count_res_fg is an internal test mode that bypasses the calc helpers. This
  # protects the direct summary-consumption boundary used by that mode.
  par <- c(1, -1)

  objective <- validation_result_witness(
    fn = function(x) matrix(sum(x^2), nrow = 1)
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = par,
    fg = objective$fg
  )
  opt$count_res_fg <- FALSE
  expect_callback_validation_error(
    function() mize_step_summary(opt, par, objective$fg, calc_fn = TRUE),
    "fg\\$fn\\(par\\).*numeric scalar.*no dimensions"
  )
  expect_identical(objective$counts(), c(fn = 1L, gr = 0L, fg = 0L))

  gradient <- validation_result_witness(
    gr = function(x) matrix(2 * x, nrow = 1)
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = par,
    fg = gradient$fg
  )
  opt$count_res_fg <- FALSE
  expect_callback_validation_error(
    function() mize_step_summary(opt, par, gradient$fg, calc_gr = TRUE),
    "fg\\$gr\\(par\\).*numeric vector.*no dimensions"
  )
  expect_identical(gradient$counts(), c(fn = 0L, gr = 1L, fg = 0L))
})

test_that("valid separate callback results retain integer and double types", {
  callbacks <- list(
    integer = list(
      fn = function(x) as.integer(sum(x^2)),
      gr = function(x) as.integer(2 * x)
    ),
    double = list(
      fn = function(x) sum(x^2),
      gr = function(x) 2 * x
    )
  )

  for (case_name in names(callbacks)) {
    callback <- callbacks[[case_name]]
    witness <- validation_result_witness(fn = callback$fn, gr = callback$gr)
    result <- mize(
      c(1, -1),
      witness$fg,
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
    )

    expect_equal(result$par, c(0.8, -0.8), info = case_name)
    expect_identical(
      witness$counts(),
      c(fn = 1L, gr = 1L, fg = 0L),
      info = case_name
    )
  }
})

test_that("valid combined callbacks allow extra components", {
  witness <- validation_result_witness(fg = function(x) {
    list(fn = sum(x^2), gr = 2 * x, shared = "allowed")
  })

  result <- do.call(mize, validation_wolfe_args(witness$fg))

  expect_equal(result$par, c(0, 0))
  expect_equal(result$f, 0)
  expect_identical(witness$counts(), c(fn = 1L, gr = 1L, fg = 2L))
})

test_that("shaped nonfinite results retain optimizer handling", {
  objective <- validation_result_witness(fn = function(x) Inf)
  fn_result <- mize(
    c(1, -1),
    objective$fg,
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    max_iter = 0,
    check_conv_every = NULL
  )
  expect_equal(fn_result$terminate$what, "fn_inf")
  expect_equal(fn_result$status, "failed")
  expect_identical(objective$counts(), c(fn = 1L, gr = 0L, fg = 0L))

  gradient <- validation_result_witness(gr = function(x) c(NaN, Inf))
  gr_result <- mize(
    c(1, -1),
    gradient$fg,
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    max_iter = 1,
    check_conv_every = NULL
  )
  expect_equal(gr_result$terminate$what, "gr_inf")
  expect_equal(gr_result$status, "failed")
  expect_identical(gradient$counts(), c(fn = 1L, gr = 1L, fg = 0L))
})

test_that("shaped nonfinite Wolfe trials remain available for recovery", {
  trial <- 0L
  witness <- validation_result_witness(fg = function(x) {
    trial <<- trial + 1L
    if (trial == 1L) {
      return(list(fn = Inf, gr = rep(Inf, length(x))))
    }
    list(fn = sum(x^2), gr = 2 * x)
  })

  result <- do.call(mize, validation_wolfe_args(witness$fg))

  expect_equal(result$par, c(0, 0))
  expect_equal(result$f, 0)
  expect_identical(witness$counts(), c(fn = 1L, gr = 1L, fg = 3L))
})

test_that("NEWTON accepts existing Hessian vector and matrix forms", {
  par <- c(1, -1, 2, -2)
  weights <- c(1, 2, 1, 2)
  fn <- function(x) sum(weights * x^2) / 2
  gr <- function(x) weights * x
  hessians <- list(
    integer_vector = as.integer(weights),
    double_vector = as.double(weights),
    integer_matrix = diag(as.integer(weights)),
    double_matrix = diag(as.double(weights)),
    integer_repeated_block = diag(c(1L, 2L)),
    double_repeated_block = diag(c(1, 2))
  )

  for (case_name in names(hessians)) {
    witness <- curvature_validation_witness(
      "hs",
      hessians[[case_name]],
      fn,
      gr
    )
    result <- do.call(
      mize,
      curvature_validation_args(witness$fg, "NEWTON", par)
    )

    expect_equal(result$par, rep(0, length(par)), info = case_name)
    expect_identical(witness$counts(), c(hs = 1L, hi = 0L), info = case_name)
  }
})

test_that("NEWTON rejects malformed Hessian results at first consumption", {
  malformed <- list(
    wrong_type = "invalid",
    short_vector = c(1, 2),
    long_vector = rep(1, 8),
    nonsquare_matrix = matrix(1, nrow = 2, ncol = 3),
    nondividing_matrix = diag(3),
    nonfinite_vector = c(1, 1, Inf, 1),
    nonfinite_matrix = diag(c(1, 1, NA_real_, 1)),
    asymmetric_matrix = matrix(
      c(2, 0, 0, 0, 0.25, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2),
      nrow = 4
    )
  )
  patterns <- c(
    wrong_type = "fg\\$hs\\(par\\).*finite numeric vector.*length 4",
    short_vector = "fg\\$hs\\(par\\).*vector.*length 4",
    long_vector = "fg\\$hs\\(par\\).*vector.*length 4",
    nonsquare_matrix = "fg\\$hs\\(par\\).*symmetric square numeric matrix",
    nondividing_matrix = "fg\\$hs\\(par\\).*dimension.*divides 4",
    nonfinite_vector = "fg\\$hs\\(par\\).*finite",
    nonfinite_matrix = "fg\\$hs\\(par\\).*finite",
    asymmetric_matrix = "fg\\$hs\\(par\\).*symmetric"
  )

  for (case_name in names(malformed)) {
    witness <- curvature_validation_witness("hs", malformed[[case_name]])

    expect_callback_validation_error(
      function() {
        do.call(mize, curvature_validation_args(witness$fg, "NEWTON"))
      },
      patterns[[case_name]],
      info = case_name
    )
    expect_identical(
      witness$counts(),
      c(hs = 1L, hi = 0L),
      info = case_name
    )
  }
})

test_that("NEWTON Hessian symmetry allows floating-point noise", {
  noise <- 10 * .Machine$double.eps
  hessian <- matrix(c(2, 0.5, 0.5 + noise, 1), nrow = 2)
  witness <- curvature_validation_witness(
    "hs",
    hessian,
    fn = function(x) sum(x^2) / 2,
    gr = function(x) x
  )

  result <- do.call(
    mize,
    curvature_validation_args(
      witness$fg,
      "NEWTON",
      par = c(1, -1),
      step0 = 0.25
    )
  )

  expect_true(all(is.finite(result$par)))
  expect_identical(witness$counts(), c(hs = 1L, hi = 0L))
})

test_that("NEWTON retains fallback for symmetric non-positive Hessians", {
  hessians <- list(
    zero = matrix(0, nrow = 2, ncol = 2),
    indefinite = diag(c(1, -1))
  )

  for (case_name in names(hessians)) {
    witness <- curvature_validation_witness("hs", hessians[[case_name]])
    result <- do.call(
      mize,
      curvature_validation_args(
        witness$fg,
        "NEWTON",
        par = c(1, -2),
        step0 = 0.25
      )
    )

    expect_equal(result$par, c(0.75, -1.5), info = case_name)
    expect_identical(witness$counts(), c(hs = 1L, hi = 0L), info = case_name)
  }
})

test_that("finite Hessians with nonfinite directions fall back safely", {
  cases <- list(
    newton_zero_vector = list(method = "NEWTON", value = c(0, 0)),
    newton_tiny_vector = list(
      method = "NEWTON",
      value = c(1e-320, 1e-320)
    ),
    newton_tiny_matrix = list(
      method = "NEWTON",
      value = diag(c(1e-320, 1e-320))
    ),
    phess_tiny_matrix = list(
      method = "PHESS",
      value = diag(c(1e-320, 1e-320))
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    witness <- curvature_validation_witness("hs", case$value)
    result <- tryCatch(
      do.call(
        mize,
        curvature_validation_args(
          witness$fg,
          case$method,
          par = c(1, -2),
          step0 = 0.25
        )
      ),
      error = identity
    )

    expect_false(inherits(result, "error"), info = case_name)
    if (!inherits(result, "error")) {
      expect_true(all(is.finite(result$par)), info = case_name)
      expect_equal(
        as.numeric(result$par),
        c(0.75, -1.5),
        info = case_name
      )
    }
    expect_identical(
      witness$counts(),
      c(hs = 1L, hi = 0L),
      info = case_name
    )
  }
})

test_that("NEWTON prefers Hessian callbacks over inverse Hessians", {
  hs_calls <- 0L
  hi_calls <- 0L
  fg <- validation_fg()
  fg$gr <- function(x) x
  fg$hs <- function(x) {
    hs_calls <<- hs_calls + 1L
    diag(length(x))
  }
  fg$hi <- function(x) {
    hi_calls <<- hi_calls + 1L
    stop("inverse Hessian callback must not be called")
  }

  result <- do.call(
    mize,
    curvature_validation_args(fg, "NEWTON", par = c(1, -1))
  )

  expect_equal(result$par, c(0, 0))
  expect_identical(c(hs = hs_calls, hi = hi_calls), c(hs = 1L, hi = 0L))
})

test_that("inverse Hessian consumers accept vector and matrix forms", {
  inverse_hessians <- list(
    integer_vector = rep(1L, 4),
    double_vector = rep(1, 4),
    integer_matrix = diag(rep(1L, 4)),
    double_matrix = diag(rep(1, 4)),
    dimnamed_matrix = structure(
      diag(rep(1, 4)),
      dimnames = list(letters[1:4], letters[5:8])
    )
  )

  for (method in c("NEWTON", "BFGS", "SR1", "L-BFGS")) {
    for (case_name in names(inverse_hessians)) {
      witness <- curvature_validation_witness(
        "hi",
        inverse_hessians[[case_name]]
      )
      result <- do.call(
        mize,
        curvature_validation_args(witness$fg, method)
      )

      expect_equal(
        as.numeric(result$par),
        rep(0, 4),
        info = paste(method, case_name)
      )
      expect_identical(
        witness$counts(),
        c(hs = 0L, hi = 1L),
        info = paste(method, case_name)
      )
    }
  }
})

test_that("matrix inverse Hessians preserve vector parameters across steps", {
  for (method in c("NEWTON", "L-BFGS")) {
    seen <- list()
    record_par <- function(x) {
      expect_null(dim(x), info = method)
      seen[[length(seen) + 1L]] <<- x
      x
    }
    fg <- list(
      fn = function(x) sum(record_par(x)^2) / 2,
      gr = function(x) record_par(x),
      hi = function(x) {
        record_par(x)
        diag(length(x))
      }
    )

    result <- mize(
      c(2, -1),
      fg,
      method = method,
      line_search = "constant",
      step0 = 0.25,
      max_iter = 3,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )

    expect_gt(length(seen), 3L)
    expect_true(all(vapply(seen, function(x) is.null(dim(x)), logical(1))))
    expect_null(dim(result$par), info = method)
    expect_equal(result$par, c(0.84375, -0.421875), info = method)
  }
})

test_that("NEWTON rejects malformed inverse Hessians at first consumption", {
  malformed <- list(
    wrong_type = "invalid",
    short_vector = c(1, 2),
    long_vector = rep(1, 8),
    wrong_matrix = diag(2),
    nonsquare_matrix = matrix(1, nrow = 2, ncol = 3),
    nonfinite_vector = c(1, 1, Inf, 1),
    nonfinite_matrix = diag(c(1, 1, NaN, 1)),
    asymmetric_matrix = matrix(
      c(1, 0, 0, 0, 0.25, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1),
      nrow = 4
    )
  )
  patterns <- c(
    wrong_type = "fg\\$hi\\(par\\).*finite numeric vector.*length 4",
    short_vector = "fg\\$hi\\(par\\).*vector.*length 4",
    long_vector = "fg\\$hi\\(par\\).*vector.*length 4",
    wrong_matrix = "fg\\$hi\\(par\\).*4 x 4",
    nonsquare_matrix = "fg\\$hi\\(par\\).*4 x 4",
    nonfinite_vector = "fg\\$hi\\(par\\).*finite",
    nonfinite_matrix = "fg\\$hi\\(par\\).*finite",
    asymmetric_matrix = "fg\\$hi\\(par\\).*symmetric"
  )

  for (case_name in names(malformed)) {
    witness <- curvature_validation_witness("hi", malformed[[case_name]])

    expect_callback_validation_error(
      function() {
        do.call(mize, curvature_validation_args(witness$fg, "NEWTON"))
      },
      patterns[[case_name]],
      info = case_name
    )
    expect_identical(
      witness$counts(),
      c(hs = 0L, hi = 1L),
      info = case_name
    )
  }
})

test_that("every quasi-Newton consumer rejects asymmetric inverse Hessians", {
  asymmetric <- matrix(
    c(1, 0, 0, 0, 0.25, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1),
    nrow = 4
  )

  for (method in c("BFGS", "SR1", "L-BFGS")) {
    witness <- curvature_validation_witness("hi", asymmetric)

    expect_callback_validation_error(
      function() {
        do.call(mize, curvature_validation_args(witness$fg, method))
      },
      "fg\\$hi\\(par\\).*symmetric.*4 x 4",
      info = method
    )
    expect_identical(
      witness$counts(),
      c(hs = 0L, hi = 1L),
      info = method
    )
  }
})

test_that("inverse Hessian symmetry allows floating-point noise", {
  noise <- 10 * .Machine$double.eps
  inverse_hessian <- matrix(c(1, 0.25, 0.25 + noise, 1), nrow = 2)
  witness <- curvature_validation_witness("hi", inverse_hessian)

  result <- do.call(
    mize,
    curvature_validation_args(
      witness$fg,
      "NEWTON",
      par = c(1, -1),
      step0 = 0.25
    )
  )

  expect_true(all(is.finite(result$par)))
  expect_identical(witness$counts(), c(hs = 0L, hi = 1L))
})

test_that("inverse Hessian validation does not require positive definiteness", {
  inverse_hessians <- list(
    zero_matrix = matrix(0, nrow = 2, ncol = 2),
    negative_vector = c(-1, -1)
  )

  for (method in c("NEWTON", "BFGS", "SR1", "L-BFGS")) {
    for (case_name in names(inverse_hessians)) {
      witness <- curvature_validation_witness(
        "hi",
        inverse_hessians[[case_name]]
      )
      result <- do.call(
        mize,
        curvature_validation_args(
          witness$fg,
          method,
          par = c(1, -2),
          step0 = 0.25
        )
      )

      expect_equal(
        as.numeric(result$par),
        c(0.75, -1.5),
        info = paste(method, case_name)
      )
      expect_identical(
        witness$counts(),
        c(hs = 0L, hi = 1L),
        info = paste(method, case_name)
      )
    }
  }
})

test_that("finite inverse Hessians with overflowing directions fall back", {
  for (method in c("NEWTON", "BFGS", "L-BFGS", "SR1")) {
    witness <- curvature_validation_witness("hi", c(1e308, 1e308))
    result <- tryCatch(
      do.call(
        mize,
        curvature_validation_args(
          witness$fg,
          method,
          par = c(1, -2),
          step0 = 0.25
        )
      ),
      error = identity
    )

    expect_false(inherits(result, "error"), info = method)
    if (!inherits(result, "error")) {
      expect_true(all(is.finite(result$par)), info = method)
      expect_equal(
        as.numeric(result$par),
        c(0.75, -1.5),
        info = method
      )
    }
    expect_identical(
      witness$counts(),
      c(hs = 0L, hi = 1L),
      info = method
    )
  }
})

test_that("PHESS accepts full and repeated-block Hessians", {
  par <- c(1, -1, 2, -2)
  weights <- c(1, 2, 1, 2)
  fn <- function(x) sum(weights * x^2) / 2
  gr <- function(x) weights * x
  hessians <- list(
    integer_full = diag(as.integer(weights)),
    double_full = diag(as.double(weights)),
    integer_repeated_block = diag(c(1L, 2L)),
    double_repeated_block = diag(c(1, 2))
  )

  for (case_name in names(hessians)) {
    witness <- curvature_validation_witness(
      "hs",
      hessians[[case_name]],
      fn,
      gr
    )
    result <- do.call(
      mize,
      curvature_validation_args(witness$fg, "PHESS", par)
    )

    expect_equal(result$par, rep(0, length(par)), info = case_name)
    expect_identical(witness$counts(), c(hs = 1L, hi = 0L), info = case_name)
  }
})

test_that("PHESS rejects malformed Hessians at first consumption", {
  malformed <- list(
    vector = rep(1, 4),
    nonsquare_matrix = matrix(1, nrow = 2, ncol = 3),
    nondividing_matrix = diag(3),
    nonfinite_matrix = diag(c(1, 1, Inf, 1)),
    asymmetric_matrix = matrix(
      c(2, 0, 0, 0, 0.25, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2),
      nrow = 4
    )
  )
  patterns <- c(
    vector = "fg\\$hs\\(par\\).*matrix.*vectors are not supported by PHESS",
    nonsquare_matrix = "fg\\$hs\\(par\\).*symmetric square numeric matrix",
    nondividing_matrix = "fg\\$hs\\(par\\).*dimension.*divides 4",
    nonfinite_matrix = "fg\\$hs\\(par\\).*finite",
    asymmetric_matrix = "fg\\$hs\\(par\\).*symmetric"
  )

  for (case_name in names(malformed)) {
    witness <- curvature_validation_witness("hs", malformed[[case_name]])

    expect_callback_validation_error(
      function() {
        do.call(mize, curvature_validation_args(witness$fg, "PHESS"))
      },
      patterns[[case_name]],
      info = case_name
    )
    expect_identical(
      witness$counts(),
      c(hs = 1L, hi = 0L),
      info = case_name
    )
  }
})

test_that("periodic PHESS validates each refreshed Hessian once", {
  # hessian_every is private and make_mize() uses the initial-only default.
  # This direct probe protects the otherwise unreachable refresh boundary.
  asymmetric <- matrix(c(1, 0, 0.25, 1), nrow = 2)
  witness <- curvature_validation_witness("hs", asymmetric)
  direction <- partial_hessian_direction(hessian_every = 1)
  opt <- list(cache = list(gr_curr = c(1, -1), gr_curr_iter = 1))

  expect_callback_validation_error(
    function() {
      direction$calculate(
        opt,
        stage = list(),
        sub_stage = direction,
        par = c(1, -1),
        fg = witness$fg,
        iter = 1
      )
    },
    "fg\\$hs\\(par\\).*symmetric"
  )
  expect_identical(witness$counts(), c(hs = 1L, hi = 0L))
})

test_that("initial parameters are validated before callbacks", {
  bad_pars <- list(
    character = "bad",
    list = list(1, 2),
    empty = numeric(),
    matrix = matrix(c(1, 2), nrow = 1),
    missing = c(1, NA_real_),
    nan = c(1, NaN),
    infinite = c(1, Inf)
  )
  apis <- list(
    mize = function(par, fg) {
      mize(
        par,
        fg,
        method = "SD",
        line_search = "constant",
        step0 = 1,
        max_iter = 0,
        check_conv_every = NULL
      )
    },
    mize_init = function(par, fg) {
      opt <- make_mize(method = "SD", line_search = "constant", step0 = 1)
      mize_init(opt, par, fg)
    },
    make_mize = function(par, fg) {
      make_mize(
        method = "SD",
        line_search = "constant",
        step0 = 1,
        par = par,
        fg = fg
      )
    }
  )

  for (api_name in names(apis)) {
    for (par_name in names(bad_pars)) {
      calls <- new.env(parent = emptyenv())
      calls$fn <- 0
      calls$gr <- 0
      fg <- validation_fg(calls)

      expect_error(
        apis[[api_name]](bad_pars[[par_name]], fg),
        "par must",
        info = paste(api_name, par_name)
      )
      expect_identical(
        c(calls$fn, calls$gr),
        c(0, 0),
        info = paste(api_name, par_name)
      )
    }
  }
})

test_that("valid integer and double parameters retain names", {
  parameters <- list(
    integer = c(left = 1L, right = -1L),
    double = c(left = 1, right = -1)
  )

  for (par in parameters) {
    calls <- new.env(parent = emptyenv())
    calls$fn <- 0
    calls$gr <- 0
    res <- mize(
      par,
      validation_fg(calls),
      method = "SD",
      line_search = "constant",
      step0 = 0,
      max_iter = 0,
      max_fn = 0,
      max_gr = 0,
      max_fg = 0,
      check_conv_every = NULL
    )

    expect_identical(res$par, par)
    expect_identical(c(calls$fn, calls$gr), c(0, 0))
  }
})

test_that("iteration and evaluation limits require integer controls", {
  invalid_max_iter <- list(
    wrong_type = "1",
    nonscalar = c(1, 2),
    missing = NA_real_,
    nan = NaN,
    fractional = 1.5,
    negative = -1,
    negative_infinite = -Inf
  )
  invalid_budgets <- invalid_max_iter

  for (case_name in names(invalid_max_iter)) {
    expect_error(
      make_mize(max_iter = invalid_max_iter[[case_name]]),
      "max_iter",
      info = case_name
    )
  }
  for (argument in c("max_fn", "max_gr", "max_fg")) {
    for (case_name in names(invalid_budgets)) {
      args <- setNames(list(invalid_budgets[[case_name]]), argument)
      expect_error(
        do.call(make_mize, args),
        argument,
        info = paste(argument, case_name)
      )
    }
  }

  expect_no_error(make_mize(max_iter = 0))
  expect_no_error(make_mize(max_iter = Inf))
  expect_no_error(make_mize(max_fn = 0, max_gr = 0, max_fg = 0))
  expect_no_error(make_mize(max_fn = Inf, max_gr = Inf, max_fg = Inf))
})

test_that("line-search evaluation limits are validated when consumed", {
  invalid_budgets <- list(
    wrong_type = "1",
    nonscalar = c(1, 2),
    missing = NA_real_,
    nan = NaN,
    fractional = 1.5,
    negative = -1,
    negative_infinite = -Inf
  )

  for (argument in c("ls_max_fn", "ls_max_gr", "ls_max_fg")) {
    for (case_name in names(invalid_budgets)) {
      args <- c(
        list(line_search = "more-thuente"),
        setNames(list(invalid_budgets[[case_name]]), argument)
      )
      expect_error(
        do.call(make_mize, args),
        argument,
        info = paste(argument, case_name)
      )
    }
  }

  expect_no_error(make_mize(
    line_search = "more-thuente",
    ls_max_fn = 0,
    ls_max_gr = Inf,
    ls_max_fg = 0
  ))
  expect_no_error(make_mize(
    line_search = "constant",
    step0 = 1,
    ls_max_fn = "ignored",
    ls_max_gr = NA,
    ls_max_fg = -1
  ))
})

test_that("cadences are positive integers when consumed", {
  invalid_cadences <- list(
    wrong_type = "1",
    nonscalar = c(1, 2),
    missing = NA_real_,
    nan = NaN,
    fractional = 1.5,
    zero = 0,
    negative = -1,
    infinite = Inf
  )

  for (argument in c("check_conv_every", "log_every")) {
    for (case_name in names(invalid_cadences)) {
      calls <- new.env(parent = emptyenv())
      calls$fn <- 0
      calls$gr <- 0
      args <- list(
        par = c(1),
        fg = validation_fg(calls),
        method = "SD",
        line_search = "constant",
        step0 = 1,
        max_iter = 0,
        store_progress = argument == "log_every"
      )
      args[[argument]] <- invalid_cadences[[case_name]]

      expect_error(
        do.call(mize, args),
        argument,
        info = paste(argument, case_name)
      )
      expect_identical(
        c(calls$fn, calls$gr),
        c(0, 0),
        info = paste(argument, case_name)
      )
    }
  }

  expect_error(
    mize(
      c(1),
      validation_fg(),
      max_iter = 0,
      check_conv_every = NULL,
      store_progress = TRUE
    ),
    "check_conv_every"
  )
  expect_error(
    mize(
      c(1),
      validation_fg(),
      max_iter = 0,
      check_conv_every = NULL,
      verbose = TRUE
    ),
    "check_conv_every"
  )
  expect_no_error(mize(
    c(1),
    validation_fg(),
    method = "SD",
    line_search = "constant",
    step0 = 1,
    max_iter = 0,
    max_fn = 0,
    max_gr = 0,
    check_conv_every = NULL,
    log_every = "ignored"
  ))
  expect_no_error(mize(
    c(1),
    validation_fg(),
    method = "SD",
    line_search = "constant",
    step0 = 1,
    max_iter = 0,
    max_fn = 0,
    max_gr = 0,
    check_conv_every = 2,
    log_every = 4,
    store_progress = TRUE
  ))
})

test_that("convergence tolerances are null or finite nonnegative scalars", {
  invalid_tolerances <- list(
    wrong_type = "0",
    nonscalar = c(0, 1),
    missing = NA_real_,
    nan = NaN,
    negative = -1,
    infinite = Inf
  )

  for (argument in c(
    "abs_tol",
    "rel_tol",
    "grad_tol",
    "ginf_tol",
    "step_tol"
  )) {
    for (case_name in names(invalid_tolerances)) {
      args <- setNames(list(invalid_tolerances[[case_name]]), argument)
      expect_error(
        do.call(make_mize, args),
        argument,
        info = paste(argument, case_name)
      )
    }
  }

  expect_no_error(make_mize(
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = 0,
    ginf_tol = 0,
    step_tol = 0
  ))

  calls <- new.env(parent = emptyenv())
  calls$fn <- 0
  calls$gr <- 0
  expect_error(
    mize_init(
      make_mize(),
      c(1),
      validation_fg(calls),
      grad_tol = Inf
    ),
    "grad_tol"
  )
  expect_identical(c(calls$fn, calls$gr), c(0, 0))
})

test_that("range-checked numeric controls reject malformed scalars", {
  cases <- list(
    list(args = list(method = "NAG", nest_q = "bad"), argument = "nest_q"),
    list(args = list(method = "NAG", nest_q = c(0, 1)), argument = "nest_q"),
    list(args = list(method = "NAG", nest_q = NA_real_), argument = "nest_q"),
    list(args = list(method = "NAG", nest_q = NaN), argument = "nest_q"),
    list(args = list(method = "NAG", nest_q = Inf), argument = "nest_q"),
    list(args = list(method = "NAG", nest_q = -0.1), argument = "nest_q"),
    list(args = list(method = "DBD", step_up = "bad"), argument = "step_up"),
    list(args = list(method = "DBD", step_up = c(1, 2)), argument = "step_up"),
    list(args = list(method = "DBD", step_up = NA_real_), argument = "step_up"),
    list(args = list(method = "DBD", step_up = NaN), argument = "step_up"),
    list(args = list(method = "DBD", step_up = Inf), argument = "step_up"),
    list(args = list(method = "DBD", step_up = 0), argument = "step_up"),
    list(
      args = list(method = "DBD", step_down = "bad"),
      argument = "step_down"
    ),
    list(
      args = list(method = "DBD", step_down = c(0.1, 0.2)),
      argument = "step_down"
    ),
    list(
      args = list(method = "DBD", step_down = NA_real_),
      argument = "step_down"
    ),
    list(args = list(method = "DBD", step_down = NaN), argument = "step_down"),
    list(args = list(method = "DBD", step_down = Inf), argument = "step_down"),
    list(args = list(method = "DBD", step_down = -0.1), argument = "step_down"),
    list(
      args = list(method = "DBD", dbd_weight = "bad"),
      argument = "dbd_weight"
    ),
    list(
      args = list(method = "DBD", dbd_weight = c(0, 1)),
      argument = "dbd_weight"
    ),
    list(
      args = list(method = "DBD", dbd_weight = NA_real_),
      argument = "dbd_weight"
    ),
    list(
      args = list(method = "DBD", dbd_weight = NaN),
      argument = "dbd_weight"
    ),
    list(
      args = list(method = "DBD", dbd_weight = Inf),
      argument = "dbd_weight"
    ),
    list(
      args = list(method = "DBD", dbd_weight = 1.1),
      argument = "dbd_weight"
    ),
    list(args = list(line_search = "mt", c1 = "bad"), argument = "c1"),
    list(args = list(line_search = "mt", c1 = c(0.1, 0.2)), argument = "c1"),
    list(args = list(line_search = "mt", c1 = NA_real_), argument = "c1"),
    list(args = list(line_search = "mt", c1 = NaN), argument = "c1"),
    list(args = list(line_search = "mt", c1 = Inf), argument = "c1"),
    list(args = list(line_search = "mt", c1 = 0), argument = "c1"),
    list(
      args = list(line_search = "mt", c1 = 0.1, c2 = "bad"),
      argument = "c2"
    ),
    list(
      args = list(line_search = "mt", c1 = 0.1, c2 = c(0.2, 0.3)),
      argument = "c2"
    ),
    list(
      args = list(line_search = "mt", c1 = 0.1, c2 = NA_real_),
      argument = "c2"
    ),
    list(args = list(line_search = "mt", c1 = 0.1, c2 = NaN), argument = "c2"),
    list(args = list(line_search = "mt", c1 = 0.1, c2 = Inf), argument = "c2"),
    list(args = list(line_search = "mt", c1 = 0.1, c2 = 0.1), argument = "c2"),
    list(
      args = list(line_search = "mt", ls_max_alpha_mult = "bad"),
      argument = "ls_max_alpha_mult"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha_mult = c(1, 2)),
      argument = "ls_max_alpha_mult"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha_mult = NA_real_),
      argument = "ls_max_alpha_mult"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha_mult = NaN),
      argument = "ls_max_alpha_mult"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha_mult = 0),
      argument = "ls_max_alpha_mult"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha = "bad"),
      argument = "ls_max_alpha"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha = c(1, 2)),
      argument = "ls_max_alpha"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha = NA_real_),
      argument = "ls_max_alpha"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha = NaN),
      argument = "ls_max_alpha"
    ),
    list(
      args = list(line_search = "mt", ls_max_alpha = 0),
      argument = "ls_max_alpha"
    ),
    list(
      args = list(line_search = "mt", step_next_init = c(1, 2)),
      argument = "step_next_init"
    ),
    list(
      args = list(line_search = "mt", step_next_init = NA_real_),
      argument = "step_next_init"
    ),
    list(
      args = list(line_search = "mt", step_next_init = NaN),
      argument = "step_next_init"
    ),
    list(
      args = list(line_search = "mt", step_next_init = Inf),
      argument = "step_next_init"
    ),
    list(
      args = list(line_search = "mt", step_next_init = 0),
      argument = "step_next_init"
    ),
    list(
      args = list(line_search = "mt", strong_curvature = "bad"),
      argument = "strong_curvature"
    ),
    list(
      args = list(line_search = "mt", strong_curvature = NA),
      argument = "strong_curvature"
    ),
    list(
      args = list(line_search = "mt", approx_armijo = c(TRUE, FALSE)),
      argument = "approx_armijo"
    ),
    list(
      args = list(line_search = "mt", approx_armijo = NA),
      argument = "approx_armijo"
    ),
    list(
      args = list(line_search = "mt", ls_safe_cubic = 1),
      argument = "ls_safe_cubic"
    ),
    list(
      args = list(line_search = "mt", ls_safe_cubic = NA),
      argument = "ls_safe_cubic"
    )
  )

  for (case in cases) {
    expect_error(
      do.call(make_mize, case$args),
      case$argument,
      info = paste(names(case$args), collapse = ", ")
    )
  }

  expect_no_error(make_mize(method = "NAG", nest_q = 0))
  expect_no_error(make_mize(method = "NAG", nest_q = 1))
  expect_no_error(make_mize(method = "DBD", step_down = 0, dbd_weight = 1))
  expect_no_error(make_mize(
    line_search = "mt",
    ls_max_alpha_mult = Inf,
    ls_max_alpha = Inf,
    step_next_init = 0.5
  ))
})

test_that("method-specific controls are validated only when relevant", {
  invalid_preconditioners <- list(
    unsupported = "diagonal",
    wrong_type = 1,
    nonscalar = c("", "l-bfgs"),
    missing = NA_character_
  )
  for (method in c("CG", "TN")) {
    for (case_name in names(invalid_preconditioners)) {
      expect_error(
        make_mize(
          method = method,
          preconditioner = invalid_preconditioners[[case_name]]
        ),
        "preconditioner",
        info = paste(method, case_name)
      )
    }
    expect_no_error(make_mize(method = method, preconditioner = ""))
    expect_no_error(make_mize(method = method, preconditioner = "L-BfGs"))
  }

  invalid_memory <- list(
    wrong_type = "1",
    nonscalar = c(1, 2),
    missing = NA_real_,
    nan = NaN,
    fractional = 1.5,
    zero = 0,
    negative = -1,
    infinite = Inf
  )
  memory_configurations <- list(
    list(method = "L-BFGS"),
    list(method = "CG", preconditioner = "L-BFGS"),
    list(method = "TN", preconditioner = "L-BFGS")
  )
  for (configuration in memory_configurations) {
    for (case_name in names(invalid_memory)) {
      expect_error(
        do.call(
          make_mize,
          c(configuration, list(memory = invalid_memory[[case_name]]))
        ),
        "memory",
        info = paste(configuration$method, case_name)
      )
    }
  }

  expect_no_error(make_mize(method = "CG", preconditioner = "", memory = 0))
  expect_no_error(make_mize(
    method = "SD",
    preconditioner = c("bad", "worse"),
    memory = NA_real_,
    tn_init = "bad"
  ))

  for (tn_init in list(0, 0L, "previous", "prev", "PrEvIoUs", "PrEv")) {
    expect_no_error(make_mize(method = "TN", tn_init = tn_init))
  }
  for (tn_init in list(1, NULL, "bad", "l-bfgs", c("prev", "previous"))) {
    expect_error(make_mize(method = "TN", tn_init = tn_init), "tn_init")
  }
})

test_that("feature-specific iteration controls are gated by consumption", {
  invalid_nonnegative_counts <- list(
    wrong_type = "1",
    nonscalar = c(1, 2),
    missing = NA_real_,
    nan = NaN,
    fractional = 1.5,
    negative = -1,
    infinite = Inf
  )

  for (case_name in names(invalid_nonnegative_counts)) {
    expect_error(
      make_mize(
        method = "NAG",
        nest_burn_in = invalid_nonnegative_counts[[case_name]]
      ),
      "nest_burn_in",
      info = case_name
    )
    expect_error(
      make_mize(
        mom_schedule = "switch",
        mom_switch_iter = invalid_nonnegative_counts[[case_name]]
      ),
      "mom_switch_iter",
      info = case_name
    )
  }

  expect_no_error(make_mize(method = "NAG", nest_burn_in = 0))
  expect_no_error(make_mize(
    mom_schedule = "switch",
    mom_init = 0.2,
    mom_final = 0.8,
    mom_switch_iter = 0
  ))
  expect_error(
    make_mize(mom_schedule = "switch", mom_switch_iter = NULL),
    "mom_switch_iter"
  )
  expect_no_error(make_mize(
    method = "SD",
    nest_q = "ignored",
    nest_burn_in = "ignored",
    mom_switch_iter = "ignored",
    restart_wait = "ignored",
    step_up = "ignored",
    step_down = "ignored",
    dbd_weight = "ignored",
    line_search = "constant",
    step0 = 1,
    c1 = "ignored",
    c2 = "ignored",
    ls_max_alpha_mult = "ignored",
    ls_max_alpha = "ignored",
    step_next_init = "ignored"
  ))
  expect_no_error(make_mize(
    method = "NAG",
    nest_convex_approx = TRUE,
    nest_q = "ignored"
  ))

  for (restart_wait in list(0, -1, 1.5, Inf, NA_real_, NaN, "1", c(1, 2))) {
    expect_error(
      make_mize(method = "NAG", restart = "fn", restart_wait = restart_wait),
      "restart_wait"
    )
  }
  expect_no_error(make_mize(method = "NAG", restart = "fn", restart_wait = 1))
})

test_that("momentum controls require finite scalar values when consumed", {
  invalid_scalars <- list(
    nonscalar = c(0.2, 0.8),
    missing = NA_real_,
    nan = NaN,
    positive_infinite = Inf,
    negative_infinite = -Inf
  )

  for (case_name in names(invalid_scalars)) {
    expect_error(
      make_mize(
        method = "MOM",
        mom_schedule = invalid_scalars[[case_name]]
      ),
      "mom_schedule",
      info = case_name
    )
  }

  for (schedule in c("ramp", "switch")) {
    schedule_args <- list(
      method = "MOM",
      mom_schedule = schedule,
      mom_init = 0.2,
      mom_final = 0.8
    )
    if (schedule == "switch") {
      schedule_args$mom_switch_iter <- 2
    }

    expect_error(
      do.call(
        make_mize,
        utils::modifyList(schedule_args, list(mom_init = NULL))
      ),
      "mom_init",
      info = paste(schedule, "missing mom_init")
    )
    expect_error(
      do.call(
        make_mize,
        utils::modifyList(schedule_args, list(mom_final = NULL))
      ),
      "mom_final",
      info = paste(schedule, "missing mom_final")
    )

    for (argument in c("mom_init", "mom_final")) {
      for (case_name in names(invalid_scalars)) {
        args <- schedule_args
        args[[argument]] <- invalid_scalars[[case_name]]
        expect_error(
          do.call(make_mize, args),
          argument,
          info = paste(schedule, argument, case_name)
        )
      }
    }
    for (argument in c("mom_init", "mom_final")) {
      args <- schedule_args
      args[[argument]] <- "bad"
      expect_error(
        do.call(make_mize, args),
        argument,
        info = paste(schedule, argument, "wrong_type")
      )
    }
  }

  invalid_flags <- list(
    wrong_type = 1,
    nonscalar = c(TRUE, FALSE),
    missing = NA
  )
  for (case_name in names(invalid_flags)) {
    expect_error(
      make_mize(method = "MOM", use_init_mom = invalid_flags[[case_name]]),
      "use_init_mom",
      info = case_name
    )
    expect_error(
      make_mize(
        method = "MOM",
        mom_linear_weight = invalid_flags[[case_name]]
      ),
      "mom_linear_weight",
      info = case_name
    )
    expect_error(
      make_mize(
        method = "NAG",
        nest_convex_approx = invalid_flags[[case_name]]
      ),
      "nest_convex_approx",
      info = case_name
    )
    expect_error(
      make_mize(method = "DBD", norm_direction = invalid_flags[[case_name]]),
      "norm_direction",
      info = case_name
    )
  }

  expect_no_error(make_mize(method = "MOM", mom_schedule = -1))
  expect_no_error(make_mize(method = "MOM", mom_schedule = 2))
})

test_that("function-valued momentum schedules return finite numeric scalars", {
  fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )
  invalid_results <- list(
    wrong_type = "bad",
    nonscalar = c(0.2, 0.8),
    missing = NA_real_,
    nan = NaN,
    positive_infinite = Inf,
    negative_infinite = -Inf
  )

  for (case_name in names(invalid_results)) {
    schedule <- local({
      value <- invalid_results[[case_name]]
      function(iter) value
    })
    expect_error(
      mize(
        c(1, -1),
        fg,
        method = "MOM",
        line_search = "constant",
        step0 = 0.1,
        mom_schedule = schedule,
        use_init_mom = TRUE,
        max_iter = 1
      ),
      "mom_schedule function result must be a finite numeric scalar",
      fixed = TRUE,
      info = case_name
    )
  }

  expect_no_error(mize(
    c(1, -1),
    fg,
    method = "MOM",
    line_search = "constant",
    step0 = 0.1,
    mom_schedule = function(iter) 0.5,
    use_init_mom = TRUE,
    max_iter = 1
  ))
})

test_that("loop and summary observation flags reject malformed values", {
  fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )
  invalid_flags <- list(
    wrong_type = 1,
    nonscalar = c(TRUE, FALSE),
    missing = NA
  )

  for (argument in c("verbose", "store_progress")) {
    for (case_name in names(invalid_flags)) {
      args <- list(
        par = c(1, -1),
        fg = fg,
        method = "SD",
        line_search = "constant",
        step0 = 0.1,
        max_iter = 0
      )
      args[[argument]] <- invalid_flags[[case_name]]
      expect_error(
        do.call(mize, args),
        argument,
        info = paste(argument, case_name)
      )
    }
  }

  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1, -1),
    fg = fg
  )
  for (argument in c("calc_fn", "calc_gr")) {
    for (case_name in names(invalid_flags)) {
      args <- list(opt = opt, par = c(1, -1), fg = fg)
      args[[argument]] <- invalid_flags[[case_name]]
      expect_error(
        do.call(mize_step_summary, args),
        argument,
        info = paste(argument, case_name)
      )
    }
    for (value in list(NULL, FALSE, TRUE)) {
      args <- list(opt = opt, par = c(1, -1), fg = fg)
      args[[argument]] <- value
      expect_no_error(do.call(mize_step_summary, args))
    }
  }
})

test_that("constant line search requires a finite numeric scalar step0", {
  invalid_steps <- list(
    missing = NULL,
    wrong_type = "one",
    nonscalar = c(1, 2),
    missing_value = NA_real_,
    nan = NaN,
    positive_infinite = Inf,
    negative_infinite = -Inf
  )

  for (case_name in names(invalid_steps)) {
    expect_error(
      make_mize(line_search = "constant", step0 = invalid_steps[[case_name]]),
      "step0",
      info = case_name
    )
  }
  for (step0 in c(-1, 0, 0.25)) {
    expect_no_error(make_mize(line_search = "CoNsTaNt", step0 = step0))
  }

  expect_no_error(make_mize(line_search = "mt", step0 = "RaSmUsSeN"))
  expect_no_error(make_mize(method = "DBD", step0 = "RaSmUsSeN"))
})

test_that("DBD requires a positive scalar or documented step0 initializer", {
  invalid_steps <- list(
    zero = 0,
    negative = -1,
    nonscalar = c(0.1, 0.2),
    missing_value = NA_real_,
    nan = NaN,
    positive_infinite = Inf,
    negative_infinite = -Inf,
    unknown = "bad",
    character_vector = c("rasmussen", "scipy")
  )

  for (case_name in names(invalid_steps)) {
    expect_error(
      make_mize(method = "DBD", step0 = invalid_steps[[case_name]]),
      "step0",
      info = case_name
    )
  }

  fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )
  for (step0 in c("RaSmUsSeN", "ScIpY", "ScHmIdT", "HZ", "Hager-Zhang")) {
    result <- mize(
      c(1, -1),
      fg,
      method = "DBD",
      step0 = step0,
      max_iter = 1,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL,
      ginf_tol = NULL,
      step_tol = NULL
    )
    expect_equal(length(result$par), 2L, info = step0)
    expect_true(all(is.finite(result$par)), info = step0)
  }
  expect_no_error(make_mize(method = "DBD", step0 = 0.1))
})

test_that("DBD safeguards a non-finite derived step0", {
  fg <- list(
    fn = function(x) 0,
    gr = function(x) rep(1e308, length(x))
  )

  result <- mize(
    c(1, -1),
    fg,
    method = "DBD",
    step0 = "scipy",
    max_iter = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_true(all(is.finite(result$par)))
})

test_that("DBD Hager-Zhang step0 uses the current parameter scale", {
  fg <- list(
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x
  )

  result <- mize(
    c(1, -1),
    fg,
    method = "DBD",
    step0 = "hz",
    max_iter = 1,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL,
    ginf_tol = NULL,
    step_tol = NULL
  )

  expect_equal(result$par, c(0.989, -0.989))
})

test_that("Wolfe line searches require a positive finite numeric step0", {
  searches <- c("more-thuente", "rasmussen", "schmidt", "hager-zhang")
  invalid_steps <- list(
    zero = 0,
    negative = -1,
    nonscalar = c(1, 2),
    missing_value = NA_real_,
    nan = NaN,
    positive_infinite = Inf,
    negative_infinite = -Inf
  )

  for (line_search in searches) {
    for (case_name in names(invalid_steps)) {
      expect_error(
        make_mize(
          line_search = line_search,
          step0 = invalid_steps[[case_name]]
        ),
        "step0",
        info = paste(line_search, case_name)
      )
    }
    expect_no_error(make_mize(line_search = line_search, step0 = 0.25))
    expect_no_error(make_mize(line_search = line_search, step0 = "RaSmUsSeN"))
  }

  witness <- validation_result_witness()
  expect_error(
    mize(
      par = c(1, -1),
      fg = witness$fg,
      line_search = "rasmussen",
      step0 = Inf
    ),
    "step0"
  )
  expect_identical(witness$counts(), c(fn = 0L, gr = 0L, fg = 0L))
})

test_that("mize_init validates effective convergence controls before callbacks", {
  cases <- list(
    list(argument = "max_fn", value = 1.5),
    list(argument = "max_gr", value = NA_real_),
    list(argument = "max_fg", value = -1),
    list(argument = "abs_tol", value = Inf)
  )

  for (case in cases) {
    calls <- new.env(parent = emptyenv())
    calls$fn <- 0
    calls$gr <- 0
    args <- list(
      opt = make_mize(),
      par = c(1),
      fg = validation_fg(calls)
    )
    args[[case$argument]] <- case$value

    expect_error(
      do.call(mize_init, args),
      case$argument,
      info = case$argument
    )
    expect_identical(c(calls$fn, calls$gr), c(0, 0), info = case$argument)
  }
})

test_that("stateful max_iter accepts explicit and default Inf", {
  calls <- new.env(parent = emptyenv())
  calls$fn <- 0
  calls$gr <- 0
  fg <- validation_fg(calls)

  opt <- make_mize(max_iter = 12)
  initialized <- mize_init(opt, c(1), fg, max_iter = Inf)
  expect_identical(initialized$convergence$max_iter, Inf)
  expect_identical(c(calls$fn, calls$gr), c(0, 0))

  custom <- make_opt(make_stages(
    gradient_stage(
      direction = sd_direction(),
      step_size = constant_step_size(1)
    )
  ))
  custom <- mize_init(custom, c(1), fg)
  expect_identical(custom$convergence$max_iter, Inf)
  expect_true(custom$is_initialized)
  expect_identical(c(calls$fn, calls$gr), c(0, 0))
})

test_that("mize rejects invalid configuration before callbacks", {
  invalid_max_iter <- list(
    wrong_type = "1",
    nonscalar = c(1, 2),
    missing = NA_real_,
    nan = NaN,
    fractional = 1.5,
    negative = -1,
    positive_infinite = Inf,
    negative_infinite = -Inf
  )
  cases <- c(
    lapply(invalid_max_iter, function(value) {
      list(argument = "max_iter", value = value)
    }),
    list(
      list(argument = "max_fn", value = NA_real_),
      list(argument = "abs_tol", value = Inf),
      list(argument = "check_conv_every", value = 0)
    )
  )

  for (case in cases) {
    calls <- new.env(parent = emptyenv())
    calls$fn <- 0
    calls$gr <- 0
    args <- list(
      par = c(1),
      fg = validation_fg(calls),
      method = "SD",
      line_search = "constant",
      step0 = 1
    )
    args[[case$argument]] <- case$value

    expect_error(do.call(mize, args), case$argument, info = case$argument)
    expect_identical(c(calls$fn, calls$gr), c(0, 0), info = case$argument)
  }
})
