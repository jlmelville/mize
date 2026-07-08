make_spd_quadratic <- function() {
  a <- matrix(
    c(
      4,
      1,
      0,
      1,
      3,
      0.5,
      0,
      0.5,
      2
    ),
    nrow = 3
  )
  b <- c(1, -2, 0.5)

  list(
    fn = function(x) {
      drop(0.5 * t(x) %*% a %*% x - t(b) %*% x)
    },
    gr = function(x) {
      as.vector(a %*% x - b)
    },
    hs = function(x) {
      a
    },
    hi = function(x) {
      solve(a)
    },
    fg = function(x) {
      list(
        fn = drop(0.5 * t(x) %*% a %*% x - t(b) %*% x),
        gr = as.vector(a %*% x - b)
      )
    },
    a = a,
    b = b,
    xmin = solve(a, b)
  )
}

exact_quadratic_step_size <- function(a) {
  # Test-only exact line search for the quadratic CG invariant.
  make_step_size(list(
    name = "exact_quadratic",
    calculate = function(opt, stage, sub_stage, par, fg, iter) {
      pm <- stage$direction$value
      gm <- opt$cache$gr_curr
      sub_stage$value <- -dot(gm, pm) / drop(t(pm) %*% a %*% pm)
      list(sub_stage = sub_stage)
    }
  ))
}

cg_update_variants <- function() {
  list(
    fr = fr_update,
    cd = cd_update,
    dy = dy_update,
    hs = hs_update,
    "hs+" = hs_plus_update,
    pr = pr_update,
    "pr+" = pr_plus_update,
    ls = ls_update,
    hz = hz_update,
    "hz+" = hz_plus_update,
    prfr = prfr_update
  )
}

exact_quadratic_cg_trace <- function(cg_update, par, fg) {
  opt <- make_opt(
    make_stages(
      gradient_stage(
        direction = cg_direction(cg_update = cg_update),
        step_size = exact_quadratic_step_size(fg$a)
      ),
      verbose = FALSE
    )
  )
  opt$count_res_fg <- FALSE
  opt <- mize_init(
    opt,
    par,
    fg,
    max_iter = length(par),
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  par_trace <- matrix(NA_real_, nrow = length(par) + 1, ncol = length(par))
  direction_trace <- matrix(NA_real_, nrow = length(par), ncol = length(par))
  alpha <- rep(NA_real_, length(par))
  par_trace[1, ] <- par

  for (iter in seq_along(par)) {
    step <- mize_step(opt, par, fg)
    opt <- step$opt
    par <- as.numeric(step$par)

    direction_trace[iter, ] <- as.numeric(
      opt$stages$gradient_descent$direction$value
    )
    alpha[iter] <- opt$stages$gradient_descent$step_size$value
    par_trace[iter + 1, ] <- par
  }

  list(
    iter = opt$iter,
    par = par,
    gradient = fg$gr(par),
    par_trace = par_trace,
    direction_trace = direction_trace,
    alpha = alpha
  )
}

expect_fg_consistent <- function(fg, par) {
  combined <- fg$fg(par)

  expect_equal(combined$fn, fg$fn(par), tolerance = 1e-12)
  expect_equal(combined$gr, fg$gr(par), tolerance = 1e-12)
}

test_that("Newton reaches SPD quadratic minimizer in one exact step", {
  quadratic_fg <- make_spd_quadratic()
  par <- c(3, -1, 2)

  hessian_variants <- list(
    hs = within(quadratic_fg, hi <- NULL),
    hi = within(quadratic_fg, hs <- NULL)
  )

  for (fg in hessian_variants) {
    res <- mize(
      par,
      fg,
      method = "NEWTON",
      line_search = "constant",
      step0 = 1,
      max_iter = 1,
      check_conv_every = NULL,
      abs_tol = NULL,
      rel_tol = NULL,
      grad_tol = NULL
    )

    expect_equal(as.numeric(res$par), quadratic_fg$xmin, tolerance = 1e-12)
    expect_equal(fg$gr(as.numeric(res$par)), rep(0, 3), tolerance = 1e-12)
  }
})

test_that("PHESS reuses the initial Hessian on SPD quadratics", {
  quadratic_fg <- make_spd_quadratic()
  par <- c(3, -1, 2)
  hs_calls <- 0

  phess_fg <- quadratic_fg
  phess_fg$hi <- NULL
  phess_fg$hs <- function(x) {
    hs_calls <<- hs_calls + 1
    quadratic_fg$a
  }

  res <- mize(
    par,
    phess_fg,
    method = "PHESS",
    line_search = "constant",
    step0 = 1,
    max_iter = 2,
    check_conv_every = NULL,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  expect_equal(as.numeric(res$par), quadratic_fg$xmin, tolerance = 1e-12)
  expect_equal(phess_fg$gr(as.numeric(res$par)), rep(0, 3), tolerance = 1e-12)
  expect_equal(hs_calls, 1)
})

test_that("CG updates reach SPD quadratic minimizer under exact line search", {
  quadratic_fg <- make_spd_quadratic()
  par <- c(3, -1, 2)
  updates <- cg_update_variants()
  reference_trace <- exact_quadratic_cg_trace(updates$fr, par, quadratic_fg)

  for (update_name in names(updates)) {
    trace <- exact_quadratic_cg_trace(updates[[update_name]], par, quadratic_fg)

    expect_equal(trace$iter, length(par), info = update_name)
    expect_equal(
      trace$par,
      quadratic_fg$xmin,
      tolerance = 1e-12,
      info = update_name
    )
    expect_equal(
      trace$gradient,
      rep(0, 3),
      tolerance = 1e-11,
      info = update_name
    )
    expect_equal(
      trace$par_trace,
      reference_trace$par_trace,
      tolerance = 1e-12,
      info = update_name
    )
    expect_equal(
      trace$direction_trace,
      reference_trace$direction_trace,
      tolerance = 1e-11,
      info = update_name
    )
    expect_equal(
      trace$alpha,
      reference_trace$alpha,
      tolerance = 1e-12,
      info = update_name
    )
  }
})

test_that("BFGS and L-BFGS directions are descent after accepted curvature", {
  quadratic_fg <- make_spd_quadratic()
  par0 <- c(3, -1, 2)

  for (method in c("BFGS", "L-BFGS")) {
    opt <- make_mize(
      method = method,
      line_search = "constant",
      step0 = 0.05,
      scale_hess = FALSE,
      memory = 3,
      par = par0,
      fg = quadratic_fg
    )

    first_step <- mize_step(opt, par0, quadratic_fg)
    par <- as.numeric(first_step$par)
    second_step <- mize_step(first_step$opt, par, quadratic_fg)
    direction <- as.numeric(
      second_step$opt$stages$gradient_descent$direction$value
    )

    expect_lt(dot(quadratic_fg$gr(par), direction), 0)
  }
})

test_that("combined fg functions match separate function and gradient", {
  expect_fg_consistent(rosenbrock_fg, c(-0.7, 0.8))
  expect_fg_consistent(make_spd_quadratic(), c(0.25, -1.5, 2))
})
