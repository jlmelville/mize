condition_quadratic_fg <- function(center = 10) {
  list(
    fn = function(x) {
      (x - center)^2
    },
    gr = function(x) {
      2 * (x - center)
    }
  )
}

condition_setup <- function(fg, x, pv = -fg$gr(x) / abs(fg$gr(x))) {
  step0 <- make_step0(fg, x, pv)
  list(
    pv = pv,
    step0 = step0,
    phi = make_phi_alpha(x, fg, pv, calc_gradient_default = TRUE)
  )
}

test_that("Schmidt Armijo controls reject invalid types and ranges", {
  expect_error(
    new_schmidt_armijo_search(armijo_constant = NA_real_),
    "armijo_constant"
  )
  expect_error(new_schmidt_armijo_search(step_down = "half"), "step_down")
  expect_error(new_schmidt_armijo_search(step_down = 2), "step_down")
  expect_error(new_schmidt_armijo_search(max_fn = -1), "max_fn")
  expect_error(
    new_schmidt_armijo_search(progress_tolerance = Inf),
    "progress_tolerance"
  )
})

test_that("Armijo backtracking accepts steps with sufficient decrease", {
  cases <- list(
    list(
      name = "f1 small initial step",
      fg = f1,
      x = 0,
      alpha = 1e-1,
      c1 = 0.001
    ),
    list(name = "f2 large initial step", fg = f2, x = 0, alpha = 1e1, c1 = 0.1),
    list(
      name = "f4 very large initial step",
      fg = f4,
      x = 0,
      alpha = 1e3,
      c1 = 0.001
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- new_schmidt_armijo_search(
      armijo_constant = case$c1,
      max_fn = 10000
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = case$alpha,
      pm = setup$pv
    )

    expect_true(res$step$alpha > 0, info = case$name)
    expect_true(
      armijo_ok_step(setup$step0, res$step, case$c1),
      info = case$name
    )
  }
})

test_that("Schmidt Armijo accepts internal phi results without parameters", {
  step0 <- list(alpha = 0, f = 1, df = -1, d = -1)
  phi <- function(alpha, calc_gradient = TRUE) {
    list(alpha = alpha, f = 0, df = -1, d = -1)
  }

  res <- new_schmidt_armijo_search(armijo_constant = 0.1, max_fn = 1)(
    phi = phi,
    step0 = step0,
    alpha = 0.5,
    pm = 1
  )

  expect_equal(res$step$alpha, 0.5)
  expect_equal(res$step$f, 0)
  expect_equal(res$nfn, 1)
  expect_equal(res$ngr, 1)
})

test_that("Schmidt cubic Armijo backs off a nonfinite objective", {
  evaluated <- numeric()
  step0 <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    if (alpha >= 1) {
      return(list(alpha = alpha, f = Inf, df = Inf, d = Inf, par = alpha))
    }
    list(alpha = alpha, f = 0, df = 0, d = 0, par = alpha)
  }

  result <- new_schmidt_armijo_search(armijo_constant = 0.1, max_fn = 2)(
    phi = phi,
    step0 = step0,
    alpha = 1,
    pm = 1
  )

  expect_equal(evaluated, c(1, 0.5))
  expect_equal(result$step$alpha, 0.5)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})

test_that("Schmidt cubic Armijo uses values when a trial gradient is nonfinite", {
  evaluated <- numeric()
  step0 <- list(alpha = 0, f = 1, df = -1, d = -1, par = 0)
  phi <- function(alpha, calc_gradient = TRUE) {
    evaluated <<- c(evaluated, alpha)
    if (alpha == 1) {
      return(list(alpha = alpha, f = 2, df = Inf, d = Inf, par = alpha))
    }
    list(alpha = alpha, f = 0, df = 0, d = 0, par = alpha)
  }

  result <- new_schmidt_armijo_search(armijo_constant = 0.1, max_fn = 2)(
    phi = phi,
    step0 = step0,
    alpha = 1,
    pm = 1
  )

  expect_equal(evaluated, c(1, 0.25))
  expect_equal(result$step$alpha, 0.25)
  expect_equal(result$nfn, 2)
  expect_equal(result$ngr, 2)
})

test_that("More-Thuente successful steps satisfy strong Wolfe conditions", {
  cases <- list(
    list(
      name = "f1 small initial step",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1
    ),
    list(
      name = "f2 small initial step",
      fg = f2,
      x = 0,
      alpha = 1e-3,
      c1 = 0.1,
      c2 = 0.1
    ),
    list(
      name = "f5 large initial step",
      fg = f5,
      x = 0,
      alpha = 1e1,
      c1 = 0.001,
      c2 = 0.001
    ),
    list(
      name = "function modification",
      fg = f6,
      x = 1,
      alpha = 1,
      c1 = 0.1,
      c2 = 0.9
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    search <- new_wolfe_line_search(
      more_thuente_core,
      armijo_constant = case$c1,
      curvature_constant = case$c2,
      max_evaluations = 10000,
      method_policy = new_more_thuente_policy()
    )
    res <- search(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = case$alpha,
      pm = setup$pv
    )
    conditions <- new_line_condition_policy(case$c1, case$c2)

    expect_true(
      conditions$wolfe(setup$step0, res$step),
      info = case$name
    )
  }
})

test_that("Rasmussen and Schmidt successful steps satisfy strong Wolfe conditions", {
  cases <- list(
    list(
      name = "rasmussen",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1,
      run = function(setup, case) {
        new_wolfe_line_search(
          rasmussen_core,
          armijo_constant = case$c1,
          curvature_constant = case$c2,
          max_evaluations = 10000,
          method_policy = new_rasmussen_bracket_zoom_policy()
        )(
          phi = setup$phi,
          step0 = setup$step0,
          alpha = case$alpha,
          pm = setup$pv
        )
      }
    ),
    list(
      name = "schmidt",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1,
      run = function(setup, case) {
        new_wolfe_line_search(
          schmidt_core,
          armijo_constant = case$c1,
          curvature_constant = case$c2,
          max_evaluations = 10000,
          method_policy = new_schmidt_bracket_zoom_policy()
        )(
          phi = setup$phi,
          step0 = setup$step0,
          alpha = case$alpha,
          pm = setup$pv
        )
      }
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- case$run(setup, case)
    conditions <- new_line_condition_policy(case$c1, case$c2)

    expect_true(
      conditions$wolfe(setup$step0, res$step),
      info = case$name
    )
  }
})

test_that("weak Wolfe configuration does not require strong curvature", {
  c1 <- 1e-4
  c2 <- 0.1
  alpha <- 12.5
  setup <- condition_setup(condition_quadratic_fg(), x = 0, pv = 1)
  weak_conditions <- new_line_condition_policy(
    c1,
    c2,
    strong_curvature = FALSE
  )
  strong_conditions <- new_line_condition_policy(c1, c2)

  results <- list(
    `more-thuente` = new_wolfe_line_search(
      more_thuente_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 100,
      strong_curvature = FALSE,
      method_policy = new_more_thuente_policy()
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    ),
    rasmussen = new_wolfe_line_search(
      rasmussen_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 100,
      strong_curvature = FALSE,
      method_policy = new_rasmussen_bracket_zoom_policy()
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    ),
    schmidt = new_wolfe_line_search(
      schmidt_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 100,
      strong_curvature = FALSE,
      method_policy = new_schmidt_bracket_zoom_policy()
    )(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    )
  )

  for (name in names(results)) {
    step <- results[[name]]$step

    expect_equal(step$alpha, alpha, info = name)
    expect_true(weak_conditions$wolfe(setup$step0, step), info = name)
    expect_false(strong_conditions$wolfe(setup$step0, step), info = name)
  }
})

test_that("Hager-Zhang returned steps satisfy approximate weak Wolfe conditions", {
  cases <- list(
    list(
      name = "f1 small initial step",
      fg = f1,
      x = 0,
      alpha = 1e-3,
      c1 = 0.001,
      c2 = 0.1
    ),
    list(
      name = "f2 large initial step",
      fg = f2,
      x = 0,
      alpha = 1e3,
      c1 = 0.1,
      c2 = 0.1
    ),
    list(
      name = "f5 large initial step",
      fg = f5,
      x = 0,
      alpha = 1e1,
      c1 = 0.001,
      c2 = 0.001
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- line_search_hz(
      alpha = case$alpha,
      step0 = setup$step0,
      phi = setup$phi,
      c1 = case$c1,
      c2 = case$c2,
      max_fn = 100,
      strong_curvature = FALSE,
      always_check_convergence = TRUE,
      approx_armijo = TRUE
    )

    eps_k <- 1e-6 * abs(setup$step0$f)
    expect_true(
      hz_ok_step(
        step = res$step,
        step0 = setup$step0,
        c1 = case$c1,
        c2 = case$c2,
        eps = eps_k,
        strong_curvature = FALSE,
        approx_armijo = TRUE
      ),
      info = case$name
    )
  }
})

test_that("Hager-Zhang initial exhaustion distinguishes acceptance from fallback", {
  quartic_start <- list(alpha = 0, f = 1, df = 4, d = -16)
  quartic_trial <- list(alpha = 1, f = 81, df = -108, d = 432)

  expect_true(curvature_ok_step(quartic_start, quartic_trial, c2 = 0.5))
  expect_false(strong_curvature_ok_step(quartic_start, quartic_trial, c2 = 0.5))
  expect_false(armijo_ok_step(quartic_start, quartic_trial, c1 = 0.1))
  expect_false(approx_armijo_ok_step(quartic_start, quartic_trial, c1 = 0.1))
  expect_false(
    hz_ok_step(
      quartic_trial,
      quartic_start,
      c1 = 0.1,
      c2 = 0.5,
      eps = 1e-6,
      strong_curvature = FALSE,
      approx_armijo = TRUE
    )
  )

  linear_start <- list(alpha = 0, f = 1, df = -1, d = -1)
  linear_trial <- list(alpha = 1, f = 0, df = -1, d = -1)

  expect_true(armijo_ok_step(linear_start, linear_trial, c1 = 0.1))
  expect_false(curvature_ok_step(linear_start, linear_trial, c2 = 0.5))
  expect_false(
    hz_ok_step(
      linear_trial,
      linear_start,
      c1 = 0.1,
      c2 = 0.5,
      eps = 1e-6,
      strong_curvature = FALSE,
      approx_armijo = TRUE
    )
  )
})

test_that("Hager-Zhang U3 bracket update bisects to a positive-slope bound", {
  hz_step <- function(alpha, f, d) {
    list(alpha = alpha, f = f, d = d, df = d)
  }

  step0 <- hz_step(alpha = 0, f = 0, d = -1)
  bracket <- list(
    hz_step(alpha = 1, f = -1, d = -1),
    hz_step(alpha = 8, f = 2, d = 1)
  )
  step_c <- hz_step(alpha = 8, f = 2, d = -1)
  phi <- function(alpha) {
    if (alpha > 4) {
      hz_step(alpha, f = 1, d = -1)
    } else if (alpha < 3) {
      hz_step(alpha, f = -1, d = -1)
    } else {
      hz_step(alpha, f = -0.5, d = 1)
    }
  }

  res <- update_bracket_hz(
    bracket = bracket,
    step_c = step_c,
    step0 = step0,
    phi = phi,
    eps = 0,
    max_fn = 5,
    theta = 0.5
  )

  expect_true(res$ok)
  expect_equal(res$nfn, 3)
  expect_equal(unname(bracket_props(res$bracket, "alpha")), c(2.75, 3.625))
  expect_equal(unname(bracket_props(res$bracket, "d")), c(-1, 1))
  expect_true(res$bracket[[1]]$f <= step0$f)
})

test_that("line searches return the initial step when evaluation budgets are exhausted", {
  c1 <- 1e-4
  c2 <- 0.1
  alpha <- 12.5
  setup <- condition_setup(condition_quadratic_fg(), x = 0, pv = 1)
  weak_conditions <- new_line_condition_policy(
    c1,
    c2,
    strong_curvature = FALSE
  )

  wolfe_searches <- list(
    `more-thuente` = new_wolfe_line_search(
      more_thuente_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      method_policy = new_more_thuente_policy(
        alpha_max = Inf,
        safeguard_cubic = FALSE
      )
    ),
    rasmussen = new_wolfe_line_search(
      rasmussen_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      method_policy = new_rasmussen_bracket_zoom_policy(
        interior_fraction = 0.1,
        expansion_factor = 3,
        relative_interval_tolerance = 1e-6
      )
    ),
    schmidt = new_wolfe_line_search(
      schmidt_core,
      armijo_constant = c1,
      curvature_constant = c2,
      max_evaluations = 0,
      method_policy = new_schmidt_bracket_zoom_policy()
    ),
    `hager-zhang` = hager_zhang(c1 = c1, c2 = c2, max_fn = 0)
  )

  for (name in names(wolfe_searches)) {
    res <- wolfe_searches[[name]](
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      pm = setup$pv
    )

    expect_equal(res$nfn, 0, info = name)
    expect_equal(res$ngr, 0, info = name)
    expect_equal(res$step$alpha, 0, info = name)
    expect_false(weak_conditions$wolfe(setup$step0, res$step), info = name)
  }

  armijo_res <- new_schmidt_armijo_search(
    armijo_constant = c1,
    max_fn = 0
  )(
    phi = setup$phi,
    step0 = setup$step0,
    alpha = alpha,
    pm = setup$pv
  )

  expect_equal(armijo_res$nfn, 0)
  expect_equal(armijo_res$ngr, 0)
  expect_equal(armijo_res$step$alpha, 0)
})

test_that("finite-value guard backs off non-finite line-search evaluations", {
  fg <- list(
    fn = function(x) {
      if (x > 1) {
        return(Inf)
      }
      (x - 0.25)^2
    },
    gr = function(x) {
      if (x > 1) {
        return(Inf)
      }
      2 * (x - 0.25)
    }
  )
  setup <- condition_setup(fg, x = 0, pv = 1)

  res <- find_finite(setup$phi, alpha = 4, min_alpha = 0, max_fn = 4)

  expect_true(res$ok)
  expect_lte(res$step$alpha, 1)
  expect_true(step_is_finite(res$step))
})
