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

test_that("Armijo backtracking accepts steps with sufficient decrease", {
  cases <- list(
    list(name = "f1 small initial step", fg = f1, x = 0, alpha = 1e-1, c1 = 0.001),
    list(name = "f2 large initial step", fg = f2, x = 0, alpha = 1e1, c1 = 0.1),
    list(name = "f4 very large initial step", fg = f4, x = 0, alpha = 1e3, c1 = 0.001)
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- ArmijoBacktrack(
      step = case$alpha,
      f = setup$step0$f,
      g = setup$step0$df,
      gtd = setup$step0$d,
      c1 = case$c1,
      LS_interp = 2,
      LS_multi = 0,
      maxLS = 10000,
      funObj = setup$phi,
      pnorm_inf = max(abs(setup$pv)),
      progTol = 1e-9,
      debug = FALSE
    )

    expect_true(res$step$alpha > 0, info = case$name)
    expect_true(
      armijo_ok_step(setup$step0, res$step, case$c1),
      info = case$name
    )
  }
})

test_that("More-Thuente successful steps satisfy strong Wolfe conditions", {
  cases <- list(
    list(name = "f1 small initial step", fg = f1, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.1),
    list(name = "f2 small initial step", fg = f2, x = 0, alpha = 1e-3, c1 = 0.1, c2 = 0.1),
    list(name = "f5 large initial step", fg = f5, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001),
    list(name = "function modification", fg = f6, x = 1, alpha = 1, c1 = 0.1, c2 = 0.9)
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- cvsrch(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = case$alpha,
      c1 = case$c1,
      c2 = case$c2
    )

    expect_equal(res$info, 1, info = case$name)
    expect_true(
      strong_wolfe_ok_step(setup$step0, res$step, case$c1, case$c2),
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
        ras_ls(
          phi = setup$phi,
          alpha = case$alpha,
          step0 = setup$step0,
          c1 = case$c1,
          c2 = case$c2,
          max_fn = 10000,
          verbose = FALSE
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
        WolfeLineSearch(
          alpha = case$alpha,
          f = setup$step0$f,
          g = setup$step0$df,
          gtd = setup$step0$d,
          c1 = case$c1,
          c2 = case$c2,
          LS_interp = 2,
          LS_multi = 0,
          maxLS = 10000,
          funObj = setup$phi,
          pnorm_inf = max(abs(setup$pv)),
          progTol = 1e-9,
          debug = FALSE
        )
      }
    )
  )

  for (case in cases) {
    setup <- condition_setup(case$fg, case$x)
    res <- case$run(setup, case)

    expect_true(
      strong_wolfe_ok_step(setup$step0, res$step, case$c1, case$c2),
      info = case$name
    )
  }
})

test_that("weak Wolfe configuration does not require strong curvature", {
  c1 <- 1e-4
  c2 <- 0.1
  alpha <- 12.5
  setup <- condition_setup(condition_quadratic_fg(), x = 0, pv = 1)
  weak_ok_step <- make_wolfe_ok_step_fn(strong_curvature = FALSE)
  strong_ok_step <- make_wolfe_ok_step_fn(strong_curvature = TRUE)

  results <- list(
    `more-thuente` = cvsrch(
      phi = setup$phi,
      step0 = setup$step0,
      alpha = alpha,
      c1 = c1,
      c2 = c2,
      wolfe_ok_step_fn = weak_ok_step
    ),
    rasmussen = ras_ls(
      phi = setup$phi,
      alpha = alpha,
      step0 = setup$step0,
      c1 = c1,
      c2 = c2,
      max_fn = 100,
      wolfe_ok_step_fn = weak_ok_step,
      verbose = FALSE
    ),
    schmidt = WolfeLineSearch(
      alpha = alpha,
      f = setup$step0$f,
      g = setup$step0$df,
      gtd = setup$step0$d,
      c1 = c1,
      c2 = c2,
      LS_interp = 2,
      LS_multi = 0,
      maxLS = 100,
      funObj = setup$phi,
      pnorm_inf = max(abs(setup$pv)),
      progTol = 1e-9,
      curvature_check_fn = curvature_ok_step,
      debug = FALSE
    )
  )

  for (name in names(results)) {
    step <- results[[name]]$step

    expect_equal(step$alpha, alpha, info = name)
    expect_true(weak_ok_step(setup$step0, step, c1, c2), info = name)
    expect_false(strong_ok_step(setup$step0, step, c1, c2), info = name)
  }
})

test_that("Hager-Zhang returned steps satisfy approximate weak Wolfe conditions", {
  cases <- list(
    list(name = "f1 small initial step", fg = f1, x = 0, alpha = 1e-3, c1 = 0.001, c2 = 0.1),
    list(name = "f2 large initial step", fg = f2, x = 0, alpha = 1e3, c1 = 0.1, c2 = 0.1),
    list(name = "f5 large initial step", fg = f5, x = 0, alpha = 1e1, c1 = 0.001, c2 = 0.001)
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
  weak_ok_step <- make_wolfe_ok_step_fn(strong_curvature = FALSE)

  wolfe_searches <- list(
    `more-thuente` = more_thuente(c1 = c1, c2 = c2, max_fn = 0),
    rasmussen = rasmussen(c1 = c1, c2 = c2, max_fn = 0),
    schmidt = schmidt(c1 = c1, c2 = c2, max_fn = 0),
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
    expect_false(weak_ok_step(setup$step0, res$step, c1, c2), info = name)
  }

  armijo_res <- schmidt_armijo_backtrack(c1 = c1, max_fn = 0)(
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
