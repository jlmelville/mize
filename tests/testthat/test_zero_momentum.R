zero_momentum_result <- function(method_args, max_iter, max_fg = Inf) {
  do.call(
    mize,
    c(
      list(
        par = rb0,
        fg = rosenbrock_fg,
        line_search = "More-Thuente",
        c2 = 0.1,
        max_iter = max_iter,
        max_fg = max_fg,
        check_conv_every = 1,
        abs_tol = NULL,
        rel_tol = NULL,
        grad_tol = NULL,
        ginf_tol = NULL,
        step_tol = NULL,
        store_progress = TRUE
      ),
      method_args
    )
  )
}

zero_momentum_methods <- list(
  Momentum = list(
    method = "Momentum",
    mom_type = "classical",
    mom_schedule = 0
  ),
  NAG = list(
    method = "NAG",
    nest_q = 1,
    nest_convex_approx = FALSE
  )
)

test_that("static zero momentum follows the SD execution path", {
  result_fields <- c(
    "par",
    "f",
    "nf",
    "ng",
    "iter",
    "terminate",
    "converged",
    "status",
    "message",
    "best_par",
    "best_f",
    "last_par",
    "last_f",
    "progress"
  )

  for (run_control in list(
    iteration_limit = list(max_iter = 3, max_fg = Inf),
    evaluation_limit = list(max_iter = 500, max_fg = 30)
  )) {
    sd <- zero_momentum_result(
      list(method = "SD"),
      max_iter = run_control$max_iter,
      max_fg = run_control$max_fg
    )

    for (method_name in names(zero_momentum_methods)) {
      result <- zero_momentum_result(
        zero_momentum_methods[[method_name]],
        max_iter = run_control$max_iter,
        max_fg = run_control$max_fg
      )

      expect_identical(
        result[result_fields],
        sd[result_fields],
        info = method_name
      )
      expect_false("mu" %in% names(result$progress), info = method_name)
    }
  }
})

test_that("static zero momentum constructs the SD stage topology", {
  omitted_stage_args <- c(
    zero_momentum_methods,
    list(
      Momentum_with_restart = list(
        method = "Momentum",
        mom_type = "classical",
        mom_schedule = 0,
        mom_linear_weight = TRUE,
        restart = "fn"
      ),
      NAG_explicit_schedule = list(
        method = "NAG",
        mom_schedule = "nsconvex",
        nest_q = 1,
        nest_convex_approx = FALSE,
        mom_linear_weight = TRUE,
        restart = "fn"
      )
    )
  )

  for (method_name in names(omitted_stage_args)) {
    opt <- do.call(make_mize, omitted_stage_args[[method_name]])
    expect_named(opt$stages, "gradient_descent", info = method_name)
  }
})

test_that("dynamic and nonzero momentum stages remain configured", {
  retained_stage_args <- list(
    nonzero_constant = list(
      method = "Momentum",
      mom_schedule = 0.5
    ),
    nesterov_zero = list(
      method = "Momentum",
      mom_type = "nesterov",
      mom_schedule = 0
    ),
    scheduled_zero = list(
      method = "Momentum",
      mom_schedule = "ramp",
      mom_init = 0,
      mom_final = 0
    ),
    function_zero = list(
      method = "Momentum",
      mom_schedule = function(iter) 0
    ),
    exact_nonzero_nag = list(
      method = "NAG",
      nest_q = 0.5,
      nest_convex_approx = FALSE
    ),
    approximate_nag = list(
      method = "NAG",
      nest_q = 1,
      nest_convex_approx = TRUE
    ),
    explicit_numeric_nag = list(
      method = "NAG",
      mom_schedule = 0,
      nest_q = 1,
      nest_convex_approx = FALSE
    )
  )

  for (method_name in names(retained_stage_args)) {
    opt <- do.call(make_mize, retained_stage_args[[method_name]])
    expect_true("momentum" %in% names(opt$stages), info = method_name)
  }
})

test_that("stateful static zero momentum transitions match SD", {
  make_stateful_opt <- function(method_args) {
    do.call(
      make_mize,
      c(
        list(
          line_search = "More-Thuente",
          c2 = 0.1,
          par = rb0,
          fg = rosenbrock_fg,
          max_iter = 3,
          max_fg = Inf,
          abs_tol = NULL,
          rel_tol = NULL,
          grad_tol = NULL,
          ginf_tol = NULL,
          step_tol = NULL
        ),
        method_args
      )
    )
  }

  for (method_name in names(zero_momentum_methods)) {
    sd_opt <- make_stateful_opt(list(method = "SD"))
    zero_opt <- make_stateful_opt(zero_momentum_methods[[method_name]])
    sd_par <- rb0
    zero_par <- rb0

    for (iter in seq_len(3)) {
      sd_step <- mize_step(sd_opt, sd_par, rosenbrock_fg)
      zero_step <- mize_step(zero_opt, zero_par, rosenbrock_fg)

      expect_identical(zero_step$par, sd_step$par, info = method_name)
      expect_identical(
        zero_step[c("nf", "ng")],
        sd_step[c("nf", "ng")],
        info = method_name
      )

      sd_opt <- sd_step$opt
      zero_opt <- zero_step$opt
      sd_par <- sd_step$par
      zero_par <- zero_step$par
    }
  }
})
