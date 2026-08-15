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
  expect_no_error(make_mize(mom_schedule = "switch", mom_switch_iter = 0))
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
