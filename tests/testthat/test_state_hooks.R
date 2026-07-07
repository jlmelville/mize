state_hook_quadratic <- function() {
  list(
    fn = function(x) {
      sum(x^2)
    },
    gr = function(x) {
      2 * x
    }
  )
}

register_test_hooks <- function(opt, hooks) {
  for (hook in hooks) {
    opt <- register_hook(opt, hook)
  }
  opt
}

make_trace_opt_hook <- function(name, event) {
  hook <- function(opt, par, fg, iter, ...) {
    opt$trace <- c(opt$trace, name)
    opt
  }
  attr(hook, "name") <- name
  attr(hook, "event") <- event
  hook
}

make_trace_stage_hook <- function(name, event) {
  hook <- function(opt, stage, par, fg, iter, ...) {
    opt$trace <- c(opt$trace, name)
    list(opt = opt)
  }
  attr(hook, "name") <- name
  attr(hook, "event") <- event
  hook
}

test_that("stateful steps run lifecycle hooks in optimizer order", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )
  opt <- register_test_hooks(
    opt,
    list(
      make_trace_opt_hook("before_step", "before step"),
      make_trace_stage_hook(
        "before_gradient_descent",
        "before gradient_descent"
      ),
      make_trace_stage_hook("after_stage", "after stage"),
      make_trace_opt_hook("before_validation", "before validation"),
      make_trace_opt_hook("during_validation", "during validation"),
      make_trace_opt_hook("after_step", "after step")
    )
  )

  res <- mize_step(opt, c(1), fg)

  expect_equal(
    res$opt$trace,
    c(
      "before_step",
      "before_gradient_descent",
      "after_stage",
      "before_validation",
      "during_validation",
      "after_step"
    )
  )
  expect_equal(res$par, 0.8)
})

test_that("stage and sub-stage hooks persist returned state mutations", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  mark_stage <- function(opt, stage, par, fg, iter, ...) {
    stage$marked_by_hook <- iter
    list(stage = stage)
  }
  attr(mark_stage, "name") <- "mark_stage"
  attr(mark_stage, "event") <- "before gradient_descent"

  double_direction <- function(opt, stage, sub_stage, par, fg, iter, ...) {
    sub_stage$value <- 2 * sub_stage$value
    sub_stage$marked_by_hook <- iter
    list(sub_stage = sub_stage)
  }
  attr(double_direction, "name") <- "double_direction"
  attr(double_direction, "event") <- "after gradient_descent direction"

  opt <- register_test_hooks(opt, list(mark_stage, double_direction))

  res <- mize_step(opt, c(1), fg)
  stage <- res$opt$stages$gradient_descent

  expect_equal(res$par, 0.6)
  expect_equal(stage$marked_by_hook, 1)
  expect_equal(stage$direction$marked_by_hook, 1)
  expect_equal(stage$direction$value, -4)
  expect_equal(stage$result, -0.4)
})

test_that("termination from a stage hook short-circuits later lifecycle hooks", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  terminate_after_stage <- function(opt, stage, par, fg, iter, ...) {
    opt$trace <- c(opt$trace, "after_stage")
    opt <- set_mize_termination(
      opt,
      list(what = "hook_stop", val = iter)
    )
    list(opt = opt)
  }
  attr(terminate_after_stage, "name") <- "terminate_after_stage"
  attr(terminate_after_stage, "event") <- "after stage"

  opt <- register_test_hooks(
    opt,
    list(
      terminate_after_stage,
      make_trace_opt_hook("before_validation", "before validation"),
      make_trace_opt_hook("after_step", "after step")
    )
  )

  res <- mize_step(opt, c(1), fg)

  expect_equal(res$par, 1)
  expect_equal(res$opt$trace, "after_stage")
  expect_equal(res$opt$terminate$what, "hook_stop")
  expect_equal(res$opt$status, "terminated")
})

test_that("validation veto rolls back parameters before after-step hooks", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  veto_validation <- function(opt, par, fg, iter, par0, update) {
    opt$trace <- c(opt$trace, "veto")
    opt$ok <- FALSE
    opt
  }
  attr(veto_validation, "name") <- "veto_validation"
  attr(veto_validation, "event") <- "during validation"

  after_step_probe <- function(opt, par, fg, iter, par0, update) {
    opt$trace <- c(
      opt$trace,
      "after_step",
      paste0("after_par=", par[1]),
      paste0("update=", update[1]),
      paste0("ok=", opt$ok)
    )
    opt
  }
  attr(after_step_probe, "name") <- "after_step_probe"
  attr(after_step_probe, "event") <- "after step"

  opt <- register_test_hooks(opt, list(veto_validation, after_step_probe))

  res <- mize_step(opt, c(1), fg)

  expect_equal(res$par, 1)
  expect_equal(
    res$opt$trace,
    c("veto", "after_step", "after_par=1", "update=-0.2", "ok=FALSE")
  )
  expect_false(res$opt$ok)
})

test_that("speed adaptive restart rolls back and resets momentum state", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "MOM",
    line_search = "constant",
    step0 = 0.1,
    mom_schedule = 0.5,
    restart = "speed",
    restart_wait = 1,
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  par <- c(1)
  first <- mize_step(opt, par, fg)
  second <- mize_step(first$opt, first$par, fg)
  third <- mize_step(second$opt, second$par, fg)
  fourth <- mize_step(third$opt, third$par, fg)

  expect_equal(first$par, 0.8)
  expect_equal(first$opt$cache$update_old, -0.2)
  expect_true(first$opt$ok)

  expect_equal(second$par, 0.54)
  expect_equal(second$opt$cache$update_old, -0.26)
  expect_true(second$opt$ok)

  expect_equal(third$par, second$par)
  expect_false(third$opt$ok)
  expect_equal(third$opt$restart_at, 3)
  expect_equal(third$opt$cache$update_old, 0)

  expect_equal(fourth$par, 0.432)
  expect_true(fourth$opt$ok)
  expect_equal(fourth$opt$restart_at, 3)
  expect_equal(fourth$opt$cache$update_old, -0.108)
})

test_that("adaptive restart replaces update_old hook behavior", {
  fg <- state_hook_quadratic()
  opt <- make_mize(
    method = "MOM",
    line_search = "constant",
    step0 = 0.1,
    mom_schedule = 0.5,
    restart = "speed",
    par = c(1),
    fg = fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  after_step_hooks <- names(opt$hooks[["step"]][["after"]])
  expect_equal(sum(after_step_hooks == "update_old"), 1)
  expect_true(
    "validate_speed" %in% names(opt$hooks[["validation"]][["during"]])
  )

  failed <- opt
  failed$ok <- FALSE
  failed$cache$update_old <- -99
  failed <- life_cycle_hook(
    "step",
    "after",
    failed,
    par = c(1),
    fg = fg,
    iter = 1,
    par0 = c(1),
    update = c(-0.2)
  )

  expect_equal(failed$restart_at, 1)
  expect_equal(failed$cache$update_old, 0)

  accepted <- opt
  accepted$ok <- TRUE
  accepted$cache$update_old <- -99
  accepted <- life_cycle_hook(
    "step",
    "after",
    accepted,
    par = c(1),
    fg = fg,
    iter = 1,
    par0 = c(1),
    update = c(-0.2)
  )

  expect_null(accepted$restart_at)
  expect_equal(accepted$cache$update_old, -0.2)
})

test_that("eager updates expose per-stage parameters and summed validation update", {
  fg <- state_hook_quadratic()

  make_fixed_direction <- function(value) {
    make_direction(list(
      calculate = function(opt, stage, sub_stage, par, fg, iter) {
        opt$trace <- c(
          opt$trace,
          paste0(stage$type, "_during_par=", par[1])
        )
        sub_stage$value <- value
        list(opt = opt, sub_stage = sub_stage)
      }
    ))
  }

  opt <- make_opt(
    make_stages(
      make_stage(
        "first_stage",
        direction = make_fixed_direction(2),
        step_size = constant_step_size(value = 1)
      ),
      make_stage(
        "second_stage",
        direction = make_fixed_direction(3),
        step_size = constant_step_size(value = 1)
      )
    )
  )
  opt$eager_update <- TRUE
  opt <- mize_init(
    opt,
    c(10),
    fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  after_stage_probe <- function(opt, stage, par, fg, iter, ...) {
    opt$trace <- c(
      opt$trace,
      paste0(stage$type, "_after_par=", par[1])
    )
    list(opt = opt)
  }
  attr(after_stage_probe, "name") <- "after_stage_probe"
  attr(after_stage_probe, "event") <- "after stage"

  validation_probe <- function(opt, par, fg, iter, par0, update) {
    opt$trace <- c(
      opt$trace,
      paste0("validation_par=", par[1]),
      paste0("validation_par0=", par0[1]),
      paste0("validation_update=", update[1])
    )
    opt
  }
  attr(validation_probe, "name") <- "validation_probe"
  attr(validation_probe, "event") <- "during validation"

  opt <- register_hook(opt, after_stage_probe)
  opt <- register_hook(opt, validation_probe)

  res <- mize_step(opt, c(10), fg)

  expect_equal(res$par, 15)
  expect_equal(
    res$opt$trace,
    c(
      "first_stage_during_par=10",
      "first_stage_after_par=12",
      "second_stage_during_par=12",
      "second_stage_after_par=15",
      "validation_par=15",
      "validation_par0=10",
      "validation_update=5"
    )
  )
})
